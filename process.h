#pragma once
#ifndef PROCESS_H
#define PROCESS_H
#include <iostream>
#include <fstream>
#include "parameters.h"
#include "calculate.h"

#define M_PI       3.14159265358979323846   // pi
using namespace std;

enum { Conc_mm3, VolFrac };
enum { MonoDisperse, PolyDisperse };

enum { Log_Normal, Guassian, Junge, CustomData };



//Initialize dynamic arrays.
void InitializeArrays(parameters* para, int p1)
{
	para->minTheta = 0;
	para->maxTheta = M_PI;
	if (p1 == Phase_DTheta0_1)
		para->nTheta = 1801;    //180/(1801-1) = 0.1degree step
	if (p1 == Phase_DTheta0_5)
		para->nTheta = 361;    //180/(361-1) = 0.5degree step
	if (p1 == Phase_DTheta_1)
		para->nTheta = 181;    //180/(181-1) = 1degree step

	para->stepTheta = (para->maxTheta - para->minTheta) / (double)(para->nTheta - 1);
	para->nWavel = (floor(para->endWavel - para->startWavel) / para->stepWavel) + 1;
	para->wavelArray = new double[para->nWavel];
	for (int i = 0; i < para->nWavel; ++i)
		para->wavelArray[i] = para->startWavel + i * para->stepWavel;

	para->phaseFunctionAve = new double* [para->nWavel];
	for (int i = 0; i < para->nWavel; i++)
		para->phaseFunctionAve[i] = new double[para->nTheta];

	para->phaseFunctionPara = new double* [para->nWavel];
	for (int i = 0; i < para->nWavel; i++)
		para->phaseFunctionPara[i] = new double[para->nTheta];

	para->phaseFunctionPerp = new double* [para->nWavel];
	for (int i = 0; i < para->nWavel; i++)
		para->phaseFunctionPerp[i] = new double[para->nTheta];

	para->S1 = new std::complex<double> * [para->nWavel];
	for (int i = 0; i < para->nWavel; i++)
		para->S1[i] = new std::complex<double>[para->nTheta];

	para->S2 = new std::complex<double> * [para->nWavel];
	for (int i = 0; i < para->nWavel; i++)
		para->S2[i] = new std::complex<double>[para->nTheta];

	para->cSca = new double[para->nWavel];
	para->cExt = new double[para->nWavel];
	para->cBack = new double[para->nWavel];
	para->SizePara = new double[para->nWavel];
	para->mus = new double[para->nWavel];
	para->g = new double[para->nWavel];
	para->forward = new double[para->nWavel];
	para->backward = new double[para->nWavel];

	//*arrayFlag = true;
}

//Delete dynamic arrays
void DeleteArrays(parameters* para)
{
	delete[] para->wavelArray;
	delete[] para->cSca;
	delete[] para->cExt;
	delete[] para->cBack;
	delete[] para->SizePara;
	delete[] para->mus;
	delete[] para->g;
	delete[] para->forward;
	delete[] para->backward;

	for (int i = 0; i < para->nWavel; i++)
		delete[] para->phaseFunctionAve[i];
	for (int i = 0; i < para->nWavel; i++)
		delete[] para->phaseFunctionPara[i];
	for (int i = 0; i < para->nWavel; i++)
		delete[] para->phaseFunctionPerp[i];
	for (int i = 0; i < para->nWavel; i++)
		delete[] para->S1[i];
	for (int i = 0; i < para->nWavel; i++)
		delete[] para->S2[i];

	//*arrayFlag = false;
}

// Run Mono disperse distribution
void ProcessMonoDisperse(parameters* para, calculate* mCalc)
{
	double tempMus1000;   //to hold the value of scattering cross section at 1000nm
	double tempG1000;           //to hold the value of G at 1000nm
//    double numberDensity = para->sphNumDensity;   //Concentration in 1mm^3


	//Run Mie simulation
	mCalc->DoSimulation(para, &tempMus1000, &tempG1000);

	//Calculate fitting constant 'A'
	para->fittedA = tempMus1000 * (1 - tempG1000);

}

// Run Poly disperse distribution
void ProcessPolyDisperse(parameters* para, calculate* mCalc)
{
	double tempMus1000;   //to hold the value of scattering cross section at 1000nm
	double tempG1000;           //to hold the value of G at 1000nm


	//Run Mie simulation for radius para->radArray[i]
	mCalc->DoSimulation(para, &tempMus1000, &tempG1000);

	//Calculate fitting constant 'A'
	para->fittedA = tempMus1000 * (1 - tempG1000);

}

//Sphere distribution in poly disperse
void ProcessDistribution(parameters* para, int distIndex, int p, calculate* mCalc)
{
	if (distIndex != CustomData)
	{
		bool flagVolOrConc = true;
		para->radArray = new double[para->nRadius];
		para->numDensityArray = new double[para->nRadius];
		para->scatRefRealArray = new double[para->nRadius];
		para->scatRefImagArray = new double[para->nRadius];

		//Find size of spheres
		if (p == VolFrac)
			flagVolOrConc = true;
		if (p == Conc_mm3)
			flagVolOrConc = false;
		mCalc->DiameterRangeSetting(para, distIndex);
		mCalc->SetSphereRadiusAndRefIndex(para, distIndex, flagVolOrConc);
	}
}

int CountLines1(const char* filename)
{
	ifstream ReadFile;
	int n = 0;
	char line[512];
	string temp;
	ReadFile.open(filename, ios::in);//ios::in 表示以只读的方式读取文件
	if (ReadFile.fail())//文件打开失败:返回0
	{
		return 0;
	}
	else//文件存在
	{
		while (getline(ReadFile, temp))
		{
			n++;
		}
		return n;
	}
	ReadFile.close();
}

void ReadCustomData(parameters* para, const char* fileName)
{
	FILE* fp = fopen(fileName, "r");
	if (fp == NULL)
	{
		cerr << "cannot open  file" << endl;
	}
	int count = CountLines1(fileName);
	double rad, numDensity, scatRefReal, scatRefImag;

	if (count < 1)
	{
		cout << "Error";
		return;
	}
	else
	{
		para->nRadius = count;
		para->radArray = new double[para->nRadius];
		para->numDensityArray = new double[para->nRadius];
		para->scatRefRealArray = new double[para->nRadius];
		para->scatRefImagArray = new double[para->nRadius];

		int idx = 0;
		double sumRad = 0;
		double minRad = 1e100;
		double maxRad = 0;
		for (int idx = 0; idx < count; idx++)
		{
			fscanf(fp, "%lf\t%lf\t%lf\t%lf", &rad, &numDensity, &scatRefReal, &scatRefImag);
			para->radArray[idx] = rad / 2.0;
			para->numDensityArray[idx] = numDensity;
			para->scatRefRealArray[idx] = scatRefReal;
			para->scatRefImagArray[idx] = scatRefImag;
			sumRad += para->radArray[idx];
			if (para->radArray[idx] > maxRad)
				maxRad = para->radArray[idx];
			if (para->radArray[idx] < minRad)
				minRad = para->radArray[idx];
		}

		fclose(fp);
		para->meanRadius = sumRad / para->nRadius;
		para->minRadius = minRad;
		para->maxRadius = maxRad;
		//Assign this values to pass the sanity check
		para->stdDev = 1.0;
	}

}








#endif 