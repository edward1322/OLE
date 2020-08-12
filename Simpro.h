#pragma once
#ifndef MCLIB_H
#define MCLIB_H
#include "utilities.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <random>
#include <iomanip>
#include "Vector3.h"
#include "light.h"
#include "roughsurface.h"
#include "Interface.h"
#include "transfer.h"
#include "calculate.h"
#include "process.h"
void InitOutputData(parameters* para)
{
	int length = para->BssArrayLength;
	para->bss = new double[length];
	para->bss1 = new double[length];
	para->bss2 = new double[length];
	para->bss3 = new double[length];
	para->bss3plus = new double[length];
	for (int i = 0; i < length; i++)
	{
		para->bss[i] = 0;
		para->bss1[i] = 0;
		para->bss2[i] = 0;
		para->bss3[i] = 0;
		para->bss3plus[i] = 0;
	}
}
void LoadseaData(parameters* para)
{
	//read water properties <bio or opt file>
	FILE* fp = NULL;
	int length = para->BssArrayLength;
	para->mua0 = new double[length];
	para->mus0 = new double[length];
	para->muat = new double[length];
	double a, s;
	fp = fopen(para->custom_in, "r");
	if (fp == NULL)
	{
		printf("cannot open opt.txt for water optical file\n");
	}
	for (int i = 0; i < length; i++)
	{
		fscanf(fp, "%lf%lf", &a, &s);

		para->mua0[i] = a + 0.05070;
		para->mus0[i] = s + 0.00170;
		/*	Out->mua[i] = a ;
			Out->mus[i] = s ;*/
		para->muat[i] = a + s;

	}
	fclose(fp);
}
void LoadseabottomData(parameters* para)
{
	if (para->fseabotttype != ruser_def)
	{
		FILE* fp = NULL;
		int length = 91;
		int wav=300;
		double r1=0, r2=0, r3=0, r4=0;
		fp = fopen(para->seabotfile, "r");
		if (fp == NULL)
		{
			printf("cannot open SEA_BOTTOM_REFLECTANCES.txt for seabott file\n");
		}
		for (int i = 0; i < length; i++)
		{
			fscanf(fp, "%d%lf%lf%lf%lf", &wav, &r1, &r2, &r3, &r4);
			if (para->laserWavel <= wav)
			{
				switch (para->fseabotttype)
				{
				case 1:
					para->fseabott = r1;
					break;
				case 2:
					para->fseabott = r2;
					break; 
				case 3:
					para->fseabott = r3;
						break;
				case 4:
					para->fseabott = r4;
					break;
				default:
					break;
				}
			}
		}
		fclose(fp);
	}	
}

int CountLines(const char* filename)
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
void LoadCDFData(parameters* para)
{
	FILE* fp_cdf = NULL;
	int len;

	fp_cdf = fopen(para->custom_cdf, "r");
	if (fp_cdf == NULL)
	{
		cerr << "cannot open cdf.inp" << endl;
	}
	len = CountLines(para->custom_cdf);

	para->i_len = len;
	para->i_m_rnd = new double[len];
	para->i_angle = new double[len];
	double rnd, ang;
	for (int i = 0; i < len; i++)
	{
		fscanf(fp_cdf, "%lf%lf", &rnd, &ang);

		para->i_m_rnd[i] = rnd;
		para->i_angle[i] = ang;
	}
	fclose(fp_cdf);
}
//void GenerateCDFData(parameters* para)
//{
//	para->imie_m_rnd = new double[para->nTheta];
//	para->imie_angle = new double[para->nTheta];
//	double tmp = (para->S1[0][0].real() * para->S1[0][0].real() + para->S2[0][0].real() * para->S2[0][0].real()) / (para->cSca[0] / (M_PI * para->meanRadius * para->meanRadius) * para->SizePara[0] * para->SizePara[0]);
//	for (int t = 0; t < para->nTheta; t++)
//		{		    
//		    para->imie_angle[t] = para->minTheta + t * para->stepTheta;
//			para->imie_m_rnd[t] = (para->S1[0][t].real() * para->S1[0][t].real() + para->S2[0][t].real() * para->S2[0][t].real()) / (para->cSca[0]/ (M_PI*para->meanRadius * para->meanRadius)*para->SizePara[0]*para->SizePara[0])/tmp;
//			//printf("%lf\n", para->imie_m_rnd[t]);
//		}
//}
void LoadSPFData(parameters* para)
{
	FILE* fp_spf = NULL;
	int len;

	fp_spf = fopen(para->custom_spf, "r");
	if (fp_spf == NULL)
	{
		cerr << "cannot open spf.inp" << endl;
	}
	len = CountLines(para->custom_spf);

	para->s_len = len;
	para->s_m_rnd = new double[len];
	para->s_angle = new double[len];
	double rnd, ang;
	for (int i = 0; i < len; i++)
	{
		fscanf(fp_spf, "%lf%lf", &rnd, &ang);

		para->s_m_rnd[i] = rnd;
		para->s_angle[i] = ang;
	}
	fclose(fp_spf);
}

void Initpara(parameters* para)
{
	RandomGen(0, -(int)time(NULL) % (1 << 15), NULL);
	// photon parameters
	para->x = 0;
	para->y = 0;
	para->z = 0;
	para->ux = 0;
	para->uy = 0;
	para->uz = 1;
	para->status = ALIVE;
	para->hit = NOHIT;
	para->totalpath = 0;
	para->sleft = 0;
	para->weight = 1.0;
	para->collision = 0;
	para->steplength = 0.0;
	InitOutputData(para);
	LoadseaData(para);
	LoadseabottomData(para);
	ocean_surface_scalar(para->windspeed, para->inc_ang, para->interface);
	if (para->it== I_CUSTOME) LoadCDFData(para);
	if (para->st == S_CUSTOME) LoadSPFData(para);

	if ((para->it == I_MIE) || (para->st == MIE_CAL))
	{
		calculate* mCalc = new calculate();
		InitializeArrays(para, para->DTheta);  //Phase_DTheta0_1, Phase_DTheta0_5
		if (para->Disperse == MonoDisperse) para->distIndex = -1;

		if (para->distIndex == CustomData)
		{
			ReadCustomData(para, para->custom_psd);
		}
		else
			ProcessDistribution(para, para->distIndex, para->m_type, mCalc);

		if (para->Disperse == MonoDisperse)
		{
			ProcessMonoDisperse(para, mCalc);
		}
		else
		{
			ProcessPolyDisperse(para, mCalc);
		}
		//GenerateCDFData(para);

		//for (int t = 0; t < para->nTheta; t++)
		//{

		//	//printf("%5.5f\t %5.5f\t\n", para->S1[0][t].real(), para->S2[0][t].real());

		//	printf("%5.5f\t %5.5f\t %5.5f\t\n", para->phaseFunctionAve[0][t], para->phaseFunctionPara[0][t], para->phaseFunctionPerp[0][t]);

		//}
		//std::cout << "Hello World!\n";	
	}
		
}
//Delete dynamic arrays
void DeleteData(parameters* para)
{
	delete[]para->bss;
	delete[]para->bss1;
	delete[]para->bss2;
	delete[]para->bss3;
	delete[]para->bss3plus;
	delete[]para->mua0;
	delete[]para->mus0;
	delete[]para->muat;

	if (para->it == I_CUSTOME)
	{
		delete[]para->i_m_rnd;
		delete[]para->i_angle;
	}
	/*if (para->it == I_MIE)
	{
		delete[]para->imie_m_rnd;
		delete[]para->imie_angle;
	}*/
	if (para->st == S_CUSTOME)
	{
		delete[]para->s_m_rnd;
		delete[]para->s_angle;
	}
}


//start the photon
void Launch(parameters* para)
{
	double TransmitCoeffcient;
	double cosIncAng;
	//GenGuaDis();
	lightLauch(para);
	//refraction on the sea surface
	Vector3 IncVector(para->ux, para->uy, para->uz);
	Vector3 SeaSurNorVector(0, 0, -1);
	cosIncAng = cos(para->inc_ang * FACTOR);
	Vector3 OutVector;
	//TransmitCoeffcient = Fresnel(para->airindex, para->seaindex, cosIncAng);

	TransmitCoeffcient = para->interface[1];

	OutVector = Refract(IncVector, SeaSurNorVector, para->airindex, para->seaindex);
	Vector3 bakVector = OutVector;
	OutVector.normalize();


	//
	para->ux = OutVector.x;
	para->uy = OutVector.y;
	para->uz = OutVector.z;

	para->status = ALIVE;
	para->hit = NOHIT;
	para->totalpath = 0;
	para->sleft = 0;
	para->weight = 1.0 * TransmitCoeffcient;
	para->collision = 0;
}


//output the simulation result to a file
void OutputToFile(parameters* para)
{
	ofstream SignalFile(para->custom_out);
	if (!SignalFile)
	{
		cerr << "cannot open signal.txt for output" << endl;
	}
	double SignalDepth;
	int lightspeed = LIGHTSPEED;
	for (int i = 0; i < para->BssArrayLength - 1; i++)
	{
		SignalDepth = 0.5 * para->ns_per_bin * 1e-9 * i * lightspeed / para->seaindex;
		SignalFile << fixed << setprecision(2) << SignalDepth << "\t";
		SignalFile << scientific << setprecision(3) << para->bss[i]/para->Nphotons << "\t" << para->bss1[i] / para->Nphotons << "\t" << para->bss2[i] / para->Nphotons << "\t" << para->bss3[i] / para->Nphotons << "\t" << para->bss3plus[i] / para->Nphotons << endl;
	}
}

//main simulation routine
int SimProc(parameters* para)
{
	std::cout << "\t\t\t\t3DOLSE" << endl;
	std::cout << "\t\t3D-OCEAN LIDAR SYSTEM EMULATOR\n";
	int i_photon; int photonspace; int tempphoton; int percent;
	_tzset();
	clock_t start, end;
	start = clock();
	Initpara(para);
	photonspace = para->Nphotons / 100;
	tempphoton = 0;
	for (i_photon = 1; i_photon <= para->Nphotons; i_photon++)
	{
		tempphoton++;
		if (tempphoton == photonspace)
		{
			percent = i_photon / photonspace;
			printf("..........%d%%\n", percent);
			tempphoton = 0;
		}
		Launch(para);
		while (para->status == ALIVE)
		{
			StepSize(para);
			HitBoundary(para);
			if (para->hit == NOHIT)
				NoHit(para);
			else if (para->hit == HITSEASURF)
				HitSeaSurf(para);
			else if (para->hit == HITSEABOT)
				HitSeaBot(para);
		}
	}
	OutputToFile(para);
	end = clock();
	double elapsed = (end - start) / CLOCKS_PER_SEC;
	printf(" %8.3lf seconds for entire simulation run \n", elapsed);
	DeleteData(para);
	DeleteArrays(para);
	return 0;
}
#endif

