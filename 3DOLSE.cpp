// 3DOLSE.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//
#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include "parameters.h"
#include "Simpro.h"
#include <fstream>

void Loadpara(parameters* para);
int main()
{
	parameters* para = new parameters();	
	Loadpara(para);
	SimProc(para);	
}
//Copy input data to 'para' variable
void Loadpara(parameters* para)
{
	// seawater parameters
	para->g0 = 0.924;
	para->ffn = 1.10;
	para->ffmu = 3.62;
	para->seaindex = 1.333;
	para->airindex = 1.0;
	para->seadepth = 30;
	para->laserWavel=532;
	para->fseabott = 0.01;
	para->fseabotttype = rgreen;
	para->bo = opt;
	para->st = FF; //H_G, MHG, TTHG, FF, MIE_CAL,S_CUSTOME
	para->it = I_HG; //I_HG, I_CUSTOME
	// lidar parameters
	para->planeheight = 300;
	para->viewfield = 0.2;
	para->receivearea = 0.09;
	para->inc_ang = 0;
	para->light_R = 0;
	para->light_typ = 0;
	para->Nphotons =  1*1e7;
	para->THRESHOLD = 1.0e-4;
	para->ns_per_bin = 1;
	para->BssArrayLength = 300;
	para->windspeed = 0;	
	para->startWavel = 532;
	para->endWavel = 532;
	para->stepWavel = 1;
	para->scatRefReal = 1.379;
	para->scatRefImag = 0.0;
	para->medRef = 1.333;
	para->Disperse = PolyDisperse;  //MonoDisperse , PolyDisperse
	para->distIndex = Junge;   //Log_Normal, Guassian, Junge, CustomData
	para->m_type = VolFrac;                         //Conc_mm3, VolFrac
	para->custom_psd = ".\\def\\psd.txt";
	para->custom_spf = ".\\def\\spf.txt";
	para->custom_cdf = ".\\def\\cdf.txt";
	para->custom_in = ".\\in\\opt.txt";
	para->custom_out = ".\\out\\echo.txt";
	para->seabotfile = ".\\def\\SEA_BOTTOM_REFLECTANCES.txt";

	if (para->m_type == Conc_mm3)
		para->sphNumDensity = 1e+08;
	if (para->m_type == VolFrac)
		para->volFraction = 1 / 1.5 * 1e-7;   //volFraction 10e-1 = 1.5mg/L

	if (para->Disperse == MonoDisperse)
	{
		para->meanRadius = 1.0 / 2.0; //meanRadius
		para->nRadius = 1;
		if (para->m_type == VolFrac)
		{
			double volume = 4.0 * M_PI * para->meanRadius * para->meanRadius * para->meanRadius / 3.0;
			para->sphNumDensity = para->volFraction * 1e9 / volume;
		}
	}
	if ((para->Disperse == PolyDisperse) && (para->distIndex != CustomData))
	{
		para->meanRadius = 2 / 2.0; //MeanDiameter
		para->stdDev = 0.25; //StdDev
		para->nRadius = 15; //NSphere
	}
	para->fRay = 0.00;
	para->bMie = 0.51;
	//para->nTheta = 101;
	para->DTheta = Phase_DTheta_1; //Phase_DTheta0_1(0.1°), Phase_DTheta0_5(0.5°), Phase_DTheta_1(1°)
}