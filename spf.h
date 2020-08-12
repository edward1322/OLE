#pragma once
#ifndef SPF_H
#define SPF_H
#include <math.h>
#include "parameters.h"
//subroutine to compute the spf function
double HG(double costheta, double g)
{
	double temp, ptheta;
	temp = 1.0 + g * g - 2.0 * g * costheta;
	ptheta = (1.0 - g * g) / 4.0 / PI / (temp * sqrt(temp));
	return ptheta;
}
double ModHG(double costheta, double g)
{
	double temp, ptheta;
	temp = 1.0 + g * g - 2.0 * g * costheta;
	ptheta = 1.5 * (1.0 - g * g) / (2.0 + g * g) * (1.0 + costheta * costheta) / (temp * sqrt(temp));
	return ptheta;
}
double TT_HG(double costheta, double g)
{
	double ptheta;
	double g1 = g;
	double g2 = -0.30614 + 1.0006 * g1 - 0.01826 * pow(g1, 2) + 0.03644 * pow(g1, 3);
	double temp1 = 1 + g1 * g1 - 2 * g1 * costheta;
	double ptheta1 = (1 - g1 * g1) / 4.0 / PI / (temp1 * sqrt(temp1));
	double temp2 = 1 + g2 * g2 + 2 * g2 * costheta;
	double ptheta2 = (1 - g2 * g2) / 4.0 / PI / (temp2 * sqrt(temp2));
	double a = g2 * (1 + g2) / (g1 + g2) / (1 + g2 - g1);
	ptheta = a * ptheta1 + (1 - a) * ptheta2;
	return ptheta;
}
double Fournier_Forand(double costheta, double g, double ffn, double ffmu)
{
	float ptheta = 0.0;
	if (costheta >= cos(2.0 / 180 * PI))
	{
		ptheta = HG(costheta, g);
		return ptheta;
	}
	
	float n = ffn;
	float mu = ffmu;
	float nu = (3.0 - mu) / 2.0;
	float rad_ang = acos(costheta);
	if (isnan(rad_ang))
	{
		rad_ang = PI;
	}

	float delta = 4.0 / 3 / pow(n - 1, 2) * pow(sin(rad_ang / 2), 2);
	float delta_180 = 4.0 / 3 / pow(n - 1, 2) * pow(sin(PI / 2), 2);
	float tm = abs(1.0 - delta);
	float t = (pow(1 - delta, 2) * pow(delta, nu)); //delta=1,t->0;此时，costheta->0.9850;故要去掉这个点
	float p1 = 1.0 / 4.0 / PI / t;

	float t1 = nu * (1 - delta) - (1 - pow(delta, nu));
	float tmp = pow(sin(rad_ang / 2), -2);
	float t2 = (delta * (1 - pow(delta, nu)) - nu * (1 - delta)) * tmp;
	float p2 = t1 + t2;

	float s1 = (1 - pow(delta_180, nu)) / 16.0 / PI / ((delta_180 - 1) * pow(delta_180, nu));
	float s2 = 3 * pow(cos(rad_ang), 2) - 1;
	float p3 = s1 * s2;

	ptheta = p1 * p2 + p3;
	if (ptheta < 0)
		ptheta = 0.0;
	/*double ang;
	if (ptheta >= 25)
		ang = rad_ang / PI * 180;*/ //测试用
	double x = abs(costheta - 0.9850);
	if (x <= 0.001)
	{
		ptheta = HG(costheta, g);
	}
	return ptheta;
}
double max(double a, double b)
{
	double p;
	p = a;
	if (b > a)p = b;
	return p;
}
double min(double a, double b)
{
	double p;
	p = a;
	if (b < a)p = b;
	return p;
}
/*
		caculated by lut slope
*/
double SPF_CUSTOME(double costheta, parameters* para)
{
	float ptheta = 0.0;
	if (costheta >= cos(2.0 / 180 * PI))
	{
		ptheta = HG(costheta, para->g0);
		return ptheta;
	}
	double theta = acos(costheta);

	if (theta <= para->s_angle[0])
	{
		ptheta = para->s_m_rnd[0];
		return ptheta;
	}
	if (theta >= para->s_angle[para->s_len - 1])
	{
		ptheta = para->s_m_rnd[para->s_len - 1];
		return ptheta;
	}
	for (int i = 0; i < para->s_len - 1; i++)
	{
		if ((para->s_angle[i] <= theta) && (theta <= para->s_angle[i + 1]))
		{
			//ptheta = 0.5 * (para->s_m_rnd[i + 1] + para->s_m_rnd[i]);
			ptheta = (para->s_angle[i + 1] - theta) / (para->s_angle[i + 1] - para->s_angle[i]) * abs(para->s_m_rnd[i + 1] - para->s_m_rnd[i])+ min(para->s_m_rnd[i + 1], para->s_m_rnd[i]);
			//ptheta = max(para->s_m_rnd[i + 1], para->s_m_rnd[i]);
			break;
		}
	}
	return ptheta;
}
/*
		caculated by mie lut slope
*/
double MIE_CALF(double costheta, parameters* para)
{
	float ptheta = 0.0;
	//if (costheta >= cos(2.0 / 180 * PI))
	//{
	//	ptheta = HG(costheta, para->g0);
	//	return ptheta;
	//}
	double theta = acos(costheta);

	if (theta <= para->minTheta)
	{
		para->phaseFunctionAve[0][0];
		return ptheta;
	}
	if (theta >= para->maxTheta)
	{
		para->phaseFunctionAve[0][para->nTheta - 1];
		return ptheta;
	}
	double p_down=0;
	double p_up = 0;
	for (int t = 0; t < para->nTheta-1; t++)
	{
		p_down = para->minTheta + t * para->stepTheta;
		p_up = para->minTheta + (t+1) * para->stepTheta;
		if ((p_down <= theta)&&(theta <= p_up))
		{
			//ptheta = 0.5 * (para->s_m_rnd[i + 1] + para->s_m_rnd[i]);
			ptheta = (p_up - theta) / (p_up - p_down) * abs(para->phaseFunctionAve[0][t+1] - para->phaseFunctionAve[0][t]) + min(para->phaseFunctionAve[0][t + 1] , para->phaseFunctionAve[0][t]);
			break;
		}			
	}	
	return ptheta;
}

#endif