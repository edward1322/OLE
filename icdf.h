#pragma once
#ifndef ICDF_H
#define ICDF_H
#include "parameters.h"
//Scatter and establish new direction
double ICDF_HG(double g, double rnd)
{
	double costheta;
	if (g == 0)
		costheta = 2.0*rnd - 1.0;
	else
	{
		double temp = (1.0 - g*g) / (1.0 - g + 2 * g*rnd);
		costheta = (1.0 + g*g - temp * temp) / (2.0*g);
	}
	return costheta;
}
/*
%%%%%%%%%%%%%%%%%%%%%LUT for find the angle%%%%%%%%%%%%%%%%%%%%%%%%%%%
*/

double ICDF_CUSTOME(double rnd, parameters* para)
{
	double m_ang = 0;
	double costheta;
	int i = 0;
	if (rnd <= para->i_m_rnd[0])
	{
		m_ang = 0;
		costheta = cos(m_ang);
		return costheta;
	}
	if (rnd >= para->i_m_rnd[para->i_len - 1])
	{
		m_ang = para->i_angle[para->i_len - 1];
		costheta = cos(m_ang);
		return costheta;
	}
	for (i = 0; i < para->i_len - 1; i++)
	{
		if ((para->i_m_rnd[i] < rnd) && (rnd < para->i_m_rnd[i + 1]))
		{
			m_ang = para->i_angle[i];
			break;
		}
	}
	costheta = cos(m_ang);
	return costheta;
}
#endif