#pragma once
#ifndef LIGHT_H
#define LIGHT_H
#include "drand48.h"
void GenGuaDis(parameters* para)
{
	double GuassDisR;
	double RandomTheta;
	double tempR;
	double rnd = RandomNum();
	GuassDisR = 0.5* para->light_R*sqrt(-2 * log(1 - rnd));
	rnd = RandomNum();
	RandomTheta = 2 * PI*rnd;
	para->x = GuassDisR * cos(RandomTheta);
	para->y = GuassDisR * sin(RandomTheta);
	para->z = 0;
	tempR = sqrt(para->x* para->x + para->y* para->y + para->planeheight* para->planeheight);
	para->ux = para->x / tempR;
	para->uy = para->y / tempR;
	para->uz = para->planeheight / tempR;
}
void lightLauch(parameters* para)
{
	double rnd = RandomNum();
	double radius = para->light_R;
	if (para->light_typ == 0) {
		/* UNIFORM COLLIMATED BEAM */
		
		para->x = radius * sqrt(rnd);
		para->y = 0;
		para->z = 0;
		/* Initial trajectory as cosines */
		para->ux = 0;
		para->uy = 0;
		para->uz = cos(para->inc_ang*FACTOR);
	}
	else if (para->light_typ == 1) {
		/* COLLIMATED GAUSSIAN BEAM*/	

		double GuassDisR;
		double RandomTheta;
		double tempR;
		double rnd = RandomNum();
		/*double x = 1/exp(1);
		GuassDisR = x*theInputData.R*sqrt(-log(rnd));*/
		GuassDisR = 0.5 * para->light_R * sqrt(-2 * log(1 - rnd));
		rnd = RandomNum();
		RandomTheta = 2 * PI*rnd;
		para->x = GuassDisR * cos(RandomTheta);
		para->y = GuassDisR * sin(RandomTheta);
		para->z = 0;
		double costmp = cos(para->inc_ang / 180 * PI)*cos(para->inc_ang / 180 * PI);
		tempR = sqrt(para->x* para->x + para->y* para->y + para->planeheight* para->planeheight/costmp);
		para->ux = para->x / tempR;
		para->uy = para->y / tempR;
		para->uz = para->planeheight / tempR;
	}
	else if (para->light_typ == 2) {
		/* ISOTROPIC POINT SOURCE*/
		/* Initial position */
		para->x = 0;
		para->y = 0;
		para->z = 0;
		/* Initial trajectory as cosines */
		double costheta = 1.0 - 2.0*rnd;
		double sintheta = sqrt(1.0 - costheta * costheta);
		rnd = RandomNum();
		double psi = 2.0*PI*rnd;
		double cospsi = cos(psi);
		double sinpsi;
		if (psi < PI)
			sinpsi = sqrt(1.0 - cospsi * cospsi);
		else
			sinpsi = -sqrt(1.0 - cospsi * cospsi);
		para->ux = sintheta * cospsi;
		para->uy = sintheta * sinpsi;
		para->uz = costheta;
	}
	else {
		printf("choose mcflag between 0 to 3\n");
	}
}
#endif