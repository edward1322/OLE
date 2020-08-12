#pragma once
#ifndef INTERFACE_H
#define INTERFACE_H

#include <stdlib.h>
#include "math.h"
#include "Vector3.h"


using namespace std;

//subroutine to compute the Fresnel transmission coefficient
double Fresnel(double n1, double n2, double CosIncAng)
{
	double r;
	CosIncAng = fabs(CosIncAng);
	if (n1 == n2) { /** matched boundary. **/
		r = 0.0;
	}
	else if (CosIncAng > (1.0 - 1.0e-12)) { /** normal incidence. **/
		r = (n2 - n1) / (n2 + n1);
		r *= r;
	}
	else if (CosIncAng < 1.0e-6) {	/** very slanted. **/
		r = 1.0;
	}
	else {			  		/** general. **/
		double sa1, sa2; /* sine of incident and transmission angles. */
		double ca2;      /* cosine of transmission angle. */
		sa1 = sqrt(1 - CosIncAng * CosIncAng);
		sa2 = n1 * sa1 / n2;
		if (sa2 >= 1.0) {
			/* double check for total internal reflection. */
			r = 1.0;
		}
		else {
			double cap, cam;	/* cosines of sum ap or diff am of the two */
								/* angles: ap = a1 + a2, am = a1 - a2. */
			double sap, sam;	/* sines. */
			ca2 = sqrt(1 - sa2 * sa2);
			cap = CosIncAng * ca2 - sa1 * sa2; /* c+ = cc - ss. */
			cam = CosIncAng * ca2 + sa1 * sa2; /* c- = cc + ss. */
			sap = sa1 * ca2 + CosIncAng * sa2; /* s+ = sc + cs. */
			sam = sa1 * ca2 - CosIncAng * sa2; /* s- = sc - cs. */
			r = 0.5*sam*sam*(cam*cam + cap * cap) / (sap*sap*cam*cam);
			/* rearranged for speed. */
		}
	}
	return(1 - r);
}
//subroutine to handle the reflection
Vector3 Reflect(Vector3 IncVector, Vector3 NorVector)
{
	Vector3 RefVector = IncVector - 2.0*NorVector*(NorVector*IncVector / (vectorMag(NorVector)*vectorMag(IncVector)));
	return RefVector;
}
//subroutine to handle the refraction
Vector3 Refract(Vector3 IncVector, Vector3 NorVector, double ni, double no)
{
	double NRation = ni / no;
	double cosI = (-IncVector)*NorVector / (vectorMag(NorVector)*vectorMag(IncVector));
	double cosT2 = 1.0 - NRation * NRation*(1 - cosI * cosI);
	Vector3 TraVector = NRation * IncVector + ((NRation*cosI - sqrt(fabs(cosT2)))*NorVector);
	return TraVector * (cosT2 > 0);
}
#endif