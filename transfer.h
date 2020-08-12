#pragma once
#ifndef TRANSFER_H
#define TRANSFER_H
#include "parameters.h"
#include "icdf.h"
#include "spf.h"
#include "roughsurface.h"

int SIGN(double t)
{
	if (t > 0)
		return 1;
	else
		return -1;
}
//generate the stepsize of the photon to move
void StepSize(parameters* para)
{
	double rnd;
	float depth_i = 0.5 * para->ns_per_bin * 1e-9 * LIGHTSPEED / para->seaindex;
	int sea_i = round(para->z / depth_i);
	/*if (sea_i == 0)
	{
		sea_i += 1;
	}*/
	if (para->sleft == 0.0)
	{
		rnd = RandomNum();
		para->steplength = -log(rnd) / para->muat[sea_i];
	}
	else
	{
		para->steplength = para->sleft;
		para->sleft = 0;
	}
}
//Move to next scattering or absorption event
void Move(parameters* para)
{
	para->x += para->steplength * para->ux;
	para->y += para->steplength * para->uy;
	para->z += para->steplength * para->uz;
	para->totalpath += para->steplength;
}
//decide the photon hit the sea surface and the sea bottom of sea or not
void HitBoundary(parameters* para)
{
	double dl2b; //length to boundary;
	if (para->uz > 0.0)
		dl2b = (para->seadepth - para->z) / para->uz;  //dl2b>0;
	else if (para->uz < 0.0)
		dl2b = -para->z / para->uz;    //dl2b>0;
	if (para->uz != 0.0 && para->steplength > dl2b)   //not horizontal &//crossing
	{
		para->sleft = para->steplength - dl2b;
		para->steplength = dl2b;
		if (para->uz > 0)
			para->hit = HITSEABOT;
		else
			para->hit = HITSEASURF;
	}
	else
		para->hit = NOHIT;
}
//The photon weight is small,and the photon packet tries to survive a roulette
void Roulette(parameters* para)
{
	if (para->weight < para->THRESHOLD) {
		double rnd = RandomNum();
		if (rnd <= CHANCE)
			para->weight /= CHANCE;
		else para->status = DEAD;
	}
}
//Absorb
void Absorb(parameters* para)
{
	double albedo;
	float depth_i = 0.5 * para->ns_per_bin * 1e-9 * LIGHTSPEED / para->seaindex;
	int sea_i = round(para->z / depth_i);
	/*if (sea_i == 0)
	{
		sea_i += 1;
	}*/
	albedo = para->mus0[sea_i] / para->muat[sea_i];
	para->weight *= albedo;
}

void ScatterOrReflect(parameters* para)
{
	double costheta;
	double sintheta;
	double cospsi;
	double sinpsi;
	double psi;
	double uxx, uyy, uzz; double rnd1;
	double rnd = RandomNum();
	if (para->hit == NOHIT)
	{
		para->collision++;
		switch (para->it)
		{
		case I_HG:
			costheta = ICDF_HG(para->g0, rnd);
			break;
		case I_CUSTOME:
			costheta = ICDF_CUSTOME(rnd,para);
			break;
		default:
			costheta = ICDF_HG(para->g0, rnd);
			break;
		}		
		sintheta = sqrt(1.0 - costheta * costheta);
	}
	else if (para->hit == HITSEABOT)
	{
		costheta = -sqrt(1 - rnd);
		sintheta = sqrt(rnd);
	}
	rnd1 = RandomNum();
	psi = 2.0 * PI * rnd1; cospsi = cos(psi);
	if (psi < PI)
		sinpsi = sqrt(1.0 - cospsi * cospsi);
	else
		sinpsi = -sqrt(1.0 - cospsi * cospsi);
	/*****new trajectory*****/
	if (1 - fabs(para->uz) <= ONE_MINUS_COSZERO)
	{
		uxx = sintheta * cospsi;
		uyy = sintheta * sinpsi;
		uzz = costheta * SIGN(para->uz);
	}
	else
	{
		double temp = sqrt(1.0 - para->uz * para->uz);
		uxx = sintheta * (para->ux * para->uz * cospsi - para->uy * sinpsi) / temp + para->ux * costheta;
		uyy = sintheta * (para->uy * para->uz * cospsi + para->ux * sinpsi) / temp + para->uy * costheta;
		uzz = -sintheta * cospsi * temp + para->uz * costheta;
	}
	/*update trajectory*/
	para->ux = uxx;
	para->uy = uyy;
	para->uz = uzz;
	if (isnan(uzz))
		printf("trace error!");
}
//state the photons received by the detector which escape from the sea surface
void StatSurface(parameters* para, const Vector3& OutVector)
{
	int bin;
	double TransmitCoefficient;
	double stemp, xtemp, ytemp, rtemp;
	//TransmitCoefficient = Fresnel(theSea.seaindex, theSea.airindex, thePhoton.uz);

	TransmitCoefficient = para->interface[3];
	double rnd = RandomNum();
	if (rnd < TransmitCoefficient)
	{
		para->status = DEAD;
		double rfieldview; double rseasurf;
		rfieldview = para->planeheight * tan(para->viewfield / 2.0);
		rseasurf = sqrt(para->x * para->x + para->y * para->y);
		if (rseasurf < rfieldview)
		{
			stemp = para->planeheight / fabs(OutVector.z);
			xtemp = para->x + stemp * OutVector.x;
			ytemp = para->y + stemp * OutVector.y;
			rtemp = xtemp * xtemp + ytemp * ytemp;
			int lightspeed = LIGHTSPEED;
			if (rtemp < para->receivearea / PI)
			{
			    //int lightspeed = LIGHTSPEED;
				bin = (para->totalpath * para->seaindex / lightspeed) * 1e9 / para->ns_per_bin;
				if (bin >= para->BssArrayLength)
					bin = para->BssArrayLength - 1;
				para->bss[bin] += para->weight;
				if (para->collision == 1)
					para->bss1[bin] += para->weight;
				if (para->collision == 2)
					para->bss2[bin] += para->weight;
				if (para->collision == 3)
					para->bss3[bin] += para->weight;
				if (para->collision > 3)
					para->bss3plus[bin] += para->weight;
			}
		}
	}
}

//state the lidar received signal at each collision position of the photon
void StatWater(parameters* para, int hit)
{
	int bin; double costheta; double ptheta = 0.0;
	double TransmitCoeffcient;
	double h;//the square of the distance between the photon and the detector
	double d;//distance between photon and sea surface
	double l;//distance surface and detector;
	double NRation = para->seaindex / para->airindex;
	Vector3 PhoVector(para->ux, para->uy, para->uz);
	double rPhoton2 = para->x * para->x + para->y * para->y;
	double rPhoton = sqrt(rPhoton2);

	float depth_i = 0.5 * para->ns_per_bin * 1e-9 * LIGHTSPEED / para->seaindex;
	int sea_i = round(para->z / depth_i);

	//when the photon is in the interaction volume
	if (rPhoton < (para->seadepth + para->planeheight) * tan(para->viewfield / 2))
	{
		if (rPhoton == 0)
		{
			Vector3 BacVector(0, 0, -1);
			costheta = BacVector * PhoVector / (vectorMag(BacVector) * vectorMag(PhoVector));
			switch (para->st)
			{
			case H_G:
				ptheta = HG(costheta, para->g0);
				break;
			case MHG:
				ptheta = ModHG(costheta, para->g0);
				break;
			case TTHG:
				ptheta = TT_HG(costheta, para->g0);
				break;
			case FF:
				ptheta = Fournier_Forand(costheta, para->g0, para->ffn, para->ffmu);
				break;
			case MIE_CAL:
				ptheta = MIE_CALF(costheta, para);
				break;
			case S_CUSTOME:
				ptheta = SPF_CUSTOME(costheta, para);
				break;
			default:
				ptheta = HG(costheta, para->g0);
				break;
			}

			TransmitCoeffcient = para->interface[3];
			h = (para->z + para->planeheight) * (para->z + para->planeheight);
			d = para->z;
			l = para->planeheight;
		}
		else
		{
			double r = para->planeheight * NRation * rPhoton / (para->z + para->planeheight * NRation);
			double rx = r * para->x / rPhoton;
			double ry = r * para->y / rPhoton;
			Vector3 BacVector(rx - para->x, ry - para->y, -para->z);
			BacVector.normalize();
			TransmitCoeffcient = para->interface[3];
			costheta = PhoVector * BacVector / (vectorMag(BacVector) * vectorMag(PhoVector));
			switch (para->st)
			{
			case H_G:
				ptheta = HG(costheta, para->g0);
				break;
			case MHG:
				ptheta = ModHG(costheta, para->g0);
				break;
			case TTHG:
				ptheta = TT_HG(costheta, para->g0);
				break;
			case FF:
				ptheta = Fournier_Forand(costheta, para->g0, para->ffn, para->ffmu);
				break;
			case MIE_CAL:
				ptheta = MIE_CALF(costheta, para);
				break;
			case S_CUSTOME:
				ptheta = SPF_CUSTOME(costheta, para);
				break;
			default:
				ptheta = HG(costheta, para->g0);
				break;
			}

			h = rPhoton2 + (para->z + para->planeheight) * (para->z + para->planeheight);
			d = sqrt(para->z * para->z + (rPhoton - r) * (rPhoton - r));
			l = sqrt(r * r + para->planeheight * para->planeheight);
		}
		//compute the returned the weight at each collison postion
		double avg = 0.0;
		for (int i = 0; i <= sea_i; i++)
		{
			avg += para->muat[i];
		}
		avg /= (sea_i+1);

		double pweight = ptheta * para->receivearea / h * exp(-avg * d) * TransmitCoeffcient * para->weight;
		para->weight -= pweight;
		int lightspeed = LIGHTSPEED;
		bin = (d * para->seaindex / lightspeed + para->totalpath * para->seaindex / lightspeed) * 1e9 / para->ns_per_bin;
		if (bin >= para->BssArrayLength)
			bin = para->BssArrayLength - 1;
		para->bss[bin] += pweight;
		if (para->collision == 1)
			para->bss1[bin] += pweight;
		if (para->collision == 2)
			para->bss2[bin] += pweight;
		if (para->collision == 3)
			para->bss3[bin] += pweight;
		if (para->collision > 3)
			para->bss3plus[bin] += pweight;

	}
}
//subroutine to handle the photon which hit the seasurface
void HitSeaSurf(parameters* para)
{
	int hit = 1;
	Vector3 IncVector(para->ux, para->uy, para->uz);
	Vector3 OutVector, BacVector;
	Vector3 SeaSurNorVector(0, 0, 1);
	double cosI, sinI;
	cosI = IncVector * SeaSurNorVector / (vectorMag(IncVector) * vectorMag(SeaSurNorVector));
	sinI = sqrt(1.0 - cosI * cosI);
	Move(para);
	Absorb(para);
	if (sinI < para->airindex / para->seaindex) //check for total internal reflection
	{
		OutVector = Refract(IncVector, SeaSurNorVector, para->seaindex, para->airindex);
		OutVector.normalize();
		StatSurface(para, OutVector);
	}
	BacVector = Reflect(IncVector, SeaSurNorVector);
	para->ux = BacVector.x;
	para->uy = BacVector.y;
	para->uz = BacVector.z;
}
//subroutine to handle the photon which didn hit
void NoHit(parameters* para)
{
	int hit = 0;
	Move(para); Absorb(para);
	ScatterOrReflect(para); //需放在StatWater前
	StatWater(para, hit);
	Roulette(para);

}
//subroutine to handle the photon which hit the seabottom
void HitSeaBot(parameters* para)
{
	int hit = 2;
	Move(para); Absorb(para);
	ScatterOrReflect(para);
	para->weight *= para->fseabott;
	StatWater(para, hit);
	
}

#endif