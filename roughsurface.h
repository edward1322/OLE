#pragma once
#ifndef ROUGHSURFACE_H
#define ROUGHSURFACE_H
/*
	求风速W，入射角θ下的反射，透射率
	返回1*4数组，[Raa,Taw,Rww,Twa]
*/

#include<stdio.h>
#include<math.h>
#define NMAX 500  // 最大高斯节点
#define RIW 1.333    //refractive index of water  海水的折射率
#define NVAL 200 //使用的高斯节点

double scalar_RFAA(double theta);  // 菲涅尔air-incident反射率
double scalar_RFWW(double theta);  // 菲涅尔air-incident透射率
double scalar_TFWA(double theta);  // 菲涅尔water-incident反射率
double scalar_TFAW(double theta);  // 菲涅尔water-incident透射率 
double scalar_SHADOW(double mu, double mu_, double sigma2);  // SHADOW function
double scalar_PROB(double mu_n, double sigma2);   // PROB function

void scalar_GUASS(int N, double* Z, double* W);  //高斯点
double scalar_braa(double sigma2, double mu_i, double phi_i, double mu_r, double phi_r);  //raa（σ2，μ',φ',μ,φ）
double scalar_brww(double sigma2, double mu_i, double phi_i, double mu_r, double phi_r);  //rww（σ2，μ',φ',μ,φ）
double scalar_btaw(double sigma2, double mu_i, double phi_i, double mu_t, double phi_t);  //taw（σ2，μ',φ',μ,φ）
double scalar_btwa(double sigma2, double mu_i, double phi_i, double mu_t, double phi_t);  //twa（σ2，μ',φ',μ,φ）


void ocean_surface_scalar(double W, double theta, double result[4])
{
	double raa = 0, taw = 0, rww = 0, twa = 0;
	double pi = acos(-1);
	double sigma2 = 0.003 + W * 0.00512;
	double mu_i = cos(theta * pi / 180), phi_i = 0;
	double x[NMAX], w[NMAX];
	double x_mu_up[NVAL + 1], w_mu[NVAL + 1], x_mu_down[NVAL + 1];
	double x_phi[NVAL + 1], w_phi[NVAL + 1];
	double c = 0, d = 2 * pi;
	double raa_phi, twa_phi, rww_phi, taw_phi;

	scalar_GUASS(NVAL, x, w);
	for (int i = 1; i <= NVAL; i++)
	{
		x_mu_up[i] = 0.5 + 0.5 * x[i];
		w_mu[i] = 0.5 * w[i];
		x_mu_down[i] = -0.5 + 0.5 * x[i];
		x_phi[i] = (d + c) / 2 + (d - c) / 2 * x[i];
		w_phi[i] = (d - c) / 2 * w[i];
	}
	for (int i = 1; i <= NVAL; i++)
	{
		raa_phi = 0, twa_phi = 0, rww_phi = 0, taw_phi = 0;
		for (int j = 1; j <= NVAL; j++)
		{
			raa_phi += w_phi[j] * scalar_braa(sigma2, mu_i, phi_i, x_mu_up[i], x_phi[j]) * fabs(x_mu_up[i]);
			twa_phi += w_phi[j] * scalar_btwa(sigma2, mu_i, phi_i, x_mu_up[i], x_phi[j]) * fabs(x_mu_up[i]);
			rww_phi += w_phi[j] * scalar_brww(sigma2, mu_i, phi_i, x_mu_down[i], x_phi[j]) * fabs(x_mu_down[i]);
			taw_phi += w_phi[j] * scalar_btaw(sigma2, mu_i, phi_i, x_mu_down[i], x_phi[j]) * fabs(x_mu_down[i]);
		}
		raa += w_mu[i] * raa_phi;
		twa += w_mu[i] * twa_phi;
		rww += w_mu[i] * rww_phi;
		taw += w_mu[i] * taw_phi;
	}
	result[0] = raa / pi, result[1] = taw / pi, result[2] = rww / pi, result[3] = twa / pi;
}

double scalar_braa(double sigma2, double mu_i, double phi_i, double mu_r, double phi_r)
{
	mu_i = -fabs(mu_i); mu_r = fabs(mu_r);
	double pi = acos(-1);
	double mu = -mu_r, mu_ = -mu_i;
	double sinr = sqrt(1 - pow(mu, 2)), sini = sqrt(1 - pow(mu_, 2));
	double cosTHETA = mu * mu_ + sini * sinr * cos(phi_r - phi_i);
	double THETA = acos(cosTHETA);
	double theta = (pi - THETA) / 2;
	double mu_n = fabs(mu - mu_) / sqrt(2 - 2 * cosTHETA);
	double r = scalar_RFAA(theta);
	double s = scalar_SHADOW(mu, mu_, sigma2);
	double p = scalar_PROB(mu_n, sigma2);

	double br = s * pi * p / (4 * fabs(mu) * mu_n * fabs(mu_)) * r;

	return br;
}

double scalar_RFAA(double theta)
{
	double pi = acos(-1);
	double m = RIW;
	double c1 = pow(m, 2) * cos(theta), c2 = sqrt(pow(m, 2) - pow(sin(theta), 2));
	double rl = (c1 - c2) / (c1 + c2);
	c1 = cos(theta); c2 = sqrt(pow(m, 2) - pow(sin(theta), 2));
	double rr = (c1 - c2) / (c1 + c2);
	double r = 0.5 * (pow(rl, 2) + pow(rr, 2));
	return r;
}

double scalar_brww(double sigma2, double mu_i, double phi_i, double mu_r, double phi_r)
{
	mu_i = fabs(mu_i); mu_r = -fabs(mu_r);
	double pi = acos(-1);
	double mu = -mu_r, mu_ = -mu_i;
	double sinr = sqrt(1 - pow(mu, 2)), sini = sqrt(1 - pow(mu_, 2));
	double cosTHETA = mu * mu_ + sini * sinr * cos(phi_r - phi_i);
	double THETA = acos(cosTHETA);
	double theta = (pi - THETA) / 2;
	double mu_n = fabs(mu - mu_) / sqrt(2 - 2 * cosTHETA);
	double r = scalar_RFWW(theta);
	double s = scalar_SHADOW(mu, mu_, sigma2);
	double p = scalar_PROB(mu_n, sigma2);

	double br = s * pi * p / (4 * fabs(mu) * mu_n * fabs(mu_)) * r;

	return br;
}

double scalar_RFWW(double theta)
{
	double pi = acos(-1);
	double m = RIW;
	double theta_lim = asin(1 / m);
	double c1 = cos(theta), c2 = m * sqrt(1 - pow(m, 2) * pow(sin(theta), 2));
	double rl = (c1 - c2) / (c1 + c2);
	c1 = m * cos(theta); c2 = sqrt(1 - pow(m, 2) * pow(sin(theta), 2));
	double rr = (c1 - c2) / (c1 + c2);
	double r = 0.5 * (pow(rl, 2) + pow(rr, 2));
	if (theta > theta_lim) r = 1;
	return r;
}

double scalar_btaw(double sigma2, double mu_i, double phi_i, double mu_t, double phi_t)
{
	mu_i = -fabs(mu_i), mu_t = -fabs(mu_t);
	double pi = acos(-1);
	double cn1 = 1, cn2 = RIW;
	double mu = -mu_t, mu_ = -mu_i;
	double sint_ = sqrt(1 - pow(mu, 2)), sini_ = sqrt(1 - pow(mu_, 2));
	double cosTHETA = mu * mu_ + sini_ * sint_ * cos(phi_t - phi_i);
	double C = sqrt(pow(cn2, 2) + 1 - 2 * cn2 * cosTHETA);
	double cosi = (cn2 * cosTHETA - 1) / C;
	double cost = (cn2 - cosTHETA) / C;
	double thetai = acos(cosi);
	double mu_n = -(mu_ - cn2 * mu) / C;

	double t = scalar_TFAW(thetai);
	if (cosi < 0) t = 0;

	double s = scalar_SHADOW(mu, mu_, sigma2);
	double g = (pow(cn2, 2) * cost * cosi) / pow((cn2 * cost - cn1 * cosi), 2);
	double gT = scalar_PROB(mu_n, sigma2) / (fabs(mu) * mu_n);
	double bt = s * pi * t * gT * g / (fabs(mu_));
	return bt;
}

double scalar_TFAW(double theta)
{
	double m = RIW;
	double pi = acos(-1);
	double c1 = pow(m, 2) * cos(theta), c2 = sqrt(pow(m, 2) - pow(sin(theta), 2));
	double rl = (c1 - c2) / (c1 + c2);
	c1 = cos(theta); c2 = sqrt(pow(m, 2) - pow(sin(theta), 2));
	double rr = (c1 - c2) / (c1 + c2);

	double tl = 1 / m * (1 + rl), tr = 1 + rr;
	double t = 0.5 * (pow(tl, 2) + pow(tr, 2));
	double c = m * sqrt(1 - pow((sin(theta) / m), 2)) / cos(theta);
	double TF = c * t;
	return TF;
}
double scalar_btwa(double sigma2, double mu_i, double phi_i, double mu_t, double phi_t)
{
	mu_i = fabs(mu_i), mu_t = fabs(mu_t);
	double pi = acos(-1);
	double cn2 = RIW;
	double mu = -mu_t, mu_ = -mu_i;
	double sint_ = sqrt(1 - pow(mu, 2)), sini_ = sqrt(1 - pow(mu_, 2));
	double cosTHETA = mu * mu_ + sini_ * sint_ * cos(phi_t - phi_i);
	double C = sqrt(pow(cn2, 2) + 1 - 2 * cn2 * cosTHETA);
	double cosi = (cn2 - cosTHETA) / C;
	double cost = (cn2 * cosTHETA - 1) / C;
	double thetai = acos(cosi);
	double mu_n = -(cn2 * mu_ - mu) / C;

	double t = scalar_TFWA(thetai);
	if (cosi < 0) t = 0;

	double s = scalar_SHADOW(mu, mu_, sigma2);
	double g = (cost * cosi) / pow((cn2 * cosi - cost), 2);
	double gT = scalar_PROB(mu_n, sigma2) / (fabs(mu) * mu_n);
	double bt = s * pi * t * gT * g / (fabs(mu_));
	return bt;
}

double scalar_TFWA(double theta)
{
	double m = RIW;
	double theta_lim = asin(1 / m);
	double c1 = cos(theta), c2 = m * sqrt(1 - pow(m, 2) * pow(sin(theta), 2));
	double rl = (c1 - c2) / (c1 + c2);
	c1 = m * cos(theta), c2 = sqrt(1 - pow(m, 2) * pow(sin(theta), 2));
	double rr = (c1 - c2) / (c1 + c2);

	double cosTheta_t = sqrt(1 - pow(m * sin(theta), 2));
	if (theta > theta_lim) cosTheta_t = 0;

	double tl = m * (1 + rl), tr = 1 + rr;
	double t = 0.5 * (pow(tl, 2) + pow(tr, 2));
	double c = cosTheta_t / (m * cos(theta));
	double TF = c * t;
	return TF;
}

double scalar_A(double mu, double sigma2)
{
	double pi = acos(-1);
	double eta = mu / (sqrt(sigma2) * sqrt(1 - pow(mu, 2)));
	double a = 0.5 * (1 / (sqrt(pi) * eta) * exp(-pow(eta, 2)) - erfc(eta));
	return a;
}
double scalar_SHADOW(double mu, double mu_, double sigma2)
{
	double s = 1 / (1 + scalar_A(fabs(mu), sigma2) + scalar_A(fabs(mu_), sigma2));
	return s;

}


double scalar_PROB(double mu_n, double sigma2)
{
	double pi = acos(-1);
	double expval = exp(-(1 - pow(mu_n, 2)) / (sigma2 * pow(mu_n, 2)));
	double a = 1 / (pi * sigma2 * pow(mu_n, 3));
	double p = a * expval;
	return p;

}

/* Cal gauss points*/
void scalar_GUASS(int N, double* Z, double* W)
{
	double A = 1, B = 2, C = 3;
	int IND = N % 2;
	int K = N / 2 + IND;
	int M;
	double F = N;
	double X, PA, PB, PC, DJ;
	for (int I = 1; I <= K; I = I + 1)
	{
		M = N + 1 - I;
		if (I == 1)X = A - B / ((F + A) * F);
		else if (I == 2) X = (Z[N] - A) * 4 + Z[N];
		else if (I == 3) X = (Z[N - 1] - Z[N]) * 1.6 + Z[N - 1];
		else X = (Z[M + 1] - Z[M + 2]) * C + Z[M + 3];
		if (I == K && IND == 1) X = 0;
		while (1)
		{
			PB = 1;
			PC = X;
			DJ = A;
			for (int J = 2; J <= N; J += 1)
			{
				DJ = DJ + A;
				PA = PB;
				PB = PC;
				PC = X * PB + (X * PB - PA) * (DJ - A) / DJ;
			}
			PA = A / ((PB - X * PC) * F);
			PB = PA * PC * (A - X * X);
			X = X - PB;
			if (fabs(PB) <= (1e-12 * fabs(X))) break;
		}
		Z[M] = X;
		W[M] = PA * PA * (A - X * X);
		W[M] *= B;
		if (I == K && IND == 1) break;
		Z[I] = -Z[M];
		W[I] = W[M];
	}

	return;
}

#endif