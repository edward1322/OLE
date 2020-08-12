#pragma once
#ifndef DRAND48_H
#define DRAND48_H

#include <stdlib.h>
#include <stdio.h>
#define mm 0x100000000LL
#define cc 0xB16
#define aa 0x5DEECE66DLL

static unsigned long long seed = 1;

double drand48(void)
{
	seed = (aa * seed + cc) & 0xFFFFFFFFFFFFLL;
	unsigned int x = seed >> 16;
	return 	((double)x / (double)mm);

}


/*********************************************************************
*      RANDOM NUMBER GENERATOR
*      A random number generator that generates uniformly
*      distributed random numbers between 0 and 1 inclusive.
*      The algorithm is based on:
*      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
*      Flannery, "Numerical Recipes in C," Cambridge University
*      Press, 2nd edition, (1992).
*      and
*      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
*      of "The Art of Computer Programming", Addison-Wesley, (1981).
*
*      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
*      When Type is 1, returns a random number.
*      When Type is 2, gets the status of the generator.
*      When Type is 3, restores the status of the generator.
*
*      The status of the generator is represented by Status[0..56].
*      Make sure you initialize the seed before you get random
*      numbers.
****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double RandomGen(char Type, long Seed, long *Status)
{
	static long i1, i2, ma[56];   /* ma[0] is not used. */
	long        mj, mk;
	short       i, ii;

	if (Type == 0) { /* set seed. */
		mj = MSEED - (Seed < 0 ? -Seed : Seed);
		mj %= MBIG;
		ma[55] = mj;
		mk = 1;
		for (i = 1; i <= 54; i++) {
			ii = (21 * i) % 55;
			ma[ii] = mk;
			mk = mj - mk;
			if (mk < MZ)
				mk += MBIG;
			mj = ma[ii];
		}
		for (ii = 1; ii <= 4; ii++)
			for (i = 1; i <= 55; i++) {
				ma[i] -= ma[1 + (i + 30) % 55];
				if (ma[i] < MZ)
					ma[i] += MBIG;
			}
		i1 = 0;
		i2 = 31;
	}
	else if (Type == 1) {       /* get a number. */
		if (++i1 == 56)
			i1 = 1;
		if (++i2 == 56)
			i2 = 1;
		mj = ma[i1] - ma[i2];
		if (mj < MZ)
			mj += MBIG;
		ma[i1] = mj;
		return (mj * FAC);
	}
	else if (Type == 2) {       /* get status. */
		for (i = 0; i < 55; i++)
			Status[i] = ma[i + 1];
		Status[55] = i1;
		Status[56] = i2;
	}
	else if (Type == 3) {       /* restore status. */
		for (i = 0; i < 55; i++)
			ma[i + 1] = Status[i];
		i1 = Status[55];
		i2 = Status[56];
	}
	else
		printf("Wrong parameter to RandomGen().");
	return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/*******  end subroutine  ******/

//UNIFORM RANDOM SETUP C LIB
std::default_random_engine random(time(NULL));
std::uniform_real_distribution<double> dis2(0.0, 1.0);
double RandomNum()
{
	double x;

	//random type 1
	//x = drand48();

	//random type 2
	x = dis2(random);

	//random type3
	//x = RandomGen(1, 0, NULL);

	if (x <= 1e-7 || x >= (1 - 1e-7)) x = RandomNum();

	return x;
}
#endif