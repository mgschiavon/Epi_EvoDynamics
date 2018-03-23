//
//	Evolutionary dynamics of epigenetic switches in fluctuating environments
//	EVOLUTIONARY MODEL
//
// Created by Mariana Gómez-Schiavon
// October 2015
//

// Libraries:
#include <stdio.h>
#include <cmath>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include "ran1.h"

using namespace std;

class parE
{
	public:
		bool k;
		bool a;
		bool nH;
		bool KD;
		bool d;
	
	parE()
	{
		k = 0;
		a = 0;
		nH = 0;
		KD = 0;
		d = 0;
	}
};

class PopAverage
{
	public:
		double k, nH, KD;			// Evolving parameters
		double A;					// Protein number
		double w;					// Fitness
		double f;					// Fraction
	
	void reset()
	{
		k = 0;
		nH = 0;
		KD = 0;
		A = 0;
		w = 0;
		f = 0;
	}
	
	void normalize(int N)
	{
		k = k/f;
		nH = nH/f;
		KD = KD/f;
		A = A/f;
		w = w/f;
		f = f/N;
	}
};

class GeneNetwork
{
	public:
		// Genotype:
		double k, nH, KD;			// Evolving parameters
		double a, d;				// Fixed parameters
		// Phenotype:
		double A;					// Protein number
		double SS[3];				// Protein number in the steady state(s)
		// Evolutionary metrics:
		double w;					// Fitness
		int parent;					// Parent cell number
		bool mut;					// Flag for new mutant individuals (0 - Parental genotype; 1 - Mutant)
		// Lineage metrics:
		int LTag;					// Tag for parental lineage
		int LMut;					// Number of accumulated mutations per cycle
		bool LBis;					// Flag for bistable lineages (1 - Bistable genotype the whole cycle; 0 - Otherwise)
		double LFit;				// Average fitness of the lineage on the cycle
		double kE[2];				// Genotype (k) at the end of each epoch [HIGH,LOW]
		double nHE[2];				// Genotype (nH) at the end of each epoch [HIGH,LOW]
		double KDE[2];				// Genotype (KD) at the end of each epoch [HIGH,LOW]
		int AS;						// Flag for the adaptation strategy used by the lineage (0 - Ma; 1 - PE; 2 - IE; 3 - BA; 4 - GA)
		int ASp;					// Flag for the adaptation strategy used by the parental lineage (0 - Ma; 1 - PE; 2 - IE; 3 - BA; 4 - GA)
		// Gene dynamics variables:
		double r[2], rT;		// Reaction propensities & total reaction propensity
		double time;			// Running time
	
	GeneNetwork()
	{
		w = 0;
		parent = 0;
		mut = 0;
		
		LMut = 0;
		LBis = 1;
		LFit = 0;
		AS = 0;
		
		time = 0;
	}
	
	// Assign initial genotype & phenotype:
	void ic(double A0, double k0, double a0, double nH0, double KD0, double d0)
	{
		A = A0;		// Protein number
		k = k0;		// Maximum synthesis rate
		a = a0;		// Basal synthesis activity (a=1 -> No feedback)
		nH = nH0;	// Hill coefficient
		KD = KD0;	// Affinity constant
		d = d0;		// Protein degradation constant
		findSS();	// Calculate steady states (SS)
	}
	
	// Find the steady state solutions:
	void findSS()
	{
		double Ai = 0;
		int ss = 0;
		double FA0, FA;
		
		// Initialize SS:
		SS[0] = 0;
		SS[1] = 0;
		SS[2] = 0;
		
		// Initial value of the equilibrium equation:
		FA0 = Ai - (k/d)*(a+((1-a)*(pow(Ai,nH)/(pow(KD,nH)+pow(Ai,nH)))));
		
		for (Ai=0.1; Ai <= 1000; Ai=Ai+0.1)
		{
			// Calculate equilibrium equation:
			FA = Ai - (k/d)*(a+((1-a)*(pow(Ai,nH)/(pow(KD,nH)+pow(Ai,nH)))));
			// If there is a change of sign in the equilibrium equation, there should be a steady state:
			if (FA*FA0 <= 0)
			{
				SS[ss] = Ai;
				ss += 1;
			}
			// It stops if the system already have 3 steady states:
			if (ss == 3)
			{
				break;
			}
			// Make the actual value the previous value and iterate:
			FA0 = FA;
		}
	}

	// Update reaction propensities:
	void updateR()
	{
		r[0] = k*(a+((1-a)*(pow(A,nH)/(pow(KD,nH)+pow(A,nH)))));	// A synthesis reaction
		r[1] = d*A;													// A degradation reaction
		rT = r[0]+r[1];												// Total propensity
	}
		
	// Perform the Gillespie algorithm:
	void updateTime(long *seed)
	{
		double rand1 = (double)ran1(seed);		// Random number between 0 and 1
		double re_time;							// Time between reactions
		
		updateR();								// Update reaction propensities
		re_time = -log(1-rand1)/rT;				// Time for the next reaction
		time += re_time;						// Actualize global time
	}
	void updateSystem(long *seed)
	{ 
		double rand2 = (double)ran1(seed);		// Random number between 0 and 1
		if (rand2 < r[0]/rT)
		{
			A++;								// Synthesis reaction
		}
		else if (rand2 <= (r[0]+r[1])/rT)
		{
			A--;								// Protein degradation
		}
		else
		{
			cerr << "Runtime Error: rT not real." << endl;
			exit(0);
		}
	}
	
	// Evaluate cell's fitness:
	void calculateFitnessNORMAL(int myE, double Aopt[2], double v)
	{
		w = exp(-(A-Aopt[myE])*(A-Aopt[myE])/(2*(v*Aopt[myE])));	// Normal function.
	}
	void calculateFitnessSTEP(int myE, double Aopt[2], double v)
	{
		w = 0;
		if (abs(A-Aopt[myE]) <= sqrt(v*Aopt[myE]))		// Step-like function around the optimum.
		{
			w = 1;
		}
	}
	
	// Random mutation in the genotype:
	void mutate(long *seed, double M, parE PAR_E)
	{
		#define PI 3.14159265
		double ratioQ = (double)ran1(seed);			// Random number between 0 and 1
		double thetaQ = (double)ran1(seed)*2*PI;	// Random number between 0 and 2*PI
		double phiQ = (double)ran1(seed)*2*PI;		// Random number between 0 and 2*PI
		
		if(PAR_E.k)
		{
			k *= pow(M,ratioQ*sin(thetaQ)*cos(phiQ));
			if(k>1000)				// Defines an upper limit for the parameter value
			{
				k=1000;
			}
			if(k<0.01)		// Defines a lower limit for the parameter value
			{
				k=0.01;
			}
		}
		if(PAR_E.nH)
		{
			nH *= pow(M,ratioQ*sin(thetaQ)*sin(phiQ));
			if(nH>16)			// Defines an upper limit for the parameter value
			{
				nH=16;
			}
			if(nH<0.01)		// Defines a lower limit for the parameter value
			{
				nH=0.01;
			}
		}
		if(PAR_E.KD)
		{
			KD *= pow(M,ratioQ*cos(thetaQ));
			if(KD>120)			// Defines an upper limit for the parameter value
			{
				KD = 120;
			}
			if(KD<0.01)		// Defines a lower limit for the parameter value
			{
				KD=0.01;
			}
		}
		mut = 1;
		LMut++;
	}
	
	// Define lineage adaptation strategy:
	void adaptS()
	{
		ASp = AS;
		if(LBis==1)
		{
			if(LMut==0)
			{
				AS = 1;			// Pure epigenetic: Bis && No mut
			}
			else if(LMut==1)
			{
				AS = 2;			// Impure epigenetic: Bis && 1 mut
			}
			else
			{
				AS = 3;			// Bistable adaptation: Bis && >1 mut
			}
		}
		else
		{
			if(LMut>1)
			{
				AS = 4;			// Genetic adaptation: No Bis && >1 mut
			}
			else
			{
				AS = 0;			// Maladapted: No Bis && <2 mut
			}
		}
	}
};

// TOURNAMENT SELECTION
void selectionTournament(GeneNetwork *myPop, int N, int t, double u, double M, long *seed, parE PAR_E)
{
	GeneNetwork newPop[N];	// The new population to be defined.
	
	// Fill in the next generation:
	int iN = 0;
	int rand1;
	while (iN < N)
	{
		// Stochastically choose the tournament cells and replicate the cell with higher fitness in the group:
		int indexT[t];
		indexT[0] = int(N*(double)ran1(seed));
		int maxT = indexT[0];
		int j = 1;
		while (j < t)
		{
			rand1 = int(N*(double)ran1(seed));
			for (int i=0; i < j; i++)
			{
				if (rand1 == indexT[i])
				{
					rand1 = -1;
					break;
				}
			}
			if (rand1 > 0)
			{
				indexT[j] = rand1;
				if(myPop[maxT].w < myPop[indexT[j]].w)
				{
					maxT = indexT[j];
				}
				j++;
			}
		}
		newPop[iN] = myPop[maxT];	// The new cell is equal to the old cell selected.
		newPop[iN].parent = maxT;	// Parent cell number
		iN++;
	}
	
	// Refill the population matrix & allow random mutation:
	for (int i = 0; i < N; i++)
	{
		myPop[i] = newPop[i];
		myPop[i].time = 0;		// Reset the time to zero.
		myPop[i].mut = 0;		// Flag for parental genotype.
		if ((double)ran1(seed) <= u)
		{
			myPop[i].mutate(seed, M, PAR_E);
		}
	}
}
