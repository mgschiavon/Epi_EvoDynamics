//
// Evolutionary dynamics of epigenetic switches in fluctuating environments
// SIMULATIONS
//
// Created by Mariana Gómez-Schiavon
// October 2015
// Modified April 2019
//

// This code simulates the evolution of a population with the gene network
// defined in the GeneNetwork class (see GeneNetwork.h). The population size
// (N), environmental fluctuation frequency (nu), selection pressure (sT),
// mutation frequency (u), mutation step size (M), gene expression dynamics'
// algorithm (alg: 0, Deterministic; 1, Gillespie), as well as the initial
// genotype (k0; nH0; KD0) are input arguments defined by the user. The
// simulation runs for GMAX generations, and the gene expression dynamics for
// each individual cell are simulated for TMAX time, which corresponds to its
// life span. The same initial protein number (A0) is assumed for all cells in
// the population. Some parameters in the gene network model are assumed to be
// fixed: the basal synthesis activity (ALPHA) and the protein degradation rate
// (DEG). From the three parameters that define the cell's genotype (k; nH; KD),
// which ones would actually allow to evolve is determined by the parE class
// (PAR_E). The environment fluctuates with a fixed frequency (nu) between two
// states (0, High; 1, Low), and the initial environmental state is defined by
// E0. Each environment is characterized by an optimal protein number (AOPT);
// and the fitness in each environment is calculated using a Lorentzian function
// centered in the optimal protein number and a width measure V. The cells to be
// cloned in the next generations are selected applying the tournament selection
// scheme, with a probability u of randomly mutate each of its parameter values.
// Every generation, the average genotype & subpopulations is printed in the AGS
// output file; every cycle, the adaptation strategy per cycle is printed in the
// ASC output file; and the lineage information per cycle is printed only in the
// last Cp cycles in theis printed in the LxC output file; and finally, the
// individuals information per epoch for the first nine generation after each
// environmental change and the last generation of each epoch is printed for
// the last Ep epochs is printed in the IxE output file.
// NOTE: A variant from EvoDynamics.cpp where the seed index needs to be 
//       specified as input by the user (from 0 to 9).

// Libraries:
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <cmath>
#include "GeneNetwork.h"

int main(int argn, char *args[])
{
	///////////////////////////////////////////////////////////////////////////////////////
	// Simulation constants:
	#define TMAX 4			// Maximum time to simulate each cell during one generation.
	int GMAX = 10000;		// Number of generations to simulate.
	int Cp = 100;			// Number of cycles to print full lineages' information.
	int Ep = 10;			// Number of epoch to print full individuals' information (only first 9 & last generation per epoch).
	// Initial protein concentration:
	int A0 = 80;
	// Initial environmental state:
	int E0 = 1;			// 0 - High; 1 - Low
	// Fixed parameters:
	double ALPHA = 0.25;		// Basal synthesis activity (a=1 -> No feedback)
	double DEG = 1;			// Protein degradation constant
	// Evolving parameters (1 - Evolve; 0 - Fixed):
	parE PAR_E;
	PAR_E.k = 1;			// Maximum synthesis rate
	PAR_E.nH = 1;			// Hill coefficient
	PAR_E.KD = 1;			// Affinity constant
	// Optimal protein number:
	double AOPT[2] = {80,20}; 	// {E_High,E_Low}
	double V = 0.2;			// Width measure for Lorentzian fitness function.
	// Replicas' seeds:
	long seeds[10] = {-17,-23,-7,-3,-5,-9,-11,-13,-15,-19};	// Seeds for the random number generator.
	///////////////////////////////////////////////////////////////////////////////////////

	if (argn < 11)
	{
		cerr << "Error! Input ten arguments: N; nu; sT; u; M; Algorithm: (0) Deterministic, (1) Gillespie; k0; nH0; KD0; seed index (0 to 9)." << endl;
		exit(0);
	}

	// Input arguments
	int N = atoi(args[1]);		// Number of cells in the population.
	double nu = atof(args[2]);	// Frequency of environment switching.
	int sT = atoi(args[3]);		// Selection strength.
	double u = atof(args[4]);	// Mutation rate (x individual x generation).
	double M = atof(args[5]);	// Mutation factor.
	bool alg = atoi(args[6]);	// Type of algorithm to simulate with (0 = Deterministic, 1 = Gillespie).
	double k0 = atof(args[7]);	// Initial value for the maximum synthesis rate (k).
	double nH0 = atof(args[8]);	// Initial value for the Hill coefficient (nH).
	double KD0 = atof(args[9]);	// Initial value for the affinity constant (KD).
	int iS = atof(args[10]);	// Replica index (i.e. value from seeds[] to use to generate random numbers).

		// OUTPUT FILES
		char myFileName[255];
		// Average genotype & subpopulations:
		sprintf(myFileName,"AGS_N%d_nu%1.4f_s%d_u%1.3f_M%1.2f_a%d_k%1.3f_nH%1.2f_KD%1.2f_%d.dat",N,nu,sT,u,M,alg,k0,nH0,KD0,iS);
		ofstream AGS(myFileName,ios::out);
		AGS << "Generation" << '\t' << "Environment" << '\t' << "Bis-Fraction" << '\t';
		AGS << "<k>_B" << '\t' << "<nH>_B" << '\t' << "<KD>_B" << '\t';
		AGS << "<A>_B" << '\t' << "<w>_B" << '\t';
		AGS << "<k>_M" << '\t' << "<nH>_M" << '\t' << "<KD>_M" << '\t';
		AGS << "<A>_M" << '\t' << "<w>_M" << endl;
		// Adaptation strategy x cycle:
		sprintf(myFileName,"ASC_N%d_nu%1.4f_s%d_u%1.3f_M%1.2f_a%d_k%1.3f_nH%1.2f_KD%1.2f_%d.dat",N,nu,sT,u,M,alg,k0,nH0,KD0,iS);
		ofstream ASC(myFileName,ios::out);
		ASC << "Generation" << '\t';
		ASC << "PE" << '\t';
		ASC << "IE" << '\t';
		ASC << "BA" << '\t';
		ASC << "GA" << endl;
		// Lineage information x cycle:
		sprintf(myFileName,"LxC_N%d_nu%1.4f_s%d_u%1.3f_M%1.2f_a%d_k%1.3f_nH%1.2f_KD%1.2f_%d.dat",N,nu,sT,u,M,alg,k0,nH0,KD0,iS);
		ofstream LxC(myFileName,ios::out);
		LxC << "Generation" << '\t';
		LxC << "LTag" << '\t' << "LMut" << '\t' << "LBis" << '\t' << "LFit" << '\t';
		LxC << "k_H" << '\t' << "nH_H" << '\t' << "KD_H" << '\t';
		LxC << "k_L" << '\t' << "nH_L" << '\t' << "KD_L" << endl;
		// Individuals information per epoch (first 9 & last generations):
		sprintf(myFileName,"IxE_N%d_nu%1.4f_s%d_u%1.3f_M%1.2f_a%d_k%1.3f_nH%1.2f_KD%1.2f_%d.dat",N,nu,sT,u,M,alg,k0,nH0,KD0,iS);
		ofstream IxE(myFileName,ios::out);
		IxE << "Generation" << '\t' << "Environment" << '\t' << "Cell tag" << '\t';
		IxE << "k" << '\t' << "nH" << '\t' << "KD" << '\t';
		IxE << "A" << '\t' << "w" << '\t';
		IxE << "SS1" << '\t' << "SS2" << '\t' << "SS3" << '\t';
		IxE << "Parent tag" << endl;

		long seed = seeds[iS];	// Initialize random seed.
		int E = 1 - E0;			// Initialize environment.

		// Create the cell population:
		GeneNetwork myPop[N];
		PopAverage myPopB;
		PopAverage myPopM;

		// Assign user-defined-parameters values:
		myPop[0].ic(A0,k0,ALPHA,nH0,KD0,DEG);
		myPop[0].LTag = 0;
		for (int i=1; i<N; i++)
		{
			myPop[i] = myPop[0];
			myPop[i].LTag = i;
		}

		// Simulate generations:
		for (int iG=1; iG<=GMAX; iG++)
		{
			// Update environmental state each epoch:
			if(iG%int(1/nu)==1)
			{
				E = (E+1)%2;
			}

			myPopB.reset();
			myPopM.reset();

			// Loop over the population:
			for (int i = 0; i<N; i++)
			{
				// Find the solutions of the system (SS) -only if a mutation has occurred:
				if(myPop[i].mut)
				{
					myPop[i].findSS();
				}

				// Gillespie Algorithm:
				if (alg)
				{
					while (myPop[i].time < TMAX)
					{
						myPop[i].updateTime(&seed);
						myPop[i].updateSystem(&seed);
					}
				}
				// Deterministic Algorithm (assuming steady state is reached in the life span):
				else
				{
					// If the cell is bistable & the cell was previously in the "high" state, stay in this state:
					if((myPop[i].SS[2] > 0) && (myPop[i].SS[1] < myPop[i].A))
					{
						myPop[i].A = myPop[i].SS[2];
					}
					else
					{
						myPop[i].A = myPop[i].SS[0];
					}
				}

				// Calculate the fitness of cell [i]:
				myPop[i].calculateFitness(E, AOPT, V);
				myPop[i].LFit += (myPop[i].w*(nu/2));

				// MONOSTABLE
				if(myPop[i].SS[2] == 0)
				{
					myPop[i].LBis = 0;
					myPopM.k += myPop[i].k;
					myPopM.nH += myPop[i].nH;
					myPopM.KD += myPop[i].KD;
					myPopM.A += myPop[i].A;
					myPopM.w += myPop[i].w;
					myPopM.f++;
				}
				// BISTABLE
				else
				{
					myPopB.k += myPop[i].k;
					myPopB.nH += myPop[i].nH;
					myPopB.KD += myPop[i].KD;
					myPopB.A += myPop[i].A;
					myPopB.w += myPop[i].w;
					myPopB.f++;
				}

				// PRINTING - Individuals information per epoch (first 9 & last generations):
				if((iG>(GMAX-(Ep/nu))) && (iG%int(1/nu)<10))
				{
					IxE << iG << '\t' << E << '\t' << i << '\t';
					IxE << myPop[i].k << '\t' << myPop[i].nH << '\t' << myPop[i].KD << '\t';
					IxE << myPop[i].A << '\t' << myPop[i].w << '\t';
					IxE << myPop[i].SS[0] << '\t' << myPop[i].SS[1] << '\t' << myPop[i].SS[2] << '\t';
					IxE << myPop[i].parent << endl;
				}
			}

			// PRINTING - Average genotype & subpopulations:
			myPopB.normalize(N);
			myPopM.normalize(N);
			AGS << iG << '\t' << E << '\t' << myPopB.f << '\t';
			AGS << myPopB.k << '\t' << myPopB.nH << '\t' << myPopB.KD << '\t';
			AGS << myPopB.A << '\t' << myPopB.w << '\t';
			AGS << myPopM.k << '\t' << myPopM.nH << '\t' << myPopM.KD << '\t';
			AGS << myPopM.A << '\t' << myPopM.w << endl;

			// Each epoch, record genotype per lineage:
			if(iG%int(1/nu)==0)
			{
				for (int i=0; i<N; i++)
				{
					myPop[i].kE[E] = myPop[i].k;
					myPop[i].nHE[E] = myPop[i].nH;
					myPop[i].KDE[E] = myPop[i].KD;
				}
			}


			// Each cycle:
			if(iG%int(2/nu)==0)
			{
				double AS[5] = {0,0,0,0,0};
				for (int i=0; i<N; i++)
				{
					myPop[i].adaptS();
					AS[myPop[i].ASp]++;

					if(iG>(GMAX-(Cp*2/nu)))
					{
						// PRINTING - Lineage information x cycle:
						LxC << iG << '\t';
						LxC << myPop[i].LTag << '\t' << myPop[i].LMut << '\t' << myPop[i].LBis << '\t' << myPop[i].LFit << '\t';
						LxC << myPop[i].kE[0] << '\t' << myPop[i].nHE[0] << '\t' << myPop[i].KDE[0] << '\t';
						LxC << myPop[i].kE[1] << '\t' << myPop[i].nHE[1] << '\t' << myPop[i].KDE[1] << endl;
					}
					// Reset lineage's statistics:
					myPop[i].LMut = 0;
					myPop[i].LBis = 1;
					myPop[i].LFit = 0;
					myPop[i].LTag = i;
				}
				// PRINTING - Adaptation strategy x cycle:
				ASC << iG << '\t';
				ASC << AS[1]/N << '\t';
				ASC << AS[2]/N << '\t';
				ASC << AS[3]/N << '\t';
				ASC << AS[4]/N << endl;
			}

			// Apply selection and evolve parameters:
			selectionTournament(myPop, N, sT, u, M, &seed, PAR_E);
		}

		// Close files
		AGS.close();
		ASC.close();
		LxC.close();
		IxE.close();

	return 0;
}
