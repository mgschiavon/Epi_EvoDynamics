# Epi_EvoDynamics
Evolutionary dynamics of an epigenetic switch in a fluctuating environment

To simulate the evolutionary dynamics of a self-regulated gene as described in [Gómez-Schiavon & Buchler](https://doi.org/10.1101/072199), three main files are required:
- `EvoDynamics.cpp` : Main instructions
- `GeneNetwork.h`   : Library
- `ran1.h`          : Random number generator

All of them in C++ language. The "Minimal random number generator of Park and Miller with Bays-Durham shuffle and added safeguards" from Numerical Recipes in C++ [1] was used to generate random numbers. A negative integer is required to initialize (i.e. *seed*).

[1] *Press, W. H., Teukolsky, S. A., and Vetterling, W. T. (2002), Numerical recipes in C++: the art of scientific computing, Cambridge University Press.*

The code `EvoDynamics.cpp` simulates the evolution of a population with the gene network defined in the `GeneNetwork` class (see `GeneNetwork.h`). The population size (`N`), environmental fluctuation frequency (`nu`), selection pressure (`sT`), mutation frequency (`u`), mutation step size (`M`), gene expression dynamics' algorithm (`alg`: `0`, *deterministic*; `1`, *stochastic/Gillespie*), as well as the initial genotype (`k0`; `nH0`; `KD0`) are input arguments defined by the user. The simulation runs for `GMAX` generations, and the gene expression dynamics for each individual cell are simulated for `TMAX` time, which corresponds to its life span. The same initial protein number (`A0`) is assumed for all cells in the population. Some parameters in the gene network model are assumed to be fixed: the basal synthesis activity (`ALPHA`) and the protein degradation rate (`DEG`). From the three parameters that define the cell's genotype (*k*; *n<sub>H</sub>*; *K<sub>D</sub>*), which ones are allowed to evolve is determined by the `parE` class (`PAR_E`). The environment fluctuates with a fixed frequency (`nu`) between two states (`0`, *HIGH*; `1`, *LOW*), and the initial environmental state is defined by `E0`. Each environment is characterized by an optimal protein number (`AOPT`); and the fitness in each environment is calculated using a Lorentzian function centered in the optimal protein number and a width measure `V`. The cells to be cloned in the next generations are selected applying the tournament selection scheme, with a probability `u` of randomly mutate each of its parameter values. Every generation, the average genotype in the subpopulations is printed in the `AGS` output file; every cycle, the adaptation strategy per cycle is printed in the `ASC` output file; and the lineage information per cycle is printed only in the last `Cp` cycles in these printed in the `LxC` output file; and finally, the individuals information per epoch for the first nine generation after each environmental change and the last generation of each epoch is printed for the last `Ep` epochs is printed in the `IxE` output file.

## Running instructions
1. Compile `EvoDynamics.cpp` in the presence of the other two files, for example:
```
g++ EvoDynamics.cpp -o ED_run.out
```

2. Run with the input arguments: population size (`N`), environmental fluctuation frequency (`nu`), selection pressure (i.e. tournament size, `sT`), mutation rate (`u`), mutation step size (`M`), algorithm to use (`0`, if *deterministic*; `1`, if *stochastic*), and the initial genotype (`k0`; `nH0`; `KD0`). For example:
```
./ED_run.out 1000 0.1 6 0.03 1.1 1 80 6 45
```

In addition of the input arguments, other simulation parameters can be easily modified by the user in the `EvoDynamics.cpp` file:
- `TMAX`: Maximum time to simulate each cell during one generation (e.g. `4`).
- `GMAX` : Number of generations to simulate (e.g. `10000`).
- `Cp` : Number of cycles to print full lineages' information (e.g. `100`).
- `Ep` : Number of epoch to print full individuals' information (only first 9 & last generation per epoch; e.g. `10`).
- `A0` : Initial protein concentration (e.g. `80`).
- `E0` : Initial environmental state (`0` if *HIGH*; `1` if *LOW*; e.g. `1`).
- `ALPHA` : Basal synthesis activity (`a=1` -> No feedback; e.g. `0.25`).
- `DEG` : Protein degradation constant (e.g. `1`).
- `PAR_E.k` : Maximum synthesis rate (`1` if *Evolve*; `0` if *Fixed*; e.g. `1`).
- `PAR_E.nH` : Hill coefficient (`1` if *Evolve*; `0` if *Fixed*; e.g. `1`).
- `PAR_E.KD` : Affinity constant (`1` if *Evolve*; `0` if *Fixed*; e.g. `1`).
- `AOPT[2]` : Optimal protein number ({A<sub>opt</sub><sup>(H)</sup>,A<sub>opt</sub><sup>(L)</sup>; e.g. `{80,20}`).
- `V` : Width measure for Lorentzian fitness function (e.g. `0.2`).
- `numRep` : Total number of replicas to run (e.g. `10`).
- `seeds[10]` : Seeds for the random number generator (e.g. `{-17, -23, -7, -3, -5, -9, -11, -13, -15, -19}`).

### Output files

The output files for each replica are:

- Average genotype & subpopulations: File `AGS` with columns (printed each generation):
  - Generation *g*
  - Environment *E* (`0` if *HIGH*; `1` if *LOW*)
  - Bistable fraction *f<sub>B</sub>*
  - *\<k><sub>B</sub>*
  - *<n<sub>H</sub>><sub>B</sub>*
  - *<K<sub>D</sub>><sub>B</sub>*
  - *\<A><sub>B</sub>*
  - *\<w><sub>B</sub>*
  - *\<k><sub>M</sub>*
  - *<n<sub>H</sub>><sub>M</sub>*
  - *<K<sub>D</sub>><sub>M</sub>*
  - *\<A><sub>M</sub>*
  - *\<w><sub>M</sub>*

NOTE: *B* and *M* subindexes stand for *bistable* and *monostable* subpopulations, respectively.
   
- Adaptation strategy per cycle: File `ASC` with columns (printed at the end of each environmental cycle):
  - Generation *g*
  - PE
  - IE
  - BA
  - GA

NOTE: *PE* corresponds to the fraction of parental lineages with completely bistable lineages and no mutations; *IE* corresponds to the fraction of parental lineages with completely bistable lineages and exactly one mutation; *BA* corresponds to the fraction of parental lineages with completely bistable lineages and 2 or more mutations; and *GA* correspond to the fraction of parental lineages with some monostable genotypes and 2 or more mutations.
    
- Lineage information per cycle: File `LxC` with columns (printed at the end of the last `Cp` environmental cycles):
  - Generation *g*
  - Tag of parental lineage
  - Number of accumulated mutations
  - Lineage bistability (`1` if all genotypes were bistable, `0` otherwise)
  - Average fitness of the lineage in the cycle
  - *k<sup>(H)</sup>*
  - *n<sub>H</sub><sup>(H)</sup>*
  - *K<sub>D</sub><sup>(H)</sup>*
  - *k<sup>(L)</sup>*
  - *n<sub>H</sub><sup>(L)</sup>*
  - *K<sub>D</sub><sup>(L)</sup>*

NOTE: *x<sup>(H)</sup>*, *k<sup>(L)</sup>* correspond to the biophysical parameter *x* at the end of the *HIGH* and *LOW* epochs, respectively. Notice that these are not parental lineages, but simply the lineages of the current cycle.
    
- Individuals information per epoch (first 9 \& last generations): File `LxC` with columns (printed at the end of the last `Ep` epochs):
  - Generation *g*
  - Environment *E* (`0` if *HIGH*; `1` if *LOW*)
  - Cell tag
  - *k*
  - *n<sub>H</sub>*
  - *K<sub>D</sub>*
  - *A*
  - *w*
  - SS1
  - SS2
  - SS3
  - Parent tag

NOTE: SS\# correspond to the phenotype steady state values given the cell's genotype (if the cell is monostable, SS2 = SS3 = 0). The cell tag simply corresponds to the numbered population, from 0 to N-1, and parent tag tells you the corresponding cell tag of the parental cell in the previous generation.

In all cases, columns are separated by tabs.

### Alternative assumptions

Substitute main instructions file (`EvoDynamics.cpp`) and/or library accordingly (`GeneNetwork.h`):

#### Alternative selection schemes:
- `EvoDynamics_ALT_SS_Trunc.cpp` : Main instructions for Truncation selection
- `EvoDynamics_ALT_SS_Prop.cpp` : Main instructions for Proportional selection
- `EvoDynamics_ALT_SS_Weight.cpp` : Main instructions for Weighted selection
- `GeneNetwork_ALT_SS.h` : Library for all alternative selection schemes

#### Random environment fluctuations:
- `EvoDynamics_ALT_RandEnv.cpp` : Main instructions

#### Moran model:
- `EvoDynamics_ALT_MoranM.cpp` : Main instructions
- `GeneNetwork_ALT_MoranM.h` : Library

#### Alternative fitness functions:
- `EvoDynamics_ALT_WFunct_Normal.cpp` : Normal Fitness function
- `EvoDynamics_ALT_WFunct_Step.cpp` : Step-like Fitness function
- `GeneNetwork_ALT_WFunct.h` : Library

#### Average A as phenotype:
- `EvoDynamics_ALT_Phen_AP.cpp` : Main instructions
- `GeneNetwork_ALT_Phen_AP.h` : Library

#### Protein distribution as phenotype:
- `EvoDynamics_ALT_Phen_PD.cpp` : Main instructions
- `GeneNetwork_ALT_Phen_PD.h` : Library

#### "Cube" mutation:
- `EvoDynamics_ALT_C3Mut.cpp` : Main instructions
- `GeneNetwork_ALT_C3Mut.h` : Library

#### Basal activity (alpha) evolving:
- `EvoDynamics_aEvo.cpp` : Main instructions
- `GeneNetwork_aEvo.h` : Library

## Referencing

If you use this code or the data associated with it please cite:

Gómez-Schiavon & Buchler (2017); https://doi.org/10.1101/072199.

## Copyright

(C) Copyright 2018 Mariana Gómez-Schiavon

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
