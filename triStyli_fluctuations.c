#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>
#define VERSION "27.06.2019"
#define DEPENDENCY "None\n"
#define MAX_NUMBER_OF_INITIAL_NTRL_ALLELES 999	// number of segregating alleles when generating the first parental population
#define RANGE 0.1	// value in [0;1] to modify the current allelic effect between [(1-RANGE) x current_value ; (1+RANGE) * current_value].
#define KRED  "\033[1m\033[31m"
#define KMAG  "\x1B[31m"
#define STOP  "\x1B[0m"

//	gcc triStyli_fluctuations.c -L/usr/local/lib -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o triStyli

// choice of mating partners:
// 1) all hermaphrodite individuals share the same probability of being sampled as a mother
// 2) only a proportion 'p' of hermaphrodites is subsampled among all in order to be "putative futur fathers".
// A number K.p of available males is sampled using a Poisson distribution.
// These K.p available males are sampled without replacement, and have an individual fitness equal to 1-f_morphe_i (negative frequency selection)
// To produce N babies, N fathers are sampled among : i) the selected males (K.p) and ii) selected males with compatible phenotypes. All compatible males have the same fitness now.
typedef struct Deme Deme;
struct Deme{
	int nIndividus;
	long* ntrlLoci;	// contains the neutral alleles for nIndividus (number of individuals within the deme) x 2 (diploids) x nNtlrLoci (number of neutral loci)
	double* quantiLoci; // contains the allelic effects for nIndividus (number of individuals within the deme) x 2 (diploids) x nQuantiLoci (number of quantitative loci)
	int* sexChro; // two values per individual: 00 = XX (or ZZ), 01 = XY (or ZW)
	int* sex; // one value per individual: 0 = heterogametic, 1 = homogametic 
	double* femaleAllocation; // =sum of the allelic effects over the quantitative loci
	double* maleAllocation; // =(1 - femaleAllocation)
	int* nOffsprings; // Poisson random

	int* neutralS; // A neutral locus with initial conditions similar to those of locusS
	int* neutralM; // A neutral locus with initial conditions similar to those of locusM
	int* locusS; // locus S determining the morphe. Ss (10) or SS (11) --> short morph (0)
	int* locusM; // locus M determining the morphe. ssMm (00 10) or ssMM (00 11) --> mid-styled morph (1). ssmm (00 00) --> long-styled morph (2)
	int* morphe; // 0: short-styled, 1: mid-styled, 2: long-styled
};

void initializePopulation(gsl_rng* r, Deme* population, const int nDemes, const int nIndividuals[], const int nNtrlLoci, const int nQuantiLoci, const double fecundity, const int initialSituation);
void libererMemoirePopulation(Deme* population, const int nDemes);
void afficherPopulation(Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int sexualSystem);
void configMetapop(gsl_rng* r, Deme* population, const int nDemes, const int nIndividuals[], const double migration, const double extinction, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]);
void setToZero(const int nDemes, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]);
void setNumberIndividualsPerDeme(int nIndividuals[], const int nDemes, const int maxIndPerDem_1);
void setNumberIndividualsPerDeme_fluctuation(gsl_rng* r, const int nDemes, int nIndividuals[], const int maxIndPerDem_1, const int maxIndPerDem_2, const double P_transition_12_size, const double P_transition_21_size);
void setFMigColToZero(double f_mig_col[]);
void initializeNewPopulation(Deme* newPopulation, const int nDemes, const int nIndividuals[], const int nProducedSeeds[], const int nNtrlLoci, const int nQuantiLoci, const double fecundity);
void panmixie(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double ntrlMutation, const double quantiMutation, const double fecundity, const double reprodMale, const double sexAvantage, const int sexualSystem, const double selfingRate, const int currentGeneration, const int generationNewAllele, const int newAllele, const int generationNewAllele2, const int newAllele2, const int initialSituation, const double homomorphe_probability, const int freqDepSel, const int sysReprod);
void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials);
void weightedSample_without_replacement(gsl_rng* r, const int* liste, const double* weights, int* target, const int sizeOfListe, const int nTrials);
void replacement(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int nIndividuals[], const int nNtrlLoci, const int nQuantiLoci, const int nImmigrants[], int nProducedSeeds[], int extinctionStatus[], const int recolonization, const int generation, const int sexualSystem, const int colonizationModel, double* f_mig_col);
void writeNindividuals(const Deme* population, const int nDemes, const double extinction, const double migration, const int seed);
void genePop(Deme* population, const int nDemes, const int nNtrlLoci, const int seed, int time);
void checkCommandLine(int argc);
void statisticsPopulations(Deme* population, const int nDemes, const int maxIndPerDem_1, const int maxIndPerDem_2, const double P_transition_12_size, const double P_transition_21_size, const int nQuantiLoci, const double fecundity, const double migration, const double extinction, const int recolonization, const int sexualSystem, const double sexAvantage, const int seed, int time, const double selfingRate, const int colonizationModel, const double global_fst_cm, const double global_fst_coal, const double global_gpst, const double global_D, const double global_Fis, const double homomorphe_probability, const int sysReprod);
void statisticsMigrantsColonizers(const int seed, const double* f_mig_col);
double fstMullon(const int maxIndPerDem_1, const double extinction, const int recolonization, const double migration);
double fstRousset(const int maxIndPerDem_1, const double extinction, const int recolonization, const double migration, const int colonizationModel);
void sexInvador(gsl_rng* r, Deme* population, const int nDemes, const int* extinctionStatus, const double sexAvantage, const int sexualSystem, const double fecundity);
void global_stat(Deme* population, const int nDemes, const long nNtrlLoci, double* diff_stats);
void nc(const Deme* population, const int nDemes, const int nNtrlLoci, const int locus, double* target);
double heteroZ(const double* cont_table_tot);
void printMorphes(Deme* population, const int nDemes, const int generation);

int main(int argc, char *argv[]){
	checkCommandLine(argc); // stop the code if the number of arguments doesn't fit with the expected one

	int i = 0;
	int j = 0;

	// Get Parameters from comamnd line
	const int nDemes = atoi(argv[1]); // number of demes
	const int maxIndPerDem_1 = atoi(argv[2]);	// carrying capacity per deme
	const int maxIndPerDem_2 = atoi(argv[3]);	// carrying capacity per deme
	const double P_transition_12_size = atof(argv[4]);	// carrying capacity per deme
	const double P_transition_21_size = atof(argv[5]);	// carrying capacity per deme
	const int nGeneration = atoi(argv[6]);	// number of generations to simulate
	const int nGenerationExtinction = atoi(argv[7]);	// generations at which extinction+recolonization starts
	const int nNtrlLoci = atoi(argv[8]);	// number of neutral loci
	const double ntrlMutation = atof(argv[9]);	// mutation rate of the ntrl loci
	const int nQuantiLoci = atoi(argv[10]);	// number of quantitative loci
	const double quantiMutation = atof(argv[11]);	// mutation rate of the quantative loci
	const double fecundity = atof(argv[12]);	// max number of offspring when femAlloc=100%
	const double reprodMale = atof(argv[13]); // proportion of individuals reproducing through the male function
	const double migration = atof(argv[14]);	// immigration rate
	const double extinction = atof(argv[15]);	// extinction rate
	const int recolonization = atoi(argv[16]);	// number of recolonizing individuals
	const int colonizationModel = atoi(argv[17]);	// 0 = migrant pool model 1 = propagule pool model
	const int sexualSystem = atoi(argv[18]);    // 0 = only hermaphrodites; 1 = XY system; 2 = ZW system
	const double sexAvantage = atof(argv[19]); // avantage confered by the Y or Z chromosome over hermaphrodites
	const double selfingRate = atof(argv[20]); // probability to have an ovule being fertilized by sperm from the same individual
	const int verbose = atoi(argv[21]); // frequency at which statistics are written in the output file
	const int seed = atoi(argv[22]); // seed of the random generator
	
	const int initialSituation = atoi(argv[23]); // 0 = 1/3 of each morphes; 1 = ss mm; 2 = ss MM; 3 = SS mm; 4 = SS MM
	const int newAllele = atoi(argv[24]); // 0 = S; 1 = s; 2 = M; 3 = m (this argument is not considered if argument 19 is == 0)
	const int generationNewAllele = atoi(argv[25]); // generation at which ONE copy of the FIRST new allele is brought into the metapop (this argument is not considered if argument 19 is == 0)
	const int newAllele2 = atoi(argv[26]); // 0 = S; 1 = s; 2 = M; 3 = m (this argument is not considered if argument 19 is == 0)
	const int generationNewAllele2 = atoi(argv[27]); // generation at which ONE copy of the SECOND new allele is brought into the metapop (this argument is not considered if argument 19 is == 0)
	const double homomorphe_probability = atof(argv[28]); // probability of homomorphic matings in a polymorphic deme.
	const int freqDepSel = atoi(argv[29]); // if 0: no biased when sampling the reproductive males; if 1: reproductive males are sampled with a rare advantage
	const int sysReprod = atoi(argv[30]); // if 0: panmixia; if 1: trisStyle
	
	// Random generator
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, seed);

	// array of statistics summaryzing the genetic differentiation
	double* diff_stats = NULL;
	diff_stats = malloc(5 * sizeof(double));
	if( diff_stats == NULL){
		exit(0);
	}	

	// Initializing the metapopulation
	Deme* population = NULL;
	int* nIndividuals = NULL;
	population = malloc(nDemes * sizeof(Deme));
	nIndividuals = malloc(nDemes * sizeof(int));

	if(population == NULL || nIndividuals == NULL){
		exit(0);
	}

	setNumberIndividualsPerDeme(nIndividuals, nDemes, maxIndPerDem_1); // at generation 0 : all demes share the same carrying capacity maxIndPerDem_1
	initializePopulation(r, population, nDemes, nIndividuals, nNtrlLoci, nQuantiLoci, fecundity, initialSituation);

	// fst
	double global_fst_cm = 0.0;
	double global_fst_coal = 0.0;
	double global_gpst = 0.0;
	double global_D = 0.0;
	double global_Fis = 0.0;

	// Evolution of the metapopulation
	for(i=0; i<=nGeneration; i++){	// start of the loop 'i' over the generations
		//printf("\n%d\n", i);
		int* nImmigrants = NULL; // stock the number of received immigrants per deme
		int* extinctionStatus = NULL;	// per deme: 0 = non-extincted; 1 = extincted
		int* nProducedSeeds = NULL; // stock the numer of produced seeds per deme
		double* f_mig_col = NULL; // vector of size 6 containing: [f_mig_short, f_mig_mid, f_mig_long, f_col_short, f_col_mid, f_col_long]
		Deme* newPopulation = NULL;
		nImmigrants = malloc(nDemes * sizeof(int));
		extinctionStatus = malloc(nDemes * sizeof(int));
		nProducedSeeds = malloc(nDemes * sizeof(int));
		newPopulation = malloc(nDemes * sizeof(Deme));
		
		f_mig_col = malloc( 6 * sizeof(double));

		if(nImmigrants == NULL || extinctionStatus == NULL || nProducedSeeds == NULL || newPopulation == NULL || f_mig_col == NULL){
			exit(0);
		}
		setToZero(nDemes, nImmigrants, extinctionStatus, nProducedSeeds);	// initialize vectors used to configure the new-population
		setNumberIndividualsPerDeme_fluctuation(r, nDemes, nIndividuals, maxIndPerDem_1, maxIndPerDem_2, P_transition_12_size, P_transition_21_size);

//		for(j=0; j<nDemes; ++j){
			//printf("%d\t", nIndividuals[j]);
//		}
//		printf("\n");

		if(i >= nGenerationExtinction){ // extinction starts
			configMetapop(r, population, nDemes, nIndividuals, migration, extinction, nImmigrants, extinctionStatus, nProducedSeeds);	// get the parameters to create the new-population
		}else{
			configMetapop(r, population, nDemes, nIndividuals, migration, 0.0, nImmigrants, extinctionStatus, nProducedSeeds);	// get the parameters to create the new-population
		}
		
		initializeNewPopulation(newPopulation, nDemes, nIndividuals, nProducedSeeds, nNtrlLoci, nQuantiLoci, fecundity);


		panmixie(r, population, newPopulation, nDemes, nNtrlLoci, nQuantiLoci,  ntrlMutation, quantiMutation, fecundity, reprodMale, sexAvantage, sexualSystem, selfingRate, i, generationNewAllele, newAllele, generationNewAllele2, newAllele2, initialSituation, homomorphe_probability, freqDepSel, sysReprod);

		
		//if( i == nGeneration ){
		if( i%verbose == 0 ){
			for(j=0; j<5; j++){
				diff_stats[j] = 0.0; // initialize the vector of genPop statistics to zero
			}

			global_stat(population, nDemes, nNtrlLoci, diff_stats);
			global_fst_cm = diff_stats[0];
			global_fst_coal = diff_stats[1];
			global_gpst = diff_stats[2];
			global_D = diff_stats[3];
			global_Fis = diff_stats[4];
			
			if( i==0 ){
				statisticsPopulations(population, nDemes, maxIndPerDem_1, maxIndPerDem_2, P_transition_12_size, P_transition_21_size, nQuantiLoci, fecundity, migration, 0.0, recolonization, sexualSystem, sexAvantage, seed, i, selfingRate, colonizationModel, global_fst_cm, global_fst_coal, global_gpst, global_D, global_Fis, homomorphe_probability, sysReprod);
			}else{
				if(i >= nGenerationExtinction){ // extinction starts
					statisticsPopulations(newPopulation, nDemes, maxIndPerDem_1, maxIndPerDem_2, P_transition_12_size, P_transition_21_size, nQuantiLoci, fecundity, migration, extinction, recolonization, sexualSystem, sexAvantage, seed, i, selfingRate, colonizationModel, global_fst_cm, global_fst_coal, global_gpst, global_D, global_Fis, homomorphe_probability, sysReprod);
				}else{
					statisticsPopulations(newPopulation, nDemes, maxIndPerDem_1, maxIndPerDem_2, P_transition_12_size, P_transition_21_size, nQuantiLoci, fecundity, migration, 0.0, recolonization, sexualSystem, sexAvantage, seed, i, selfingRate, colonizationModel, global_fst_cm, global_fst_coal, global_gpst, global_D, global_Fis, homomorphe_probability, sysReprod);
				}
			}
//			genePop(newPopulation, nDemes, nNtrlLoci, seed, i);
		}
		// remplacer population par newPopulation
		// Nettoyer memoire
		libererMemoirePopulation(population, nDemes);
		free(population);
		population = malloc(nDemes * sizeof(Deme));

		if(population == NULL){
			exit(0);
		}
		
		setFMigColToZero(f_mig_col);
		replacement(r, population, newPopulation, nDemes, nIndividuals, nNtrlLoci, nQuantiLoci, nImmigrants, nProducedSeeds, extinctionStatus, recolonization, i, sexualSystem, colonizationModel, f_mig_col); // replace the parents (population) by the offspring (newPopulation)

		if( i%verbose == 0 ){
			statisticsMigrantsColonizers(seed, f_mig_col);
		}
	
//		if( i%verbose == 0 ){	
//		printMorphes(population, nDemes, i);
//		}
		free(nImmigrants);
		free(extinctionStatus);
		free(nProducedSeeds);
		libererMemoirePopulation(newPopulation, nDemes);
		free(newPopulation);
		free(f_mig_col);
	}	// end of the loop 'i' over the generations
	gsl_rng_free (r);
	libererMemoirePopulation(population, nDemes);
	free(population);
	free(nIndividuals);
	free(diff_stats);
	return(0);
}


void initializePopulation(gsl_rng* r, Deme* population, const int nDemes, const int nIndividuals[], const int nNtrlLoci, const int nQuantiLoci, const double fecundity, const int initialSituation){
	int i = 0;
	int j = 0;
	int k = 0;
	int valuesNtrlAlleles = MAX_NUMBER_OF_INITIAL_NTRL_ALLELES;
	double minQuanti = 0.0;
	double maxQuanti = 0.0;
	maxQuanti = 1.0/2/nQuantiLoci;
//	maxQuanti = 0.5/2/nQuantiLoci;	// uncomment to fix sex allocation to 0.5
	for(i=0; i<nDemes; i++){	// loop along the demes
		population[i].nIndividus = nIndividuals[i];
		population[i].ntrlLoci = malloc(2 * nIndividuals[i] * nNtrlLoci * sizeof(long));
		population[i].quantiLoci = malloc(2 * nIndividuals[i] * nQuantiLoci * sizeof(long));
		population[i].sexChro = malloc(2 * nIndividuals[i] * sizeof(int));
		population[i].sex = malloc(nIndividuals[i] * sizeof(int));
		population[i].femaleAllocation = malloc(nIndividuals[i] * sizeof(long));
		population[i].maleAllocation = malloc(nIndividuals[i] * sizeof(long));
		population[i].nOffsprings = malloc(nIndividuals[i] * sizeof(int));
	
		population[i].neutralS = malloc(2 * nIndividuals[i] * sizeof(int));  
		population[i].neutralM = malloc(2 * nIndividuals[i] * sizeof(int)); 
		population[i].locusS = malloc(2 * nIndividuals[i] * sizeof(int)); // locus S determining the morphe. Ss (10) or SS (11) --> short morph (0)
		population[i].locusM = malloc(2 * nIndividuals[i] * sizeof(int)); // locus M determining the morphe. ssMm (00 10) or ssMM (00 11) --> mid-styled morph (1). ssmm (00 00) --> long-styled morph (2)
		population[i].morphe = malloc(nIndividuals[i] * sizeof(int)); // 0: short-styled, 1: mid-styled, 2: long-styled
		
		if(population[i].ntrlLoci == NULL || population[i].quantiLoci == NULL || population[i].sexChro == NULL || population[i].sex == NULL || population[i].femaleAllocation == NULL || population[i].maleAllocation == NULL || population[i].nOffsprings == NULL || population[i].neutralS == NULL || population[i].neutralM == NULL || population[i].locusS == NULL || population[i].locusM == NULL || population[i].morphe == NULL){
			exit(0);
		}

		int cnt = -1;
		for(j=0; j<nIndividuals[i]; j++){	// loop along the individuals
			cnt += 1;
			population[i].femaleAllocation[j] = 0.0;
			population[i].maleAllocation[j] = 1.0;
			for(k=0; k<(2*nNtrlLoci); k++){	// loop along the {2: diploid} x {nNtrlLoci: number of neutral loci} positions
				population[i].ntrlLoci[j*2*nNtrlLoci+k] = gsl_rng_uniform_int(r, valuesNtrlAlleles) + 1;
			}
			for(k=0; k<(2*nQuantiLoci); k++){	// loop along the {2: diploid} x {nQuantiLoci: number of quantitative loci} positions
				population[i].quantiLoci[j*2*nQuantiLoci+k] = gsl_ran_flat(r, minQuanti, maxQuanti);
//				population[i].quantiLoci[j*2*nQuantiLoci+k] = maxQuanti; // uncomment to fix sex allocation to 0.5
				population[i].femaleAllocation[j] += population[i].quantiLoci[j*2*nQuantiLoci+k];
				population[i].maleAllocation[j] -= population[i].quantiLoci[j*2*nQuantiLoci+k];
			}

				population[i].sexChro[2*j] = 0;
				population[i].sexChro[2*j + 1] = 0;
				population[i].sex[j] = 1;
				
				// 0 = 1/3 of each morphes; 1 = ss mm; 2 = ss MM; 3 = SS mm; 4 = SS MM
				if( initialSituation == 0 ){
					if( j%3 == 0 ){ // short : all SsMm at generation 0
						population[i].neutralS[2*j] = 1;
						population[i].neutralS[2*j + 1] = 0;
						population[i].neutralM[2*j] = 1;
						population[i].neutralM[2*j + 1] = 0;
						
						population[i].locusS[2*j] = 1;
						population[i].locusS[2*j + 1] = 0;
						population[i].locusM[2*j] = 1;
						population[i].locusM[2*j + 1] = 0;
						population[i].morphe[j] = 0;
					}
					
					if( j%3 == 1 ){ // mid : all ssMm at generation 0
						population[i].neutralS[2*j] = 0;
						population[i].neutralS[2*j + 1] = 0;
						population[i].neutralM[2*j] = 1;
						population[i].neutralM[2*j + 1] = 0;
						
						population[i].locusS[2*j] = 0;
						population[i].locusS[2*j + 1] = 0;
						population[i].locusM[2*j] = 1;
						population[i].locusM[2*j + 1] = 0;
						population[i].morphe[j] = 1;
					}
					
					if( j%3 == 2 ){ // long : all ssmm at generation 0
						population[i].neutralS[2*j] = 0;
						population[i].neutralS[2*j + 1] = 0;
						population[i].neutralM[2*j] = 0;
						population[i].neutralM[2*j + 1] = 0;
						
						population[i].locusS[2*j] = 0;
						population[i].locusS[2*j + 1] = 0;
						population[i].locusM[2*j] = 0;
						population[i].locusM[2*j + 1] = 0;
						population[i].morphe[j] = 2;
					}
				}
				
				if( initialSituation == 1 ){ // all are ssmm (long)
					population[i].neutralS[2*j] = 0;
					population[i].neutralS[2*j + 1] = 0;
					population[i].neutralM[2*j] = 0;
					population[i].neutralM[2*j + 1] = 0;
					
					population[i].locusS[2*j] = 0;
					population[i].locusS[2*j + 1] = 0;
					population[i].locusM[2*j] = 0;
					population[i].locusM[2*j + 1] = 0;
					population[i].morphe[j] = 2;
				}
			
				if( initialSituation == 2 ){ // all are ssMM (mid)
					population[i].neutralS[2*j] = 0;
					population[i].neutralS[2*j + 1] = 0;
					population[i].neutralM[2*j] = 1;
					population[i].neutralM[2*j + 1] = 1;
					
					population[i].locusS[2*j] = 0;
					population[i].locusS[2*j + 1] = 0;
					population[i].locusM[2*j] = 1;
					population[i].locusM[2*j + 1] = 1;
					population[i].morphe[j] = 1;
				}
				
				if( initialSituation == 3 ){ // all are SSmm (short)
					population[i].neutralS[2*j] = 1;
					population[i].neutralS[2*j + 1] = 1;
					population[i].neutralM[2*j] = 0;
					population[i].neutralM[2*j + 1] = 0;
					
					population[i].locusS[2*j] = 1;
					population[i].locusS[2*j + 1] = 1;
					population[i].locusM[2*j] = 0;
					population[i].locusM[2*j + 1] = 0;
					population[i].morphe[j] = 0;
				}
				
				if( initialSituation == 4 ){ // all are SSMM (short)
					population[i].neutralS[2*j] = 1;
					population[i].neutralS[2*j + 1] = 1;
					population[i].neutralM[2*j] = 1;
					population[i].neutralM[2*j + 1] = 1;
					
					population[i].locusS[2*j] = 1;
					population[i].locusS[2*j + 1] = 1;
					population[i].locusM[2*j] = 1;
					population[i].locusM[2*j + 1] = 1;
					population[i].morphe[j] = 0;
				}
				
				//printf("%d\n%d\n", population[i].sexChro[2*j], population[i].sexChro[2*j+1]);
            			//population[i].nOffsprings[j] = floor(fecundity * population[i].femaleAllocation[j]) + gsl_ran_binomial(r, (fecundity * population[i].femaleAllocation[j]) - floor(fecundity * population[i].femaleAllocation[j]), 1);	// nOffs = floor(fecundity x femaleAllocation) + 1 according to a random Binomial integer
            			population[i].nOffsprings[j] = gsl_ran_poisson(r, fecundity); // for tristyli project: number of babies is the same for all

        }	// end of loop along the individuals
    }	// end of loop along the demes
}

void setToZero(const int nDemes, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]){
	// function to set all to demes to zero for the #of immigrants, the extinction status and the #of produced seeds.
	// before configMetapop()
	int i = 0;
	for(i=0; i<nDemes; i++){
		nImmigrants[i] = 0;
		extinctionStatus[i] = 0;
		nProducedSeeds[i] = 0;
	}
}

void setNumberIndividualsPerDeme(int nIndividuals[], const int nDemes, const int maxIndPerDem_1){
	int i = 0;
	for(i=0; i<nDemes; i++){
		nIndividuals[i] = maxIndPerDem_1;
	}
}

void setNumberIndividualsPerDeme_fluctuation(gsl_rng* r, const int nDemes, int nIndividuals[], const int maxIndPerDem_1, const int maxIndPerDem_2, const double P_transition_12_size, const double P_transition_21_size){
	int i = 0;
	unsigned int transition = 0;
	for(i=0; i<nDemes; i++){
		transition = 0;
		// if K1
		if(nIndividuals[i] == maxIndPerDem_1){ // if carrying capacity is K1
			transition = gsl_ran_binomial(r, P_transition_12_size, 1); // proba to switch from K1 to K2
			if( transition == 1 ){ // if transition
				nIndividuals[i] = maxIndPerDem_2; // switch from K1 to K2
			}
		}else{
			if(nIndividuals[i] == maxIndPerDem_2){
				transition = gsl_ran_binomial(r, P_transition_21_size, 1); // proba to switch from K2 to K1
				if( transition == 1 ){
					nIndividuals[i] = maxIndPerDem_1;
				}
			}
		}
	}
}

void setFMigColToZero(double f_mig_col[]){
	int i = 0;
	for(i=0; i<6; i++){
		f_mig_col[i] = 0.0;
	}
}

void configMetapop(gsl_rng* r, Deme* population, const int nDemes, const int nIndividuals[], const double migration, const double extinction, int nImmigrants[], int extinctionStatus[], int nProducedSeeds[]){
	// function to get per deme the #of immigrants received, the extinction status and the #of produced seeds.
	// after setToZero()
	// modifies nImmigrants[]; extinctionStatus[] and nProducedSeeds[]
	int i = 0;
	int j = 0;
	int nShort = 0;
	int nMid = 0;
	int nLong = 0;
	
	const unsigned int binomTrials = 1;

	for(i=0; i<nDemes; i++){	// start of the loop over the demes
		int nHeterogametic = 0;  // number of heterogametic individuals in the deme

		for(j=0; j<population[i].nIndividus; j++){
			if(population[i].sex[j] == 0){
				nHeterogametic += 1;
			}
		}
		
		if(nHeterogametic == population[i].nIndividus){
			printf("Deme with only one sex: %d\n", i);
			extinctionStatus[i] = 1;
		}
		
		nShort = 0;
		nMid = 0;
		nLong = 0;	
		for(j=0; j<population[i].nIndividus; j++){
			if(population[i].morphe[j] == 0){
				nShort = nShort + 1;
			}
			if(population[i].morphe[j] == 1){
				nMid = nMid + 1;
			}
			if(population[i].morphe[j] == 2){
				nLong = nLong + 1;
			}
			
		}
				
		nImmigrants[i] = gsl_ran_poisson(r, migration);
		if(nImmigrants[i] > nIndividuals[i]){
			nImmigrants[i] = nIndividuals[i];
		}
		extinctionStatus[i] = gsl_ran_binomial(r, extinction, binomTrials);	// 0: non-extincted; 1: extincted

		if(extinctionStatus[i] == 1){	// comment this block if you allow recolonization and migration occuring at the same time
			nImmigrants[i] = 0;	// no migrant if a deme is extinct
		}
		for(j=0; j<population[i].nIndividus; j++){	// start of the loop over the individuals
			nProducedSeeds[i] += population[i].nOffsprings[j];
		}	// end of the loop over the individuals
	}	// end of the loop over the demes
}

void initializeNewPopulation(Deme* newPopulation, const int nDemes, const int nIndividuals[], const int nProducedSeeds[], const int nNtrlLoci, const int nQuantiLoci, const double fecundity){
	// function to initialize the new population: allocate memory and set values to 0
	int i = 0;
	int j = 0;
	int k = 0;
	for(i=0; i<nDemes; i++){ // start of the loop along demes
		int taille = 0;
		taille = nProducedSeeds[i];	// size of the deme = nProducedSeeds within the parental population
		if(taille < fecundity){	// if not enough seeds are produced ==> deme is considered as extincted
			taille = fecundity;
		}
		if(taille > nIndividuals[i]){	// if too many individuals have to be present in the deme ==> cutoff to the nIndividuals[i] (=carrying capacity)
			taille = nIndividuals[i];
		}

		newPopulation[i].nIndividus = taille;
		newPopulation[i].ntrlLoci = NULL;
		newPopulation[i].quantiLoci = NULL;
		newPopulation[i].sexChro = NULL;
		newPopulation[i].sex = NULL;
		newPopulation[i].femaleAllocation = NULL;
		newPopulation[i].maleAllocation = NULL;
		newPopulation[i].nOffsprings = NULL;
		newPopulation[i].neutralS = NULL;
		newPopulation[i].neutralM = NULL;
		newPopulation[i].locusS = NULL;
		newPopulation[i].locusM = NULL;
		newPopulation[i].morphe = NULL;

		newPopulation[i].ntrlLoci = malloc(2 * taille * nNtrlLoci * sizeof(long));
		newPopulation[i].quantiLoci = malloc(2 * taille * nQuantiLoci * sizeof(long));
		newPopulation[i].sexChro = malloc(2 * taille * sizeof(int));
		newPopulation[i].sex = malloc(taille * sizeof(int));
		newPopulation[i].femaleAllocation = malloc(taille * sizeof(long));
		newPopulation[i].maleAllocation = malloc(taille * sizeof(long));
		newPopulation[i].nOffsprings = malloc(taille * sizeof(int));
		newPopulation[i].neutralS = malloc(2 * taille * sizeof(int));
		newPopulation[i].neutralM = malloc(2 * taille * sizeof(int));
		newPopulation[i].locusS = malloc(2 * taille * sizeof(int));
		newPopulation[i].locusM = malloc(2 * taille * sizeof(int));
		newPopulation[i].morphe = malloc(taille * sizeof(int));

		if(newPopulation[i].ntrlLoci == NULL || newPopulation[i].quantiLoci == NULL || newPopulation[i].sexChro == NULL || newPopulation[i].sex == NULL || newPopulation[i].femaleAllocation == NULL || newPopulation[i].maleAllocation == NULL || newPopulation[i].nOffsprings == NULL || newPopulation[i].neutralS == NULL || newPopulation[i].neutralM == NULL ||  newPopulation[i].locusS == NULL || newPopulation[i].locusM == NULL || newPopulation[i].morphe == NULL){
			exit(0);
		}

		for(j=0; j<taille; j++){ // start the loop along individuals
			newPopulation[i].femaleAllocation[j] = 0.0;
			newPopulation[i].maleAllocation[j] = 0.0;
			newPopulation[i].nOffsprings[j] = 0;

			newPopulation[i].sexChro[2*j] = 0;
			newPopulation[i].sexChro[2*j + 1] = 0;
			newPopulation[i].sex[j] = 0;

			newPopulation[i].neutralS[2*j] = 0;
			newPopulation[i].neutralS[2*j + 1] = 0;
			newPopulation[i].neutralM[2*j] = 0;
			newPopulation[i].neutralM[2*j + 1] = 0;

			newPopulation[i].locusS[2*j] = 0;
			newPopulation[i].locusS[2*j + 1] = 0;
			newPopulation[i].locusM[2*j] = 0;
			newPopulation[i].locusM[2*j + 1] = 0;
			newPopulation[i].morphe[j] = 0;
			
			for(k=0; k<(2*nNtrlLoci); k++){ // loop along the {2: diploid} x {nNtrlLoci: number of neutral loci} positions
				newPopulation[i].ntrlLoci[j*2*nNtrlLoci + k] = 0;
			}
			for(k=0; k<(2*nQuantiLoci); k++){       // loop along the {2: diploid} x {nQuantiLoci: number of quantitative loci} positions
				newPopulation[i].quantiLoci[j*2*nQuantiLoci + k] = 0;
			}
		}	// end of the loop along individuals
	} // end of the loop along demes

/*	for(i=0; i<nDemes; ++i){
		printf("%d\t", newPopulation[i].nIndividus);
	}
	printf("\n");*/
}

void panmixie(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const double ntrlMutation, const double quantiMutation, const double fecundity, const double reprodMale, const double sexAvantage, const int sexualSystem, const double selfingRate, const int currentGeneration, const int generationNewAllele, const int newAllele, const int generationNewAllele2, const int newAllele2, const int initialSituation, const double homomorphe_probability, const int freqDepSel, const int sysReprod){
	// function returning a new deme after a run of panmixia from the old deme
	// here: only "deme specific" events are simulated (meiosis + mutation)
	double currentAllelicEffect = 0.0;	// current allelic effect of a quantitative allele
	double minEffect = 0.0;	// minimum value that a quantitative mutation can take
	double maxEffect = 0.0;	// maximum value that a quantitative mutation can take
	int i = 0;
	int j = 0;
	int tmp = 0;
	int nMutation = 0;
	int autofec = 0;
	
	int nShort = 0;
	int nMid = 0;
	int nLong = 0;
	
	int n_non_short = 0;
	int n_non_mid = 0;
	int n_non_long = 0;
	
	int K = 0; // #of_individuals in the parental deme.
	int N = 0; // #of_babies to produce
	int N_reprodMales = 0; // number of fathers contributing to the next generation because of sexual selection among males
	
	int test_polymorphic = 0; // test whether there is only one morphe (=0), or multiple morphes in the deme (=1)
	
		for(i=0; i<nDemes; i++){	// start the loop over the nDemes

//			for(j=0; j<population[i].nIndividus; j++){
//				printf("%d %d %d %d %d\n", population[i].locusS[2*j], population[i].locusS[2*j+1], population[i].locusM[2*j], population[i].locusM[2*j+1], population[i].morphe[j]);
//			}

			K = 0;	// #of_individuals in the parental deme.
			N = 0;	// #of_babies to produce
			N_reprodMales = 0; // number of reproductie males 
			
			int* available_males = NULL; // males having won the sexual competition
			int* paternalIndexes = NULL; // competiting males
			int* parentalIndexes = NULL; // contain de the parental indexes (i.e:{0, 1, 2, 3, 4} if K == 5)
			int* mothers = NULL;	// array containing the mothers ID. Ex: {2, 19, 3, 3, 1} means that the first baby has individual#2 as mother, second baby has individual#19 as mother, babies 3 and 4 have individual#3
			int* fathers = NULL;	// array containing the fathers ID.
			double* fitnesses = NULL; // array containing the phenotypes of individuals
			
			K = population[i].nIndividus;	// #of_individuals in the parental deme
			N = newPopulation[i].nIndividus; 	// #of_produced_seeds by the K parents
			
			// Vector of all possible males, than can be putatively fathers
			paternalIndexes = malloc(K * sizeof(int));
			fitnesses = malloc(K * sizeof(double));
			if( paternalIndexes == NULL || fitnesses == NULL ){
				exit(0);
			}else{
				for(j=0; j<K; ++j){
					paternalIndexes[j] = j; // get the list of possible fathers
					fitnesses[j] = 0.0; // initialization
				}
			}
			
			// Get the number of reprodMales and draw a list of possible fathers (in available_males)
			if( reprodMale >= 1.0 ){
				N_reprodMales = K;
			}else{
				N_reprodMales = gsl_ran_poisson(r, reprodMale * K);
				if( N_reprodMales == 0 ){
					N_reprodMales = 1;
				}else{
					if( N_reprodMales > K ){
						N_reprodMales = K;
					}
				}
			}
			
			// chose the parents
			parentalIndexes = malloc(K * sizeof(int));
			mothers = malloc(N * sizeof(int));	// indexes of mothers of the N autochtones within the deme of size newPopulation[i].nIndividus
			fathers = malloc(N * sizeof(int));
			
			int* non_short = NULL; // indexes of non-short styled individuals
			int* non_mid = NULL;
			int* non_long = NULL;
			int* short_individuals = NULL; // indexes of short styled individuals
			int* mid_individuals = NULL;
			int* long_individuals = NULL;
			
			// block to bring a first new allele if needed	
			if( currentGeneration == generationNewAllele && initialSituation >0 ){ // start the block responsible of the new mutation at S and M loci
				if( i == 0 ){ // in the first deme
					j = 0;
					j = gsl_rng_uniform_int( r, population[i].nIndividus);
					
					if( newAllele == 0 ){
						population[i].neutralS[2*j] = 1;
						population[i].locusS[2*j] = 1;
					}
					
					if( newAllele == 1 ){
						population[i].neutralS[2*j] = 0;
						population[i].locusS[2*j] = 0;
					}
					
					if( newAllele == 2 ){
						population[i].neutralM[2*j] = 1;
						population[i].locusM[2*j] = 1;
					}
					
					if( newAllele == 3 ){
						population[i].neutralM[2*j] = 0;
						population[i].locusM[2*j] = 0;
					}
					
					if( population[i].locusS[2*j] + population[i].locusS[2*j+1] != 0){
						population[i].morphe[j] = 0;
					}else{
						if( population[i].locusM[2*j] + population[i].locusM[2*j+1] == 0){
							population[i].morphe[j] = 2;
						}else{
							population[i].morphe[j] = 1;
						}
					}
				}
			} // end of the block responsible of the new mutation at S and M loci
			
			// block to bring a second new allele if needed	
			if( currentGeneration == generationNewAllele2 && initialSituation >0 ){ // start the block responsible of the new mutation at S and M loci
				if( i == 0 ){ // in the first deme
					j = 0;
					j = gsl_rng_uniform_int( r, population[i].nIndividus);
					
					if( newAllele2 == 0 ){
						population[i].neutralS[2*j] = 1;
						population[i].locusS[2*j] = 1;
					}
					
					if( newAllele2 == 1 ){
						population[i].neutralS[2*j] = 0;
						population[i].locusS[2*j] = 0;
					}
					
					if( newAllele2 == 2 ){
						population[i].neutralM[2*j] = 1;
						population[i].locusM[2*j] = 1;
					}
					
					if( newAllele2 == 3 ){
						population[i].neutralM[2*j] = 0;
						population[i].locusM[2*j] = 0;
					}
					
					if( population[i].locusS[2*j] + population[i].locusS[2*j+1] != 0){
						population[i].morphe[j] = 0;
					}else{
						if( population[i].locusM[2*j] + population[i].locusM[2*j+1] == 0){
							population[i].morphe[j] = 2;
						}else{
							population[i].morphe[j] = 1;
						}
					}
				}
			} // end of the block responsible of the new mutation at S and M loci
			
			// loop to count the number of individuals for each morphe in the parental pop
			// nShort, nMid, nLong count the number of S, M and L morphes in the deme : among all individuals
			nShort = 0;
			nMid = 0;
			nLong = 0;
			for(j=0; j<K; ++j){
				// counting the number of different morphes
				if(population[i].morphe[j] == 0){
					nShort = nShort + 1;
				}
				if(population[i].morphe[j] == 1){
					nMid = nMid + 1;
				}
				if(population[i].morphe[j] == 2){
					nLong = nLong + 1;
				}
			}
			
			// test the floral polymorphism
			test_polymorphic = 0;
			if( nShort == K || nMid == K || nLong == K ){
				test_polymorphic = 0;
			}else{
				test_polymorphic = 1;
			}
			
			// loop to get the fitnesses depending on the frequency of individual phenotypes
			for( j=0; j<K; ++j){
				if( test_polymorphic == 0 || freqDepSel == 0){ // if monomorphic or no frequency-dependent subsampling, then all individuals have the same fitness
					fitnesses[j] = 1;
				}else{
					if( population[i].morphe[j] == 0 ){
						fitnesses[j] = 1 - nShort / (1.0 * K);
					}
					if( population[i].morphe[j] == 1 ){
						fitnesses[j] = 1 - nMid / (1.0 * K);
					}
					if( population[i].morphe[j] == 2 ){
						fitnesses[j] = 1 - nLong / (1.0 * K);
					}
				}
			}

			// get the reproductive males
			available_males = malloc(N_reprodMales * sizeof(int));
			if( available_males == NULL ){
				exit(0);
			}else{
				//gsl_ran_choose(r, available_males, N_reprodMales, paternalIndexes, K, sizeof(int)); // all males have the same proba of being father
				weightedSample_without_replacement(r, paternalIndexes, fitnesses, available_males, K, N_reprodMales); //
			}
			
/*			printf("nShort=%d\tnMid=%d\tnLong=%d\n", nShort, nMid, nLong);
			for( j=0; j<N_reprodMales; ++j){
				printf("%d\t", available_males[j]);
			}
			printf("\n");
			for( j=0; j<N_reprodMales; ++j){
				printf("%d\t", population[i].morphe[available_males[j]]);
			}
			printf("\n");
			for( j=0; j<N_reprodMales; ++j){
				printf("%lf\t", fitnesses[available_males[j]]);
			}
			printf("\n");
*/
			// nShort, nMid, nLong count the number of S, M and L morphes in the deme : among all available males
			nShort = 0;
			nMid = 0;
			nLong = 0;
			
			for( j=0; j<N_reprodMales; ++j){
				if( population[i].morphe[ available_males[j] ] == 0 ){
					nShort = nShort + 1;
				}
				
				if( population[i].morphe[ available_males[j] ] == 1 ){
					nMid = nMid + 1;
				}
				
				if( population[i].morphe[ available_males[j] ] == 2 ){
					nLong = nLong + 1;
				}
			}
			
			// test the floral polymorphism among reproductive males
			test_polymorphic = 0;
			if( nShort == N_reprodMales || nMid == N_reprodMales || nLong == N_reprodMales ){
				test_polymorphic = 0;
			}else{
				test_polymorphic = 1;
			}
			
			n_non_short = N_reprodMales - nShort;
			n_non_mid = N_reprodMales - nMid;
			n_non_long = N_reprodMales - nLong;
			non_short = malloc( n_non_short * sizeof(int) );
			non_mid = malloc( n_non_mid * sizeof(int) );
			non_long = malloc( n_non_long * sizeof(int) );
			short_individuals = malloc( nShort * sizeof(int) );
			mid_individuals = malloc( nMid * sizeof(int) );
			long_individuals = malloc( nLong * sizeof(int) );
		
				
			if(parentalIndexes == NULL || mothers == NULL || fathers == NULL || non_short == NULL || non_mid == NULL || non_long == NULL || short_individuals == NULL || mid_individuals == NULL || long_individuals == NULL){
				exit(0);
			}
			// get the short individuals among reproductive males
			tmp = 0;
			for(j=0; j<N_reprodMales; j++){
				if(population[i].morphe[ available_males[j] ] == 0){
					short_individuals[tmp] = available_males[j];
					tmp++;
				}
			}
			
			// get the mid individuals among reproductive males
			tmp = 0;
			for(j=0; j<N_reprodMales; j++){
				if(population[i].morphe[ available_males[j] ] == 1){
					mid_individuals[tmp] = available_males[j];
					tmp++;
				}
			}
			
			// get the long individuals among reproductive males
			tmp = 0;
			for(j=0; j<N_reprodMales; j++){
				if(population[i].morphe[ available_males[j] ] == 2){
					long_individuals[tmp] = available_males[j];
					tmp++;
				}
			}
			
			
			// get the non_short individuals among reproductive males
			tmp = 0;
			for(j=0; j<N_reprodMales; j++){
				if(population[i].morphe[ available_males[j] ] != 0){
					non_short[tmp] = available_males[j];
					tmp++;
				}
			}
			
			// get the non_mid individuals among reproductive males
			tmp = 0;
			for(j=0; j<N_reprodMales; j++){
				if(population[i].morphe[ available_males[j] ] != 1){
					non_mid[tmp] = available_males[j];
					tmp++;
				}
			}
			
			// get the non_long individuals among reproductive males
			tmp = 0;
			for(j=0; j<N_reprodMales; j++){
				if(population[i].morphe[ available_males[j] ] != 2){
					non_long[tmp] = available_males[j];
					tmp++;
				}
			}
			
			
			for(j=0; j<K; j++){
				parentalIndexes[j] = j;
			}
			
			// Sampling the parents
			for(j=0; j<N; j++){
				mothers[j] = gsl_rng_uniform_int(r, K); // get the mothers uniformally over the deme
				fathers[j] = 0;
			}

			
			// get the fathers
			// const int sysReprod = atoi(argv[30]); // if 0: panmixia; if 1: trisStyle
			if( sysReprod == 1 ){ // if triStyle system (i.e: mating between different morphes, versus, panmixia
				if( test_polymorphic == 1 ){
					unsigned int test_homomorphic_cross = 0; // test if a given cross occur among individuals of same morphes
					int n_i = 0;
					int n_tot = 0;
					float q1 = 0.0;
					float q2 = 0.0;
					float r1 = 0.0;
					
					// crosses between different morphes if the deme is polymorphic
					for(j=0; j<N; j++){
						if( population[i].morphe[mothers[j]] == 0 ){
							//
							n_i = N_reprodMales - n_non_short;
							n_tot = N_reprodMales;
						
							q1 = ( n_i * homomorphe_probability * 1.0) / n_tot;
							q2 = ( n_tot - n_i  * 1.0) / n_tot;
							r1 = (n_i * homomorphe_probability * 1.0) / n_tot / (q1 + q2);
							//
							
							test_homomorphic_cross = 0;
	//						test_homomorphic_cross = gsl_ran_binomial(r, (N_reprodMales-n_non_short)*(1-homomorphe_probability)/((N_reprodMales-n_non_short) + n_non_short), 1);
							test_homomorphic_cross = gsl_ran_binomial(r, r1, 1);
							if ( test_homomorphic_cross == 1){
								//fathers[j] = gsl_rng_uniform_int(r, N_reprodMales);
								fathers[j] = short_individuals[ gsl_rng_uniform_int(r, nShort) ];
							}else{
								fathers[j] = non_short[ gsl_rng_uniform_int(r, n_non_short ) ];
							}
						}
						
						if( population[i].morphe[mothers[j]] == 1 ){
							//
							n_i = N_reprodMales - n_non_mid;
							n_tot = N_reprodMales;
						
							q1 = ( n_i * homomorphe_probability * 1.0) / n_tot;
							q2 = ( n_tot - n_i  * 1.0) / n_tot;
							r1 = (n_i * homomorphe_probability * 1.0) / n_tot / (q1 + q2);
							//
							
							test_homomorphic_cross = 0;
	//						test_homomorphic_cross = gsl_ran_binomial(r, (N_reprodMales-n_non_short)*(1-homomorphe_probability)/((N_reprodMales-n_non_short) + n_non_short), 1);
							test_homomorphic_cross = gsl_ran_binomial(r, r1, 1);
							if ( test_homomorphic_cross == 1){
								//fathers[j] = gsl_rng_uniform_int(r, N_reprodMales);
								fathers[j] = mid_individuals[ gsl_rng_uniform_int(r, nMid) ];
							}else{	
								fathers[j] = non_mid[ gsl_rng_uniform_int(r, n_non_mid ) ];
							}
						}
						
						if( population[i].morphe[mothers[j]] == 2 ){
							//
							n_i = N_reprodMales - n_non_long;
							n_tot = N_reprodMales;
						
							q1 = ( n_i * homomorphe_probability * 1.0) / n_tot;
							q2 = ( n_tot - n_i  * 1.0) / n_tot;
							r1 = (n_i * homomorphe_probability * 1.0) / n_tot / (q1 + q2);
							//
							
							test_homomorphic_cross = 0;
	//						test_homomorphic_cross = gsl_ran_binomial(r, (N_reprodMales-n_non_short)*(1-homomorphe_probability)/((N_reprodMales-n_non_short) + n_non_short), 1);
							test_homomorphic_cross = gsl_ran_binomial(r, r1, 1);
							if ( test_homomorphic_cross == 1){
								//fathers[j] = gsl_rng_uniform_int(r, N_reprodMales);
								fathers[j] = long_individuals[ gsl_rng_uniform_int(r, nLong) ];
							}else{
						
								fathers[j] = non_long[ gsl_rng_uniform_int(r, n_non_long ) ];
							}
						}
						//printf("polymorphic mother %d father %d\n", population[i].morphe[mothers[j]], population[i].morphe[fathers[j]]);
					}
				}else{
				// all individuals share the same probability of 1/N of being the father if the deme is monomorphic
					for(j=0; j<N; j++){
						fathers[j] = available_males[ gsl_rng_uniform_int(r, N_reprodMales) ];
						//printf("monomorphic mother %d father %d\n", population[i].morphe[mothers[j]], population[i].morphe[fathers[j]]);
					}
				}
			}else{ // if panmixia
				for(j=0; j<N; j++){
					fathers[j] = gsl_rng_uniform_int(r, K); // get the fathers uniformally over the deme
				}
			}
			
			if(selfingRate > 0){
				// loop over hermaphroditic parents for dealing with self-fertilization at rates "selfing"
				for(j=0; j<N; j++){ // loop over the N mother
					autofec = gsl_ran_binomial(r, selfingRate, 1);
					if(population[i].sex[j] == 1){ // if the mother is a hermaphrodite
						if(autofec == 1){ // if self-fertilization
							fathers[j] = mothers[j];
						}
					}
				}
			}
			// print parents
/*			for( j=0; j<N; ++j){
				printf("%dx%d\t", mothers[j], fathers[j]);
			}
			printf("\n");
			for( j=0; j<N; ++j){
				printf("%d%d%d%d\t", population[i].locusS[mothers[j]*2], population[i].locusS[mothers[j]*2+1],population[i].locusM[mothers[j]*2], population[i].locusM[mothers[j]*2+1]);
			}
			printf("\n");
			for( j=0; j<N; ++j){
				printf("%d%d%d%d\t", population[i].locusS[fathers[j]*2], population[i].locusS[fathers[j]*2+1],population[i].locusM[fathers[j]*2], population[i].locusM[fathers[j]*2+1]);
			}
			printf("\n");
*/

			// Meiosis and transmission of gametes
			int pos = 0;	// position in deme.ntrlLoci (or deme.quantiLoci) of offsprings
			int pos2 = 0;	// position in deme.ntrlLoci (or deme.quantiLoci) of parents
			for(j=0; j<N; j++){	// loop over the babies.  Transmission of alleles, meiosis = gsl_ran_binomial
				for(tmp=0; tmp<nNtrlLoci; tmp++){	// loop along the neutral loci
					pos = j*nNtrlLoci*2 + tmp*2 + 0;	// position where to put the allele from the mother
					pos2 = mothers[j]*nNtrlLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // binomial transmission of the parental alleles with p=0.5 (=recombination at meiosis)
					newPopulation[i].ntrlLoci[pos] = population[i].ntrlLoci[pos2];

					pos += 1;	// position where to put the allele from the father = position of mother_allele + 1
					pos2 = fathers[j]*nNtrlLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[i].ntrlLoci[pos] = population[i].ntrlLoci[pos2];
				}	// end of loop along the neutral loci

				for(tmp=0; tmp<nQuantiLoci; tmp++){	// loop along the quantitative loci
					pos = j*nQuantiLoci*2 + tmp*2 + 0;
					pos2 = mothers[j]*nQuantiLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[i].quantiLoci[pos] = population[i].quantiLoci[pos2];

					pos += 1;
					pos2 = fathers[j]*nQuantiLoci*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1);
					newPopulation[i].quantiLoci[pos] = population[i].quantiLoci[pos2];
				}	// end of loop along the quantitative loci

				// transmission of the sex chromosomes
				for(tmp=0; tmp<1; tmp++){
					pos = j*2 + tmp*2 + 0;
	               			pos2 = mothers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // sex chromosome from the mother
			                newPopulation[i].sexChro[pos] = population[i].sexChro[pos2];

					pos += 1;
                			pos2 = fathers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // sex chromosome from the father
		        	        newPopulation[i].sexChro[pos] = population[i].sexChro[pos2];
				}
				
				// transmission of neutralS	
				for(tmp=0; tmp<1; tmp++){
					pos = j*2 + tmp*2;
					pos2 = mothers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // neutralS from the mother
					newPopulation[i].neutralS[pos] = population[i].neutralS[pos2];
					
					pos += 1;
					pos2 = fathers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // neutralS from the father
					newPopulation[i].neutralS[pos] = population[i].neutralS[pos2];	
				}
				
				// transmission of neutralM
				for(tmp=0; tmp<1; tmp++){
					pos = j*2 + tmp*2;
					pos2 = mothers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // neutralM from the mother
					newPopulation[i].neutralM[pos] = population[i].neutralM[pos2];
					
					pos += 1;
					pos2 = fathers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // neutralM from the father
					newPopulation[i].neutralM[pos] = population[i].neutralM[pos2];	
				}
			
				// transmission of locusS	
				for(tmp=0; tmp<1; tmp++){
					pos = j*2 + tmp*2;
					pos2 = mothers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // locusS from the mother
					newPopulation[i].locusS[pos] = population[i].locusS[pos2];
					
					pos += 1;
					pos2 = fathers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // locusS from the father
					newPopulation[i].locusS[pos] = population[i].locusS[pos2];	
				}
				
				// transmission of locusM
				for(tmp=0; tmp<1; tmp++){
					pos = j*2 + tmp*2;
					pos2 = mothers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // locusM from the mother
					newPopulation[i].locusM[pos] = population[i].locusM[pos2];
					
					pos += 1;
					pos2 = fathers[j]*2 + tmp*2 + gsl_ran_binomial(r, 0.5, 1); // locusM from the father
					newPopulation[i].locusM[pos] = population[i].locusM[pos2];	
				}
			}	// end of loop along the individuals

			for(j=0; j<N; j++){	// loop over babies. Determine their sex 0 (heterogametic) or 1 (homogametic)
				if(newPopulation[i].sexChro[2*j] == newPopulation[i].sexChro[2*j + 1]){ // if homogametic, sex = 1
					newPopulation[i].sex[j] = 1;
				}else{
					newPopulation[i].sex[j] = 0; // if heterogametic, sex = 0
				}
			}
			
			for(j=0; j<N; j++){	// loop over babies. Determine their morph
				if( (newPopulation[i].locusS[2*j] + newPopulation[i].locusS[2*j + 1 ]) == 0){ // if individual j has no big S allele
					if( (newPopulation[i].locusM[2*j] + newPopulation[i].locusM[2*j + 1 ]) == 0){ // if individual j has no big M allele
						newPopulation[i].morphe[j] = 2; // long-styled morphe
					}else{
						newPopulation[i].morphe[j] = 1; // mid-styled morphe
					}
				}else{ // if individual j carriers at least one big S allele
					newPopulation[i].morphe[j] = 0; // then short-styled morphe
				}
			}

			// Put mutations
			nMutation = 0;
			nMutation = gsl_ran_binomial(r, ntrlMutation, 2*N*nNtrlLoci);
			if(nMutation > 0 ){
				for(j=0; j<nMutation; j++){
					pos = rand()%(2*N*nNtrlLoci);
					newPopulation[i].ntrlLoci[pos] =  gsl_rng_uniform_int(r, MAX_NUMBER_OF_INITIAL_NTRL_ALLELES-1) + 1; // ntrl alleles in [1, MAX_NUMBER_OF_INITIAL_NTRL_ALLELES[
				}
			}
			
			nMutation = 0;
			nMutation = gsl_ran_binomial(r, quantiMutation, 2*N*nQuantiLoci);
			if(nMutation > 0){
				for(j=0; j<nMutation; j++){
					pos = rand()%(2*N*nQuantiLoci);
					currentAllelicEffect = newPopulation[i].quantiLoci[pos];
					minEffect = (1 - RANGE) * currentAllelicEffect;
					if(minEffect<0){
						minEffect = 0.0;
					}
					maxEffect = (1 + RANGE) * currentAllelicEffect;
					if(maxEffect>1.0/2/nQuantiLoci){
					maxEffect = 1.0/2/nQuantiLoci;
					}
					newPopulation[i].quantiLoci[pos] = gsl_ran_flat(r, minEffect, maxEffect);
				}
			}

			pos = 0;
			for(j=0; j<N; j++){	// loop along the individuals to calculate the new femaleAllocation and the number of offsprings
				tmp = 0;
				do{
					newPopulation[i].femaleAllocation[j] += newPopulation[i].quantiLoci[pos];	// sum of allelic effects within individuals
					pos +=1;
					tmp +=1;
				}while(tmp<2*nQuantiLoci);

				newPopulation[i].maleAllocation[j] = 1 - newPopulation[i].femaleAllocation[j];	// set the male allocation as (1 - femaleAllocation)
//				newPopulation[i].nOffsprings[j] = floor(fecundity * newPopulation[i].femaleAllocation[j]) + gsl_ran_binomial(r, (fecundity * newPopulation[i].femaleAllocation[j]) - floor(fecundity * newPopulation[i].femaleAllocation[j]), 1);	// set the #of_offsprings produced for each individual
				newPopulation[i].nOffsprings[j] = gsl_ran_poisson(r, fecundity);



				if(newPopulation[i].sex[j] == 0){ // if heterogametic individual
					if(sexualSystem == 1){  // if XY system
						newPopulation[i].femaleAllocation[j] = 0;
						newPopulation[i].maleAllocation[j] = sexAvantage;
						newPopulation[i].nOffsprings[j] = 0;
					}
					if(sexualSystem == 2){  // if ZW system
						newPopulation[i].femaleAllocation[j] = sexAvantage;
						newPopulation[i].maleAllocation[j] = 0;
						double nBabies = fecundity * sexAvantage * 1.0;
						newPopulation[i].nOffsprings[j] = floor(nBabies) + gsl_ran_binomial(r, nBabies - floor(nBabies), 1);
					}
				}
			}
			free(paternalIndexes);
			free(parentalIndexes);
			free(fitnesses);
			free(available_males);
			free(mothers);
			free(fathers);
			free(non_short);
			free(non_mid);
			free(non_long);
			free(short_individuals);
			free(mid_individuals);
			free(long_individuals);
	}	// end of loop over the nDemes
}

void replacement(gsl_rng* r, Deme* population, Deme* newPopulation, const int nDemes, const int nIndividuals[], const int nNtrlLoci, const int nQuantiLoci, const int nImmigrants[], int nProducedSeeds[], int extinctionStatus[], const int recolonization, const int generation, const int sexualSystem, const int colonizationModel, double* f_mig_col){
	int i = 0;
	int j = 0;
	int k = 0;
	int compteur = 0;
	int taille = 0;
	int demeTMP = 0;
	int indTMP = 0;
	int posDonneur = 0;	// first position in the emigrant population
	int posReceveur = 0;	// first position in the imigrant population
	int nExtinctedDemes = 0;	// number of extincted demes.
	double* indexOfDemes = NULL;
	double* nProducedSeedsDouble = NULL;

	// count the number of extincted demes
	for(i=0; i<nDemes; i++){
		if(extinctionStatus[i] == 1){
			nExtinctedDemes += 1;
		}
	}

	if(nExtinctedDemes == nDemes){
		printf("Demes are dead at generation: %d\n", generation);
		exit(1);
	}

	indexOfDemes = malloc((nDemes - nExtinctedDemes) * sizeof(double)); // there are (nDemes - nExtinctedDemes) allow to receive or send migrants/recolonizer
	nProducedSeedsDouble = malloc((nDemes - nExtinctedDemes) * sizeof(double));

	if(indexOfDemes == NULL || nProducedSeedsDouble == NULL){
		exit(0);
	}

	for(i=0; i<nDemes; i++){
		if(extinctionStatus[i] != 1){
			indexOfDemes[compteur] = i;
			nProducedSeedsDouble[compteur] = (double)nProducedSeeds[i];
			compteur += 1;
		}
	}

	for(i=0; i<nDemes; i++){
		taille = newPopulation[i].nIndividus + nImmigrants[i]; // total number of individuals in a deme: autochtones + migrants
		if(extinctionStatus[i] == 1){
			taille = recolonization;
		}
		if(taille > nIndividuals[i]){
			taille = nIndividuals[i];
		}

		population[i].nIndividus = taille;
		population[i].ntrlLoci = malloc(2 * taille * nNtrlLoci * sizeof(long));
		population[i].quantiLoci = malloc(2 * taille * nQuantiLoci * sizeof(long));
		population[i].sexChro = malloc(2 * taille * sizeof(int));
		population[i].sex = malloc(taille * sizeof(int));
		population[i].femaleAllocation = malloc(taille * sizeof(long));
		population[i].maleAllocation = malloc(taille * sizeof(long));
		population[i].nOffsprings = malloc(taille * sizeof(int));
		population[i].neutralS = malloc(2 * taille * sizeof(int));
		population[i].neutralM = malloc(2 * taille * sizeof(int));
		population[i].locusS = malloc(2 * taille * sizeof(int));
		population[i].locusM = malloc(2 * taille * sizeof(int));
		population[i].morphe = malloc(taille * sizeof(int));


		if(population[i].ntrlLoci == NULL || population[i].quantiLoci == NULL || population[i].sexChro == NULL || population[i].sex == NULL || population[i].femaleAllocation == NULL || population[i].maleAllocation == NULL || population[i].nOffsprings == NULL || population[i].neutralS == NULL || population[i].neutralM == NULL || population[i].locusS == NULL || population[i].locusM == NULL || population[i].morphe == NULL){
			exit(0);
		}

		// copy paste newPopulation into population of non extincted demes
		if(extinctionStatus[i] == 0){
			for(j=0; j<(2 * (population[i].nIndividus - nImmigrants[i]) * nNtrlLoci); j++){
				population[i].ntrlLoci[j] = newPopulation[i].ntrlLoci[j];
			}

			for(j=0; j<(2 * (population[i].nIndividus - nImmigrants[i]) * nQuantiLoci); j++){
				population[i].quantiLoci[j] = newPopulation[i].quantiLoci[j];
			}

			for(j=0; j<(2 * (population[i].nIndividus - nImmigrants[i])); j++){
                		population[i].sexChro[j] = newPopulation[i].sexChro[j];
                		population[i].neutralS[j] = newPopulation[i].neutralS[j];
                		population[i].neutralM[j] = newPopulation[i].neutralM[j];
                		population[i].locusS[j] = newPopulation[i].locusS[j];
                		population[i].locusM[j] = newPopulation[i].locusM[j];
			}

			for(j=0; j<(population[i].nIndividus - nImmigrants[i]); j++){
				population[i].femaleAllocation[j] = newPopulation[i].femaleAllocation[j];
				population[i].maleAllocation[j] = newPopulation[i].maleAllocation[j];
				population[i].nOffsprings[j] = newPopulation[i].nOffsprings[j];
				population[i].sex[j] = newPopulation[i].sex[j];
				population[i].morphe[j] = newPopulation[i].morphe[j];
			}
			// add migrants
			if(nImmigrants[i] > 0){
				double* emigrantDemes = NULL;
				emigrantDemes = malloc(nImmigrants[i] * sizeof(double)); // indexes of demes where immigrants come from
				if(emigrantDemes == NULL){
					exit(0);
				}

				weightedSample(r, indexOfDemes, nProducedSeedsDouble, emigrantDemes, nDemes-nExtinctedDemes, nImmigrants[i]); // return the indexes of nImmigrants[i] emigrant demes from nDemes

				for(j=0; j<nImmigrants[i]; j++){	// loop over the nImmigrants for the deme "i"
					demeTMP = emigrantDemes[j];	// get the emigrant deme
					indTMP = gsl_ran_flat(r, 0, newPopulation[demeTMP].nIndividus);	// randomly choose the emigrant individual from the emigrant deme
					
					// neutral loci
					posDonneur = 2 * nNtrlLoci * indTMP + 0;
					posReceveur = (2* nNtrlLoci * population[i].nIndividus) - (2 * nNtrlLoci * nImmigrants[i]) + 2 * nNtrlLoci * j;

					for(k=0; k < 2 * nNtrlLoci; k++){	// loop over positions to bring migrants through copy-pasting
						population[i].ntrlLoci[posReceveur + k] = newPopulation[demeTMP].ntrlLoci[posDonneur + k];
					}

					// quantitative loci
					posDonneur = 2 * nQuantiLoci * indTMP + 0;
					posReceveur = (2*nQuantiLoci * population[i].nIndividus) - (2 * nQuantiLoci * nImmigrants[i]) + 2 * nQuantiLoci * j ;

					for(k=0; k < 2 * nQuantiLoci; k++){
						population[i].quantiLoci[posReceveur + k] = newPopulation[demeTMP].quantiLoci[posDonneur + k];
					}

					// sex chromosomes
					posDonneur = 2 * indTMP + 0;
					posReceveur = (2 * population[i].nIndividus) - (2 * nImmigrants[i]) + 2 * j;

					for(k=0; k<2; k++){
						population[i].sexChro[posReceveur + k] = newPopulation[demeTMP].sexChro[posDonneur +k];
					}
					
					// neutralS
					posDonneur = 2 * indTMP + 0;
					posReceveur = (2 * population[i].nIndividus) - (2 * nImmigrants[i]) + 2 * j;

					for(k=0; k<2; k++){
						population[i].neutralS[posReceveur + k] = newPopulation[demeTMP].neutralS[posDonneur +k];
					}
					
					// neutralM
					posDonneur = 2 * indTMP + 0;
					posReceveur = (2 * population[i].nIndividus) - (2 * nImmigrants[i]) + 2 * j;

					for(k=0; k<2; k++){
						population[i].neutralM[posReceveur + k] = newPopulation[demeTMP].neutralM[posDonneur +k];
					}
					
					
					// locusS
					posDonneur = 2 * indTMP + 0;
					posReceveur = (2 * population[i].nIndividus) - (2 * nImmigrants[i]) + 2 * j;

					for(k=0; k<2; k++){
						population[i].locusS[posReceveur + k] = newPopulation[demeTMP].locusS[posDonneur +k];
					}
					
					// locusM
					posDonneur = 2 * indTMP + 0;
					posReceveur = (2 * population[i].nIndividus) - (2 * nImmigrants[i]) + 2 * j;

					for(k=0; k<2; k++){
						population[i].locusM[posReceveur + k] = newPopulation[demeTMP].locusM[posDonneur +k];
					}
					
					
					// copy paste femaleAllocation, maleAllocation and nOffsprings, morphe
					posDonneur = indTMP;
					posReceveur = population[i].nIndividus - nImmigrants[i] + j;

					population[i].femaleAllocation[posReceveur] = newPopulation[demeTMP].femaleAllocation[posDonneur];
					population[i].maleAllocation[posReceveur] = newPopulation[demeTMP].maleAllocation[posDonneur];
					population[i].nOffsprings[posReceveur] = newPopulation[demeTMP].nOffsprings[posDonneur];

					population[i].sex[posReceveur] = newPopulation[demeTMP].sex[posDonneur];
					
					population[i].morphe[posReceveur] = newPopulation[demeTMP].morphe[posDonneur];
					
					// compute frequencies of migrant morphes
					if( population[i].morphe[posReceveur] == 0 ){
						f_mig_col[0]++;
					}else{
						if( population[i].morphe[posReceveur] == 1){
							f_mig_col[1]++;
						}else{
							if( population[i].morphe[posReceveur] == 2){
								f_mig_col[2]++;
							}
						}
					}
					
				}	// end of loop over the n immigrants to the deme "i"

				free(emigrantDemes);
			}	// end of migrant traitement
		}	// end of treatment of non-extincted demes
		if(extinctionStatus[i] == 1){
			double* emigrantDemes = NULL;
			emigrantDemes = malloc(population[i].nIndividus * sizeof(double));
			if(emigrantDemes == NULL){
				exit(0);
			}

			weightedSample(r, indexOfDemes, nProducedSeedsDouble, emigrantDemes, (nDemes - nExtinctedDemes), population[i].nIndividus); // return the indexes of nImmigrants[i] emigrant demes from nDemes. When extinctionStatus==1, all individuals making the deme are "imigrants". Here, emigrant demes are different --> migrant pool model
			if(colonizationModel == 1 && recolonization > 1){ // --> if the user specified a propagule pool model
				for(j = 1; j<population[i].nIndividus; j++){
					emigrantDemes[j] = emigrantDemes[0]; // all emigrant demes are the same for a propagule pool model
				}
			}

			for(j=0; j<population[i].nIndividus; j++){	// loop over the individuals to put into extincted demes
				demeTMP = emigrantDemes[j];

				// if we only have hermaphrodites (sexualSystem == 0): all individuals have the same probability of being colonizers 
				if(sexualSystem == 0){
					indTMP = gsl_ran_flat(r, 0, newPopulation[demeTMP].nIndividus);	// randomly choose the colonizer individual from the emigrant deme
					// DEBUG debug printf("%d %d %d\n", demeTMP, indTMP, newPopulation[demeTMP].morphe[indTMP]);
				}else{
					// choose the colonizer individual by avoiding unisexuals. Only cosexuals can recolonize 
					double* indexOfIndividuals = NULL; // indexOfIndividuals = [0, 1, 2, ..., N-1] if there are N individuals in the deme i
					double* recolonizer = NULL; // vector of size 1 containing the sampled individual contributing to recolonization
					double* sexDouble = NULL; // convert the vector of (int) sexes [0, 0, 1, 0, 1, ... ] in doubles [0.0, 0.0, 1.0, 1.0, ... ]
					indexOfIndividuals = malloc(newPopulation[demeTMP].nIndividus * sizeof(double)); 
					recolonizer = malloc(1 * sizeof(double));
					sexDouble = malloc(newPopulation[demeTMP].nIndividus * sizeof(double));
					if(indexOfIndividuals == NULL || recolonizer == NULL){
						exit(0);
					}
					for(k=0; k<newPopulation[demeTMP].nIndividus; k++){
						indexOfIndividuals[k] = (double)k;
						sexDouble[k] = (double)newPopulation[demeTMP].sex[k];
					}
				
					weightedSample(r, indexOfIndividuals, sexDouble, recolonizer, newPopulation[demeTMP].nIndividus, 1);
					indTMP = recolonizer[0];
					free(indexOfIndividuals);
					free(recolonizer);
					free(sexDouble);
				}

				// neutral loci
				posDonneur = 2 * nNtrlLoci * indTMP + 0;
				for(k=0; k< 2 * nNtrlLoci; k++){
					population[i].ntrlLoci[2 * nNtrlLoci * j + k] = newPopulation[demeTMP].ntrlLoci[2 * nNtrlLoci * indTMP + k];
				}

				// quantitative loci
				posDonneur = 2 * nQuantiLoci * indTMP + 0;
				for(k=0; k< 2 * nQuantiLoci; k++){
					population[i].quantiLoci[2 * nQuantiLoci * j + k] = newPopulation[demeTMP].quantiLoci[2 * nQuantiLoci * indTMP + k];
				}

		                // sex chromosomes
				posDonneur = 2 * indTMP + 0;
				for(k=0; k< 2; k++){
					population[i].sexChro[2 * j + k] = newPopulation[demeTMP].sexChro[2 * indTMP + k];
				}

				// neutralS
				posDonneur = 2 * indTMP + 0;
				for(k=0; k< 2; k++){
					population[i].neutralS[2 * j + k] = newPopulation[demeTMP].neutralS[2 * indTMP + k];
				}
					
				// neutralM
				posDonneur = 2 * indTMP + 0;
				for(k=0; k< 2; k++){
					population[i].neutralM[2 * j + k] = newPopulation[demeTMP].neutralM[2 * indTMP + k];
				}
				
				// locusS
				posDonneur = 2 * indTMP + 0;
				for(k=0; k< 2; k++){
					population[i].locusS[2 * j + k] = newPopulation[demeTMP].locusS[2 * indTMP + k];
				}
					
				// locusM
				posDonneur = 2 * indTMP + 0;
				for(k=0; k< 2; k++){
					population[i].locusM[2 * j + k] = newPopulation[demeTMP].locusM[2 * indTMP + k];
				}
				
				// copy paste femaleAllocation, maleAllocation and nOffsprings
				population[i].femaleAllocation[j] = newPopulation[demeTMP].femaleAllocation[indTMP];
				population[i].maleAllocation[j] = newPopulation[demeTMP].maleAllocation[indTMP];
				population[i].nOffsprings[j] = newPopulation[demeTMP].nOffsprings[indTMP];

				population[i].sex[j] = newPopulation[demeTMP].sex[indTMP];
				
				population[i].morphe[j] = newPopulation[demeTMP].morphe[indTMP];

				// compute frequencies of colonizing morphes
				if( population[i].morphe[j] == 0 ){
					f_mig_col[3]++;
				}else{
					if( population[i].morphe[j] == 1){
						f_mig_col[4]++;
					}else{
						if( population[i].morphe[j] == 2){
							f_mig_col[5]++;
						}
					}
				}
			}	// end of loop over the individuals to put into extincted demes
			free(emigrantDemes);
		}	// end of treatment of extincted demes
	}	// end of loop over demes

	free(indexOfDemes);
	free(nProducedSeedsDouble);
}

void weightedSample(gsl_rng* r, const double* liste, const double* weights, double* target, const int sizeOfListe, const int nTrials){
	// function that fills the vector 'target' of size 'nTrials' containing the weighted-sampled 'sizeOfListe' elements of the 'liste':
	// weightedSample(gsl_rng* r, {2, 4, 6, 8, 10}, {1.2, 0.6, 0.3, 0.15, 0.05}, target, 5, 20)
	// target = {6, 6, 6, 2, 2, 2, 2, 2, 8, 2, 8, 6, 4, 2, 2, 4, 4, 4, 4, 2}
	// but can also be used for boolean sampling (pile ou face) using:
	// weightedSample(gsl_rng* r, {0, 1}, {1, 1}, target, 2, 1)
	int i = 0;
	unsigned int* n = NULL;
	int* sampledListe = NULL;
	n = malloc(sizeOfListe * sizeof(double));	// will contain the number of succes after K nTrials for each of the sizeOfListe elements of liste
	sampledListe = malloc(nTrials * sizeof(int));	// if n={0, 3, 1, 1}, sampledListe={1, 1, 1, 4, 5}

	gsl_ran_multinomial(r, sizeOfListe, nTrials, weights, n);	// return in 'n' the number of success for the sizeOfListe elements of liste

	int nValues = 0;
	int tmp = 0;
	int tmp2 = 0;
	for(i=0; i<sizeOfListe; i++){ // loop along the list called 'n' resulting from gsl_ran_multinomial
		nValues = n[i];
		if(nValues != 0){
			tmp2 = 0;
			do{
				sampledListe[tmp] = i;
				tmp++;
				tmp2++;
			}while(tmp2 < nValues);
		}
	}

	// shuffle values of the sampledListe
	gsl_permutation* p = gsl_permutation_alloc (nTrials);
	gsl_permutation_init (p);
	gsl_ran_shuffle(r, p -> data, nTrials, sizeof(size_t));

	tmp = 0;
	for(i=0; i<nTrials; i++){
		tmp=gsl_permutation_get(p, i);
		target[i]=liste[sampledListe[tmp]];
	}
	gsl_permutation_free(p);
	free(n);
	free(sampledListe);
}

void weightedSample_without_replacement(gsl_rng* r, const int* liste, const double* weights, int* target, const int sizeOfListe, const int nTrials){
	// function that fills the vector 'target' of size 'nTrials' containing the weighted-sampled 'sizeOfListe' elements of the 'liste':
	// weightedSample(gsl_rng* r, {2, 4, 6, 8, 10}, {1.2, 0.6, 0.3, 0.15, 0.05}, target, 5, 20)
	// target = {6, 6, 6, 2, 2, 2, 2, 2, 8, 2, 8, 6, 4, 2, 2, 4, 4, 4, 4, 2}
	// but can also be used for boolean sampling (pile ou face) using:
	// weightedSample(gsl_rng* r, {0, 1}, {1, 1}, target, 2, 1)
	int i = 0;
	double* liste_2 = NULL;
	double* weights_tmp = NULL;
	double* target_tmp = NULL;
	liste_2 = malloc( sizeOfListe * sizeof(double) );
	weights_tmp = malloc( sizeOfListe * sizeof(double) );
	target_tmp = malloc( 1 * sizeof(double) );
	
	if( liste_2 == NULL || weights_tmp == NULL || target_tmp == NULL ){
		exit(0);
	}else{
		for(i=0; i<sizeOfListe; ++i){
			liste_2[i] = liste[i];
			weights_tmp[i] = weights[i];
		}
		
		for(i=0; i<nTrials; ++i){
			weightedSample(r, liste_2, weights_tmp, target_tmp, sizeOfListe, 1);
			target[i] = target_tmp[0];
			weights_tmp[ (int)target_tmp[0] ] = 0.0; // remove the sampled individual from the list by atributing a proba of 0 to be sampled at the next rounds

		}
		free(weights_tmp);
		free(target_tmp);
		free(liste_2);
	}
}

void libererMemoirePopulation(Deme* population, const int nDemes){
	// free memory taken by the population at the end of each generation.
	int i = 0;
	for(i=0; i<nDemes; i++){
		free(population[i].ntrlLoci);
		free(population[i].quantiLoci);
		free(population[i].sexChro);
		free(population[i].sex);
		free(population[i].femaleAllocation);
		free(population[i].maleAllocation);
		free(population[i].nOffsprings);
		free(population[i].neutralS);
		free(population[i].neutralM);
		free(population[i].locusS);
		free(population[i].locusM);
		free(population[i].morphe);
	}
}

void writeNindividuals(const Deme* population, const int nDemes, const double extinction, const double migration, const int seed){
	int i = 0;

	char nameOfFileNInd[100];
	char nameOfFileFemAlloc[100];
	sprintf(nameOfFileNInd, "indOverTime_e%f_i%f_seed%d.txt", extinction, migration, seed);
	sprintf(nameOfFileFemAlloc, "femAllocOverTime_e%f_i%f_seed%d.txt", extinction, migration, seed);

	FILE* fichierNInd = NULL;
	FILE* fichierFemAlloc = NULL;
	fichierNInd =  fopen(nameOfFileNInd, "a");
	fichierFemAlloc =  fopen(nameOfFileFemAlloc, "a");
	if(fichierNInd != NULL && fichierFemAlloc != NULL){
		for(i=0; i<nDemes; i++){
			fprintf(fichierNInd, "%d ", population[i].nIndividus);
			fprintf(fichierFemAlloc, "%f ", gsl_stats_mean(population[i].femaleAllocation, 1, population[i].nIndividus));
		}
		fprintf(fichierNInd, "\n");
		fclose(fichierNInd);
		fprintf(fichierFemAlloc, "\n");
		fclose(fichierFemAlloc);
	}
}

void afficherPopulation(Deme* population, const int nDemes, const int nNtrlLoci, const int nQuantiLoci, const int sexualSystem){
	// called to print some informations about population in a debug mode
	int i = 0;
	int j = 0;
	int k = 0;
	int sexGenotype = 0;
	char sex = 'H';
	for(i=0; i<nDemes; i++){
		for(j=0; j<population[i].nIndividus; j++){
			if(population[i].sexChro[2*j] == population[i].sexChro[2*j + 1]){
				sexGenotype = 0;	// homogametic
			}
			if(population[i].sexChro[2*j] != population[i].sexChro[2*j + 1]){
				sexGenotype = 1;	// heterogametic
			}
			if(sexualSystem == 0){
				sex = 'H';
			}
			if(sexualSystem == 1){
				if(sexGenotype == 0){
					sex = 'H';
				}
				if(sexGenotype == 1){
					sex = 'M';
				}
			}
			if(sexualSystem == 2){
				if(sexGenotype == 0){
					sex = 'H';
				}
				if(sexGenotype == 1){
					sex = 'F';
				}
			}
			
			printf("Deme: %d Ind: %d Ntrl: ", i, j);
			for(k=0; k<(2*nNtrlLoci); k++){
				printf("%ld ", population[i].ntrlLoci[2*j*nNtrlLoci+k]);
			}
			printf(" Quanti: ");
			for(k=0; k<(2*nQuantiLoci); k++){
				printf("%.4lf ", population[i].quantiLoci[2*j*nQuantiLoci+k]);
			}
		printf("femAlloc: %.4lf malAlloc: %4.lf nOffs: %d sex: %d sex2: %c\n", population[i].femaleAllocation[j], population[i].maleAllocation[j], population[i].nOffsprings[j], population[i].sex[j], sex);
		}
	}
}

void genePop(Deme* population, const int nDemes, const int nNtrlLoci, const int seed, int time){
	//	generates an input for genepop (Rousset), launch it and clean the output to save space on the cluster
	int i = 0;
	int j = 0;
	int k = 0;
	int allele = 0;
	int cntLoci = 0;
	char nameOfGenePopFile[100];
	char nameOfSettingFile[100];
	char nameOfROutputFile[100];
//	char commandLineOne[100]; // launch genepop
//	char commandLineTwo[200]; // treat the genepop's output
	char commandLineThree[200]; // clean the tmp files
	char commandLineDiveRsity[200]; // R command calling diveRsity 
	sprintf(nameOfGenePopFile, "genepop_%d_%d.txt", time, seed);
	sprintf(nameOfSettingFile, "setting_%d_%d.txt", time, seed);
	sprintf(nameOfROutputFile, "output_diveRsity_%d_%d", time, seed);

//	sprintf(commandLineOne, "Genepop settingsFile=%s Mode=Batch >/dev/null", nameOfSettingFile); // external call of genepop (Rousset)
//	sprintf(commandLineTwo, "tail -n%d %s.FST | grep 'Locus'>tmp_%d.txt; tail -n%d %s.FST | grep 'All' >>tmp_%d.txt; mv tmp_%d.txt %s.FST", nNtrlLoci + 5, nameOfGenePopFile, seed, nNtrlLoci + 5, nameOfGenePopFile, seed, seed, nameOfGenePopFile); // external formating of genepop's output.
	sprintf(commandLineThree, "rm -rf cmdline.txt fichier.in %s %s", nameOfGenePopFile, nameOfSettingFile);

	sprintf(commandLineDiveRsity, "diveRsity.R input=%s output=%s", nameOfGenePopFile, nameOfROutputFile);

	FILE* fichierGenePop = NULL;
	FILE* settingFile = NULL;

	fichierGenePop = fopen(nameOfGenePopFile, "a");	// creating genepop file with data
	if(fichierGenePop != NULL){	// writing in file

		fprintf(fichierGenePop, "Simulated data\n");
	
		for(i=0; i<nNtrlLoci; i++){
			fprintf(fichierGenePop, "Locus%d\n", i);
		}

		for(i=0; i<nDemes; i++){	// loop over demes: start
			fprintf(fichierGenePop, "Pop\n");
			for(j=0; j<population[i].nIndividus; j++){	// loop over individuals of deme 'i': start
				fprintf(fichierGenePop, "Ind%d, ", j);
				for(k=0; k<(2*nNtrlLoci); k++){
					cntLoci += 1;
					allele = population[i].ntrlLoci[2*j*nNtrlLoci+k];
					if(allele < 10){
						fprintf(fichierGenePop, "00%d", allele);
					}
					if(allele >= 10 && allele < 100){
						fprintf(fichierGenePop, "0%d", allele);
					}
					if(allele >= 100){
						fprintf(fichierGenePop, "%d", allele);
					}
					if(cntLoci == 2){
						cntLoci = 0;
						fprintf(fichierGenePop, " ");
					}
				}
				fprintf(fichierGenePop, "\n");
			}	// loop over individuals of deme 'i': end
		}	// loop over demes: end
	}	// end of writing in file
	fclose(fichierGenePop);	// input file for genepop is generated

	settingFile = fopen(nameOfSettingFile, "a");	// creating setting file for genepop
	if(settingFile != NULL){
		fprintf(settingFile, "GenepopInputFile=%s\nMenuOptions=6.1\n", nameOfGenePopFile);
	}
	fclose(settingFile);

//	int testCommandLineOne = system(commandLineOne); // call genepop
//	if(testCommandLineOne == -1){
//		exit(1); // check the returned value
//	}

	int testNameOfROutputFile = system(commandLineDiveRsity); // call diveRsity R
	if(testNameOfROutputFile == -1){
		exit(1); // check the returned value
	}

//	int testCommandLineTwo = system(commandLineTwo); // reformat the genepop's output
//	if(testCommandLineTwo == -1){
//		exit(1); // check the returned value
//	}
	
	int testCommandLineThree = system(commandLineThree); // remove genePop formated file and setting file
	if(testCommandLineThree == -1){
		exit(1);
	}

}

void statisticsPopulations(Deme* population, const int nDemes, const int maxIndPerDem_1, const int maxIndPerDem_2, const double P_transition_12_size, const double P_transition_21_size, const int nQuantiLoci, const double fecundity, const double migration, const double extinction, const int recolonization, const int sexualSystem, const double sexAvantage, const int seed, int time, const double selfingRate, const int colonizationModel, const double global_fst_cm, const double global_fst_coal, const double global_gpst, const double global_D, const double global_Fis, const double homomorphe_probability, const int sysReprod){
	// function that calculates the mean female allocation, its standard deviation and the percentage of cosexuals in the metapopulation
	char breedingSystem[50];
	// if 0: panmixia; if 1: trisStyle
	if( sysReprod == 0 ){
		sprintf(breedingSystem, "panmixia");
	}
	if( sysReprod == 1 ){
		sprintf(breedingSystem, "heteromorphy");
	}


	// also computes FST_var = var(p)/(1-p) and FST_coal = (Htot - Hs) / Htot
	int i = 0;
	int j = 0;
	int cnt = 0;
	int cnt2 = 0;
	double fstValue = 0.0;
	double fstRoussetValue = 0.0;
	double fstValueDensity = 0.0;
	double meanAllocFemale = 0.0;
	double sdAllocFemale = 0.0;
	double meanAllocFemaleCosexual = 0.0;
	double sdAllocFemaleCosexual = 0.0;
	double cosexualProportion = 0.0;
	int nIndividusTotal = 0;

	// frequencies of the different morphes (short, mid, long) and alleles (S, s, M, m) averaged over the WHOLE metapopulation
	double f_short = 0.0;
	double f_mid = 0.0;
	double f_long = 0.0;
	
	double f_S_neutral = 0.0;
	double f_M_neutral = 0.0;
	
	double f_S = 0.0;
	double f_s = 0.0;
	double f_M = 0.0;
	double f_m = 0.0;

	double f_SSMM = 0.0;
	double f_SSMm = 0.0;
	double f_SSmm = 0.0;
	double f_SsMM = 0.0;
	double f_SsMm = 0.0;
	double f_Ssmm = 0.0;
	double f_ssMM = 0.0;
	double f_ssMm = 0.0;
	double f_ssmm = 0.0;

	// vectors of double of size 'nDemes' containing frequencies of different morphes for each deme
	double* freq_short_demes = NULL;
	double* freq_mid_demes = NULL;
	double* freq_long_demes = NULL;
	freq_short_demes = malloc(nDemes * sizeof(double));
	freq_mid_demes = malloc(nDemes * sizeof(double));
	freq_long_demes = malloc(nDemes * sizeof(double));
	
	if( freq_short_demes == NULL || freq_mid_demes == NULL || freq_long_demes == NULL ){
		exit(0);
	}
	
	for(i=0; i<nDemes; i++){
		freq_short_demes[i] = 0.0;
		freq_mid_demes[i] = 0.0;
		freq_long_demes[i] = 0.0;
	}

	// file name
	char nomFichierSortie[200];
	sprintf(nomFichierSortie, "output_%d.txt", seed);
	FILE* fichierSortie = NULL;

	fichierSortie = fopen(nomFichierSortie, "r");
	if(fichierSortie == NULL){
		fichierSortie = fopen(nomFichierSortie, "a");
		fprintf(fichierSortie, "nDemes\tnIndMaxPerDeme_1\tIndMaxPerDeme_2\tP_fluctuation_K1_K2\tP_fluctuation_K2_k1\tNtot\tnQuantiLoci\tselfingRate\tproba_homomorphic_pairing\tfecundity\tmigRate\textRate\tcolonizationModel\trecolonization\tmatingSystem\tatGeneration\tsexSystem\tsexAvantage\tseed\tavg_f_short_metapop\tsd_f_short_demes\tavg_f_mid_metapop\tsd_f_mid_demes\tavg_f_long_metapop\tsd_f_long_demes\tf_S_neutral\tf_M_neutral\tf_S\tf_s\tf_M\tf_m\tf_SSMM\tf_SSMm\tf_SSmm\tf_SsMM\tf_SsMm\tf_Ssmm\tf_ssMM\tf_ssMm\tf_ssmm\tnDemes_shortMidLong_polymorphism\tnDemes_short_lostOnly\tnDemes_mid_lostOnly\tnDemes_long_lostOnly\tnDemes_short_fixed\tnDemes_mid_fixed\tnDemes_long_fixed\tnDemes_allele_Sneutral_lost\tnDemes_allele_Sneutral_polym\tnDemes_allele_Sneutral_fixed\tnDemes_allele_Mneutral_lost\tnDemes_allele_Mneutral_polym\tnDemes_allele_Mneutral_fixed\tnDemes_allele_S_lost\tnDemes_allele_S_polym\tnDemes_allele_S_fixed\tnDemes_allele_M_lost\tnDemes_allele_M_polym\tnDemes_allele_M_fixed\tmeanFemAlloc\tsdFemAlloc\tmeanFemAllocCosexual\tsdFemAllocCosexual\tcosexualProportion\tobsFST_var\tobsFST_coal\tobsGST_p\tobsJostD\tobsFIS\texpFST_Nmax\texpFST_Nobs\texpFST_Rousset_Nmax\tf_mig_short\tf_mig_mid\tf_mig_long\tf_col_short\tf_col_mid\tf_col_long\n");
		fclose(fichierSortie);
	}else{
		fclose(fichierSortie);
	}
	
	fichierSortie = fopen(nomFichierSortie, "a");

	if(fichierSortie != NULL){	
	
		for(i=0; i<nDemes; i++){
			nIndividusTotal += population[i].nIndividus;
		}	
	
		double* allocFemale = NULL; // female allocation in the whole metapopulation
		double* allocFemaleCosexual = NULL; // female allocation of cosexuals only. allocFemale = allocFemaleCosexual if sexualSystem == 0
		
		int nShort_deme_i = 0;
		int nMid_deme_i = 0;
		int nLong_deme_i = 0;
		int nDemes_shortMidLong_polymorphism = 0;
		int nDemes_short_lostOnly = 0;
		int nDemes_mid_lostOnly = 0;
		int nDemes_long_lostOnly = 0;
		int nDemes_short_fixed = 0;
		int nDemes_mid_fixed = 0;
		int nDemes_long_fixed = 0;
	
		int f_S_neutral_i = 0;
		int f_S_i = 0;
		int f_M_neutral_i = 0;
		int f_M_i = 0;
	
		int nDemes_Sneutral_lost = 0;
		int nDemes_Sneutral_polym = 0;
		int nDemes_Sneutral_fix = 0;
		
		int nDemes_Mneutral_lost = 0;
		int nDemes_Mneutral_polym = 0;
		int nDemes_Mneutral_fix = 0;
	
		int nDemes_S_lost = 0;
		int nDemes_S_polym = 0;
		int nDemes_S_fix = 0;
		
		int nDemes_M_lost = 0;
		int nDemes_M_polym = 0;
		int nDemes_M_fix = 0;
		
		for(i=0; i<nDemes; i++){
			f_S_neutral_i = 0;
			f_S_i = 0;
			f_M_neutral_i = 0;
			f_M_i = 0;
	
			nShort_deme_i = 0;
			nMid_deme_i = 0;
			nLong_deme_i = 0;
			for(j=0; j<population[i].nIndividus; j++){
				cosexualProportion += population[i].sex[j]; // sex[j] = 0 if unisexual; sex[j] + 1 if cosexual
				// morphe
				if(population[i].morphe[j] == 0){
					nShort_deme_i += 1;
					freq_short_demes[i] += 1.0/population[i].nIndividus;
					f_short += 1;
				}else{ 
					if(population[i].morphe[j] == 1){
						nMid_deme_i += 1;
						freq_mid_demes[i] += 1.0/population[i].nIndividus;
						f_mid += 1;
					}else{ 
						if(population[i].morphe[j] == 2){
							nLong_deme_i += 1;
							freq_long_demes[i] += 1.0/population[i].nIndividus;
							f_long += 1;
						}
					}
				} 
				
			
				// neutralS
				if(population[i].neutralS[2*j] == 1){ f_S_neutral += 1; f_S_neutral_i += 1;}
				if(population[i].neutralS[2*j + 1] == 1){ f_S_neutral += 1; f_S_neutral_i += 1;}

				// neutralM
				if(population[i].neutralM[2*j] == 1){ f_M_neutral += 1; f_M_neutral_i += 1;}
				if(population[i].neutralM[2*j + 1] == 1){ f_M_neutral += 1; f_M_neutral_i += 1;}
				
				// locusS
				if(population[i].locusS[2*j] == 0){ f_s += 1;}else{ f_S += 1; f_S_i += 1;}
				if(population[i].locusS[2*j + 1] == 0){ f_s += 1;}else{ f_S += 1; f_S_i += 1;}

				// locusM
				if(population[i].locusM[2*j] == 0){ f_m += 1;}else{ f_M += 1; f_M_i += 1;}
				if(population[i].locusM[2*j + 1] == 0){ f_m += 1;}else{ f_M += 1; f_M_i += 1;}
				
				// genotypes
				// SS
				if(population[i].locusS[2*j] + population[i].locusS[2*j+1] == 2){
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 2){f_SSMM += 1;} // MM
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 1){f_SSMm += 1;} // Mm
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 0){f_SSmm += 1;} // mm
					}
				// Ss
				if(population[i].locusS[2*j] + population[i].locusS[2*j+1] == 1){
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 2){f_SsMM += 1;} // MM
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 1){f_SsMm += 1;} // Mm
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 0){f_Ssmm += 1;} // mm
					}
				// ss
				if(population[i].locusS[2*j] + population[i].locusS[2*j+1] == 0){
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 2){f_ssMM += 1;} // MM
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 1){f_ssMm += 1;} // Mm
					if(population[i].locusM[2*j] + population[i].locusM[2*j+1] == 0){f_ssmm += 1;} // mm
					}

				}
			// nDemes morphe
			if( nShort_deme_i == population[i].nIndividus ){
				// f_short == 1
				nDemes_short_fixed++; // 1/0/0
			}else{
				if( nShort_deme_i == 0 ){
					// f_short == 0
					if( nMid_deme_i == population[i].nIndividus ){
						// f_mid == 1
						nDemes_mid_fixed++; // 0/1/0
					}else{
						if( nMid_deme_i == 0 ){
							// f_mid == 0
							nDemes_long_fixed++; // 0/0/1
						}else{
							// f_mid in ]0, 1[
							nDemes_short_lostOnly++; // 0/x/x
						}
					}
				}else{
					// f_short in ]0, 1[
					if( nMid_deme_i == 0 ){
						// f_mid == 0
						nDemes_mid_lostOnly++; // x/0/x
					}else{
						// f_mid in ]0, 1[
						if( nLong_deme_i == 0 ){
							// f_long == 0
							nDemes_long_lostOnly++; // x/x/0
						}else{
							nDemes_shortMidLong_polymorphism++; // x/x/x
						}
					}
				}
			}

			// nDemes Sneutral
			if( f_S_neutral_i == 2*population[i].nIndividus ){
				nDemes_Sneutral_fix += 1;
			}else{
				if( f_S_neutral_i == 0 ){
					nDemes_Sneutral_lost += 1;
				}else{nDemes_Sneutral_polym += 1;}
			}
			
			// nDemes S
			if( f_S_i == 2*population[i].nIndividus ){
				nDemes_S_fix += 1;
			}else{
				if( f_S_i == 0 ){
					nDemes_S_lost += 1;
				}else{nDemes_S_polym += 1;}
			}
			
			// nDemes Mneutral
			if( f_M_neutral_i == 2*population[i].nIndividus ){
				nDemes_Mneutral_fix += 1;
			}else{
				if( f_M_neutral_i == 0 ){
					nDemes_Mneutral_lost += 1;
				}else{nDemes_Mneutral_polym += 1;}
			}
			
			// nDemes M
			if( f_M_i == 2*population[i].nIndividus ){
				nDemes_M_fix += 1;
			}else{
				if( f_M_i == 0 ){
					nDemes_M_lost += 1;
				}else{nDemes_M_polym += 1;}
			}
		} // end of loop over demes
		
		f_short = f_short / (1.0 * nIndividusTotal);
		f_mid = f_mid / (1.0 * nIndividusTotal);
		f_long = f_long / (1.0 * nIndividusTotal);
		f_S_neutral = f_S_neutral / (2.0 * nIndividusTotal);
		f_M_neutral = f_M_neutral / (2.0 * nIndividusTotal);
		f_S = f_S / (2.0 * nIndividusTotal);
		f_s = f_s / (2.0 * nIndividusTotal);
		f_M = f_M / (2.0 * nIndividusTotal);
		f_m = f_m / (2.0 * nIndividusTotal);
		
		f_SSMM /= (1.0 * nIndividusTotal);
		f_SSMm /= (1.0 * nIndividusTotal);
		f_SSmm /= (1.0 * nIndividusTotal);
		f_SsMM /= (1.0 * nIndividusTotal);
		f_SsMm /= (1.0 * nIndividusTotal);
		f_Ssmm /= (1.0 * nIndividusTotal);
		f_ssMM /= (1.0 * nIndividusTotal);
		f_ssMm /= (1.0 * nIndividusTotal);
		f_ssmm /= (1.0 * nIndividusTotal);
		
		allocFemale = malloc(nIndividusTotal * sizeof(double));
		allocFemaleCosexual = malloc(cosexualProportion * sizeof(double));

		for(i=0; i<nDemes; i++){
			for(j=0; j<population[i].nIndividus; j++){
				allocFemale[cnt] = population[i].femaleAllocation[j];
				if(population[i].sex[j] == 1){
					allocFemaleCosexual[cnt2] = population[i].femaleAllocation[j];
					cnt2 += 1;
				}
				cnt += 1;
			}
		}

		meanAllocFemale = gsl_stats_mean(allocFemale, 1, nIndividusTotal);
		sdAllocFemale = gsl_stats_sd(allocFemale, 1, nIndividusTotal);
		meanAllocFemaleCosexual = gsl_stats_mean(allocFemaleCosexual, 1, (int) cosexualProportion);
		sdAllocFemaleCosexual = gsl_stats_sd(allocFemaleCosexual, 1, (int) cosexualProportion);
		
		cosexualProportion = cosexualProportion / nIndividusTotal;
	
		fstValue = fstMullon(maxIndPerDem_1, extinction, recolonization, migration); // expected fst assuming that all demes are full
		fstValueDensity = fstMullon((int) nIndividusTotal/(1.0*nDemes), extinction, recolonization, migration); // expected fst, using the the average number of individuals in demes (for cases with high extinction, low fecundity)
		fstRoussetValue = fstRousset(maxIndPerDem_1, extinction, recolonization, migration, colonizationModel);

		free(allocFemale);
		free(allocFemaleCosexual);
		
		char colonizationModelTMP[50];
		
		if(colonizationModel == 0){
			sprintf(colonizationModelTMP, "migrationPool");
//			colonizationModelTMP = "migrationPool";
		}
		if(colonizationModel == 1){
			sprintf(colonizationModelTMP, "propagulePool");
//			colonizationModelTMP = "propagulePool";
		}

/*nDemes\tnIndMaxPerDeme_1\tIndMaxPerDeme_2\tP_fluctuation_K1_K2\tP_fluctuation_K2_k1\tNtot\tnQuantiLoci\tselfingRate\tproba_homomorphic_pairing\tfecundity\tmigRate\textRate\tcolonizationModel\trecolonization\tatGeneration\tsexSystem\tsexAvantage\tseed\tavg_f_short_metapop\tsd_f_short_demes\tavg_f_mid_metapop\tsd_f_mid_demes\tavg_f_long_metapop\tsd_f_long_demes\tf_S_neutral\tf_M_neutral\tf_S\tf_s\tf_M\tf_m\tf_SSMM\tf_SSMm\tf_SSmm\tf_SsMM\tf_SsMm\tf_Ssmm\tf_ssMM\tf_ssMm\tf_ssmm\tnDemes_Sfixed\tnDemes_Mfixed\tnDemes_Lfixed\tnDemes_allele_Sneutral_lost\tnDemes_allele_Sneutral_polym\tnDemes_allele_Sneutral_fixed\tnDemes_allele_Mneutral_lost\tnDemes_allele_Mneutral_polym\tnDemes_allele_Mneutral_fixed\tnDemes_allele_S_lost\tnDemes_allele_S_polym\tnDemes_allele_S_fixed\tnDemes_allele_M_lost\tnDemes_allele_M_polym\tnDemes_allele_M_fixed\tmeanFemAlloc\tsdFemAlloc\tmeanFemAllocCosexual\tsdFemAllocCosexual\tcosexualProportion\tobsFST_var\tobsFST_coal\tobsGST_p\tobsJostD\tobsFIS\texpFST_Nmax\texpFST_Nobs\texpFST_Rousset_Nmax\tf_mig_short\tf_mig_mid\tf_mig_long\tf_col_short\tf_col_mid\tf_col_long\n");*/
		fprintf(fichierSortie,"%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%s\t%d\t%s\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t", nDemes, maxIndPerDem_1, maxIndPerDem_2, P_transition_12_size, P_transition_21_size, nIndividusTotal, nQuantiLoci, selfingRate, homomorphe_probability, fecundity, migration, extinction, colonizationModelTMP, recolonization, breedingSystem, time, sexualSystem, sexAvantage, seed, f_short, gsl_stats_sd(freq_short_demes, 1, nDemes), f_mid, gsl_stats_sd(freq_mid_demes, 1, nDemes), f_long, gsl_stats_sd(freq_long_demes, 1, nDemes), f_S_neutral, f_M_neutral, f_S, f_s, f_M, f_m, f_SSMM, f_SSMm, f_SSmm, f_SsMM, f_SsMm, f_Ssmm, f_ssMM, f_ssMm, f_ssmm, nDemes_shortMidLong_polymorphism, nDemes_short_lostOnly, nDemes_mid_lostOnly, nDemes_long_lostOnly, nDemes_short_fixed, nDemes_mid_fixed, nDemes_long_fixed, nDemes_Sneutral_lost, nDemes_Sneutral_polym, nDemes_Sneutral_fix, nDemes_Mneutral_lost, nDemes_Mneutral_polym, nDemes_Mneutral_fix, nDemes_S_lost, nDemes_S_polym, nDemes_S_fix, nDemes_M_lost, nDemes_M_polym, nDemes_M_fix, meanAllocFemale, sdAllocFemale, meanAllocFemaleCosexual, sdAllocFemaleCosexual, cosexualProportion, global_fst_cm, global_fst_coal, global_gpst, global_D, global_Fis, fstValue, fstValueDensity, fstRoussetValue);
//%d\t%d\t%d\t%lf\t%lf\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%s\t%d\t%d\t%d\t%lf\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf
		fclose(fichierSortie);
	}
	
	free(freq_short_demes);
	free(freq_mid_demes);
	free(freq_long_demes);
}

void statisticsMigrantsColonizers(const int seed, const double* f_mig_col){
	// printf("%lf %lf %lf %lf %lf %lf\n", f_mig_col[0], f_mig_col[1], f_mig_col[2], f_mig_col[3], f_mig_col[4], f_mig_col[5]); // print frequencies of short/mid/long among migrants/colonizers
	char nomFichierSortie[200];
	sprintf(nomFichierSortie, "output_%d.txt", seed);
	FILE* fichierSortie = NULL;

	fichierSortie = fopen(nomFichierSortie, "a");
	if(fichierSortie != NULL){	
		double n_mig = 0.0;
		double n_col = 0.0;
		
		n_mig = f_mig_col[0] + f_mig_col[1] + f_mig_col[2];
		n_col = f_mig_col[3] + f_mig_col[4] + f_mig_col[5];

		fprintf(fichierSortie, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", f_mig_col[0]/n_mig, f_mig_col[1]/n_mig, f_mig_col[2]/n_mig, f_mig_col[3]/n_col, f_mig_col[4]/n_col, f_mig_col[5]/n_col);
		fclose(fichierSortie);
	}
}


void checkCommandLine(int argc){
	if(argc != 31){
		printf("\n%sThe number of provided arguments is not correct.%s\nYou provided %s%d%s argument while %s29%s are expected:\n\t\
		%s1.%s Number of demes (>0)\n\n\t\
		%s2.%s Carrying capacity 1 (int) : all demes share this K1 at generation 0\n\t\
		%s3.%s Carrying capacity 2 (int)\n\t\
		%s4.%s Proba of transition from K1 to K2 (float)\n\t\
		%s5.%s Proba of transition from K2 to K1 (float)\n\n\t\
		%s6.%s Number of generations (>0)\n\t\
		%s7.%s Number of generations when extinction starts to occur\n\n\t\
		%s8.%s Number of neutral loci (>=0)\n\t\
		%s9.%s Neutral mutation rate (in [0-1]))\n\t\
		%s10.%s Number of quantitative loci (>0; USELESS FOR THE MOMENT, SET IT TO 1)\n\t\
		%s11.%s Quantitative mutation rate (in [0-1]; USELESS FOR THE MOMENT, SET IT TO 0)\n\n\t\
		%s12.%s Max number of offsprings per hermaphrodite (>0)\n\t\
		%s13.%s Proportion of individuals reproducing through the male function (in [0, 1])\n\n\t\
		%s14.%s Immigration rate (Poisson distributed; >=0)\n\t\
		%s15.%s Extinction rate (Binomialy distributed; in [0-1])\n\t\
		%s16.%s Number of individuals recolonizing an extincted deme (>0)\n\n\t\
		%s17.%s Colonization model, 'migration pool' (=0) or 'propagule pool' (=1) models\n\n\t\
		%s18.%s sexualSystem is equal to 0 if autosomal, equal to 1 if XY and equal to 2 if ZW (USELESS FOR THE MOMENT, SET IT TO 0)\n\t\
		%s19.%s Sexual effects of heterogametic sex (if equal to 1.5 in XY system, males have a 50 percent advantage to sire available ovules). Required but neglected if sexualSystem == 0 (USELESS FOR THE MOMENT, SET IT TO 1)\n\n\t\
		%s20.%s Selfing rate of hermaphrodites, fixed over time (in [0-1])\n\n\t\
		%s21.%s frequency of statistics calculation, all X generations (positive integer)\n\n\t\
		%s22.%s Seed for the random generator (>0)\n\n\t\
		%s23.%s initial situation is the phenotypic state of the metapopulation at generation 0 (0: freq is 1/3 of each morphes. 1: only ss mm. 2: only ss MM. 3: only SS mm. 4: only SS MM)\n\n\t\
		%s24.%s FIRST newly introduced allele in deme 0 at generation generationNewAllele in a single individual (0 = S; 1 = s; 2 = M; 3 = m), this argument is not considered if argument 19 is equal to 0.\n\n\t\
		%s25.%s generation at which ONE copy of the FIRST new allele is brought into the metapop. This argument is required but not considered if argument 20 is equal to 0.\n\n\t\
		%s26.%s SECOND newly introduced allele in deme 0 at generation generationNewAllele2 in a single individual (0 = S; 1 = s; 2 = M; 3 = m), this argument is not considered if argument 19 is equal to 0.\n\n\t\
		%s27.%s generation at which ONE copy of the SECOND new allele is brought into the metapop. This argument is required but not considered if argument 19 is equal to 0.\n\n\t\
		%s28.%s Probability of homo-morphic pairing\n\n\t\
		%s29.%s 0: all males can be in the pool of reproductive males from which fathers will be sampled (fitness = 1 for all individuals); 1: frequency dependent sampling of a subset of malesfrom which fathers will be sampled (fitness = 1 - morphe_frequency)\n\n\t\
		%s30.%s 0: panmixia; if 1: trisStyle\n\n\
		%s\tExample:%s ./triStyli ${nDemes} ${K1} ${K2} ${pK1_to_K2} ${pK2_to_K1} ${nGenerations} ${nGenerations_extinction} ${nNtrlLoci} ${mu_ntrl} ${nQuantiLoci} ${mu_quanti} ${nBabies_max} ${prop_male_breeding} ${I} ${E} ${k} ${col_model} ${sexualSystem} ${sexEffects} ${selfingRate} ${verbose} ${seed} ${initial_situation} ${first_introduced_allele} ${time_first_introduction} ${second_introduced_allele} ${time_second_introduction} ${P_homo_pairing} ${male_choice} ${panmixia_or_trisStyly}\n\
		version: %s\n\t\tdependencies: \t%s\n", KRED, STOP, KRED, argc-1, STOP, KRED, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KMAG, STOP, KRED, STOP, VERSION, DEPENDENCY);

		exit(0);
	}
}

double fstMullon(const int maxIndPerDem_1, const double extinction, const int recolonization, const double migration){
	double res = 0.0;
	double numQ = 0.0;
	double denomQ = 0.0;
//	double phi = 0.0;
	double qr = 0.0;
	double migrationProportion = migration / maxIndPerDem_1;
	
	numQ = 1/(2.0 * maxIndPerDem_1) + extinction/(2.0 * recolonization) - extinction/(2.0 * recolonization * 2.0 * maxIndPerDem_1);
	denomQ = 1 - (1 - 1/(2.0 * maxIndPerDem_1)) * ((1 - migrationProportion)*(1 - migrationProportion) * (1 - extinction) + extinction * (1 - 1/(2.0 * recolonization)) * 1/(2.0 * recolonization -1));
//	phi = 1/(2.0 * recolonization -1);
	qr = numQ/denomQ;
	res = (qr - 1/(2.0 * maxIndPerDem_1)) * (2 * maxIndPerDem_1)/(2.0 * maxIndPerDem_1 -1);

	return(res);
}

double fstRousset(const int maxIndPerDem_1, const double extinction, const int recolonization, const double migration, const int colonizationModel){
	double res = 0.0;
	double numQ = 0.0;
	double denomQ = 0.0;
	double phi = 0.0;
	double qr = 0.0;
	double migrationProportion = migration / maxIndPerDem_1;
	
	if(colonizationModel == 1){
		phi = pow((1 - migrationProportion), 2);
	}

	numQ = 1/(2.0 * maxIndPerDem_1) + extinction/(2.0 * recolonization) - extinction/(2.0 * recolonization * 2.0 * maxIndPerDem_1);
	denomQ = 1 - (1 - 1/(2.0 * maxIndPerDem_1)) * (pow((1 - migrationProportion), 2) * (1 - extinction) + extinction * phi * (1 - 1/(2.0 * recolonization)) * 1/(2.0 * recolonization -1));

	qr = numQ/denomQ;

	res = (qr - 1/(2.0 * maxIndPerDem_1)) * (2 * maxIndPerDem_1)/(2.0 * maxIndPerDem_1 -1);

	return(res);
}

void sexInvador(gsl_rng* r, Deme* population, const int nDemes, const int* extinctionStatus, const double sexAvantage, const int sexualSystem, const double fecundity){
	int i = 0;
	for(i=0; i<nDemes; i++){
		if(extinctionStatus[i] != 1){
			if(population[i].nIndividus>1){
					population[i].sexChro[0] = 0; // add a heterogametous indivual which can be super male or female
					population[i].sexChro[1] = 1;
					population[i].sex[0] = 0;
        			        if(sexualSystem == 1){  // if XY system
        					population[i].femaleAllocation[0] = 0;
						population[i].maleAllocation[0] = sexAvantage;
						population[i].nOffsprings[0] = 0;
					}
					if(sexualSystem == 2){  // if ZW system
						population[i].femaleAllocation[0] = sexAvantage;
						population[i].maleAllocation[0] = 0;
						population[i].nOffsprings[0] = floor(fecundity * sexAvantage) + gsl_ran_binomial(r, (fecundity * sexAvantage) - floor(fecundity * sexAvantage), 1);
            			}
			}
		}
	}
}

/*# diveRsityR
#!/usr/bin/env Rscript
./diveRsity.R input=nameOfGenePopFile output=nameOfROutputFile
library(diveRsity)
options(warn=-1)
for(i in commandArgs()){
	tmp = strsplit(i, "=")
	if(tmp[[1]][1] == "input"){input = tmp[[1]][2]}
	if(tmp[[1]][1] == "output"){output = tmp[[1]][2]}
}

a=diffCalc(input, fst=T, pairwise=F, outfile=output)

*/


void nc(const Deme* population, const int nDemes, const int nNtrlLoci, const int locus, double* target){
	// fills the array target[nc, Hs, Htot]
	int i = 0;
	int j = 0;
	long allele1 = 0;
	long allele2 = 0;

	double* cont_table_tot = NULL;
	cont_table_tot = malloc((MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1)*sizeof(double));
	
	if(cont_table_tot == NULL){
		exit(0);
	}

	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		cont_table_tot[i] = 0;
	}	

	for(i=0; i<nDemes; i++){
		for(j=0; j<population[i].nIndividus; j++){
			allele1 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2];
			allele2 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2 + 1];
			
			cont_table_tot[allele1]++;
			cont_table_tot[allele2]++;
		}
	}

	// compute Ho and Hs
	double nInd = 0.0; // number of individuals in the deme_i
	double Ho = 0.0;
	double Ho_deme_i = 0.0;
	double Hs = 0.0;
	for(i=0; i<nDemes; i++){
		nInd = 0.0;
		Ho_deme_i = 0.0;
		double* cont_table_pop = NULL;
		cont_table_pop = malloc((MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1)*sizeof(double));
		
		if(cont_table_pop == NULL){
			exit(0);
		}
		
		for(j=0; j<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; j++){
			cont_table_pop[j]=0;
		}


		for(j=0; j<population[i].nIndividus; j++){
			nInd++;
			allele1 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2];
			allele2 = population[i].ntrlLoci[2*nNtrlLoci*j + locus*2 + 1];
			if(allele1 != allele2){
				Ho_deme_i++;
			}
				
			cont_table_pop[allele1]++;
			cont_table_pop[allele2]++;
		}
		Hs += heteroZ(cont_table_pop);
		free(cont_table_pop);
		Ho += Ho_deme_i/nInd;
	}
	Hs /= nDemes; // mean Hs
	Ho /= nDemes; // observed heterozygoty in the metapop

	// compute n, S1 and S2 (section 7.6 de la doc de la version 4.7 de genePop)
	int n = 0;
	int S1 = 0;
	double S2 = 0.0;
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		if(cont_table_tot[i] > 0){
			++n;
			S1 = S1 + cont_table_tot[i];
			S2 = S2 + pow(cont_table_tot[i], 2); 
		}
	}

	double Htot = 0.0;
	Htot = heteroZ(cont_table_tot);
	
	free(cont_table_tot);
	
	if( n<=1 ){
		target[0] = 0.0;
	}else{
		target[0] = (S1-S2/S1)/(n-1.0);
	}
	target[1] = Hs; // expected heterozygoty averaged over demes = 1 - sum( p_i^2 )
	target[2] = Htot; // expected heterozygoy in the whole metapopulation
	target[3] = Ho; // observed proportion of heterozygotes
}


double heteroZ(const double* cont_table){
	int i = 0;
	double n_ind = 0.0;
	double Hs = 1.0;
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		n_ind += 1.0 * cont_table[i];
	} 
	
	for(i=0; i<MAX_NUMBER_OF_INITIAL_NTRL_ALLELES+1; i++){
		Hs -= pow(cont_table[i]/n_ind, 2);
	} 

	return(Hs);	
	
}

void global_stat(Deme* population, const int nDemes, const long nNtrlLoci, double* diff_stats){
	double* nc_Hs_Htot_Ho = NULL;
	nc_Hs_Htot_Ho = malloc(4 * sizeof(double)); // contains [nc, Hs, Htot]
	if( nc_Hs_Htot_Ho == NULL ){
		exit(0);
	}

	int k = 0;
	int j = 0;
	int i = 0;
	int cnt = 0;
	double z_bar_j = 0.0;	
	
	double* z_bar = NULL;
	z_bar = malloc(nNtrlLoci * sizeof(double));
	if(z_bar == NULL){
		exit(0);
	}

	double* var_tot = NULL;
	var_tot = malloc(nNtrlLoci * sizeof(double));
	if(var_tot == NULL){
		exit(0);
	}

	double* var_among_patches = NULL;
	var_among_patches = malloc(nNtrlLoci * sizeof(double));
	if(var_among_patches == NULL){
		exit(0);
	}

	// global statistics
	double nc_k = 0.0;
	double Hs = 0.0;
	double Htot = 0.0;
	double Ho = 0.0;

	// global fst Charles
	double num_fst_cm = 0.0;
	double denom_fst_cm = 0.0;

	// global fst coal
	double num_fst_coal = 0.0;
	double denom_fst_coal = 0.0;

	// global g'st
	double num_gpst = 0.0;
	double denom_gpst = 0.0;

	// global Jost's D
	double num_D = 0.0;
	double denom_D = 0.0;

	// global Fis
	double num_Fis = 0.0;
	double denom_Fis = 0.0;
	
	// z_bar
	for(k=0; k<nNtrlLoci; k++){ // loop over loci k
		z_bar[k] = 0.0;
		cnt = 0;
		for(j=0; j<nDemes; j++){ // loop over demes j
			for(i=0; i<population[j].nIndividus; i++){ // loop over individuals i
				z_bar[k] += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 0]; // allele 1
				z_bar[k] += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 1]; // allele 2
				cnt += 2;
			}
		}
		z_bar[k] = z_bar[k]/cnt;
	}

	// variance among patches
	for(k=0; k<nNtrlLoci; k++){ // loop over loci k
		var_among_patches[k] = 0;

		cnt = 0;
		
		// compute z_bar_j
		for(j=0; j<nDemes; j++){ // loop over demes j
			z_bar_j = 0.0;
			cnt = 0;
			for(i=0; i<population[j].nIndividus; i++){ // loop over individuals i
				z_bar_j += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 0]; // allele 1
				z_bar_j += population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 1]; // allele 2
				cnt += 2;
			}
			z_bar_j /= cnt;
			
			var_among_patches[k] += pow(z_bar_j - z_bar[k], 2);
		}
		var_among_patches[k] /= nDemes;
		
	}

	// total variance in the population
	for(k=0; k<nNtrlLoci; k++){ // loop over loci k
		nc_k = 0.0;
		Hs = 0.0; 
		Htot = 0.0;
		Ho = 0.0;
		for(j=0; j<3; j++){
			nc_Hs_Htot_Ho[j] = 0.0;
		}

		var_tot[k] = 0.0;
		cnt = 0;
		for(j=0; j<nDemes; j++){ // loop over demes j
			for(i=0; i<population[j].nIndividus; i++){ // loop over individuals i
				var_tot[k] += pow(population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 0] - z_bar[k], 2); // allele 1
				var_tot[k] += pow(population[j].ntrlLoci[2*nNtrlLoci*i + 2*k + 1] - z_bar[k], 2); // allele 2
				cnt += 2;
			}
		}
		var_tot[k] = var_tot[k]/cnt;

		nc(population, nDemes, nNtrlLoci, k, nc_Hs_Htot_Ho);
		nc_k = nc_Hs_Htot_Ho[0];
		Hs = nc_Hs_Htot_Ho[1];
		Htot = nc_Hs_Htot_Ho[2];
		Ho = nc_Hs_Htot_Ho[3];
	
		// global Fst charles
		num_fst_cm += var_among_patches[k] * nc_k;
		denom_fst_cm += var_tot[k] * nc_k;
		
		// global Fst coal
		num_fst_coal += (Htot - Hs) * nc_k;
		denom_fst_coal += Htot * nc_k;
	
		// global g'st
		if( Htot != 0 ){
			num_gpst += ((Htot - Hs)/Htot) * nc_k; 
			denom_gpst += ((nDemes-1)*(1-Hs)/(nDemes-1+Hs)) * nc_k;
		}

		// global Jost's D
		num_D += (Htot - Hs) * nDemes * nc_k;
		denom_D += (1 - Hs) * (nDemes-1) * nc_k;
		
		// global Fis
		if( Hs != 0 ){
			num_Fis += (Hs - Ho) * nc_k;
			denom_Fis += Hs * nc_k;
		}
	}

	double global_fst_cm = 0.0;
	if( denom_fst_cm == 0.0){
		global_fst_cm = -9;
	}else{
		global_fst_cm = num_fst_cm / denom_fst_cm;
	}
	diff_stats[0] = global_fst_cm;

	double global_fst_coal = 0.0;
	if( denom_fst_coal == 0.0){
		global_fst_coal = -9;
	}else{
		global_fst_coal = num_fst_coal / denom_fst_coal;
	}
	diff_stats[1] = global_fst_coal;

	double global_gpst = 0.0;
	if( denom_gpst == 0.0){
		global_gpst = -9;
	}else{
		global_gpst = num_gpst / denom_gpst;
	}
	diff_stats[2] = global_gpst;

	double global_D = 0.0;
	if( denom_D == 0.0){
		global_D = -9;
	}else{
		global_D = num_D / denom_D;
	}
	diff_stats[3] = global_D;

	double global_Fis = 0.0;
	if( denom_Fis == 0.0){
		global_Fis = -9;
	}else{
		global_Fis = num_Fis / denom_Fis;
	}
	diff_stats[4] = global_Fis;

//	printf("Fst_CM\tFst_coal\tGst_p\tJostD\n");
//	printf("%f\t%f\t%f\t%f\n", global_fst_cm, global_fst_coal, global_gpst, global_D); 

	// free memory
	//free(fst);
	free(z_bar);
	free(var_tot);
	free(var_among_patches);
	free(nc_Hs_Htot_Ho);
}

void printMorphes(Deme* population, const int nDemes, const int generation){
	// plot the genotypes at the S and M loci, as well as the corresponding morphe, for all individuals
	int i=0;
	int j=0;
	printf("generation\tdeme\tindividual\tS1\tS2\tM1\tM2\tmorphe\n");	
	for(i=0; i<nDemes; i++){ // loop over the nDemes
		for(j=0; j<population[i].nIndividus; j++){
			printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", generation, i, j, population[i].locusS[2*j], population[i].locusS[2*j+1], population[i].locusM[2*j], population[i].locusM[2*j+1], population[i].morphe[j]);
		}
	}
}

