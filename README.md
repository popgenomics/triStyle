# triStyle
trisStyle is a metapopulation simulator to study the evolution of a heterostyle reproduction system and neutral unlinked loci.  
Are allowed:  
	Different models of migration are allowed (migrant _*versus*_ propagule pool).  
	Different initial situations (isoplethy, monomorphic for short, mid or long morph).  
	Different migration rates, extinction rates and number of recolonizers.  
	Different selfing rates. 
	Temporal fluctuations of carrying capacities by specifying the carrying capacities K1, K2 and the probabilities of transition from K1 to K2, and from K2 to K1.  


To compil the code:  
gcc triStyli.c -L/usr/local/lib -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o triStyli  
  
The arguments are :  
	1. Number of demes (>0)
	2. Carrying capacity 1 (int) : all demes share this K1 at generation 0
	3. Carrying capacity 2 (int)
	4. Proba of transition from K1 to K2 (float)
	5. Proba of transition from K2 to K1 (float)
	6. Number of generations (>0)
	7. Number of generations when extinction starts to occur
	8. Number of neutral loci (>=0)
	9. Neutral mutation rate (in [0-1]))
	10. Number of quantitative loci (>0; USELESS FOR THE MOMENT, SET IT TO 1)
	11. Quantitative mutation rate (in [0-1]; USELESS FOR THE MOMENT, SET IT TO 0)
	12. Max number of offsprings per hermaphrodite (>0)
	13. Proportion of individuals reproducing through the male function (in [0, 1])
	14. Immigration rate (Poisson distributed; >=0)
	15. Extinction rate (Binomialy distributed; in [0-1])
	16. Number of individuals recolonizing an extincted deme (>0)
	17. Colonization model, 'migration pool' (=0) or 'propagule pool' (=1) models
	18. sexualSystem is equal to 0 if autosomal, equal to 1 if XY and equal to 2 if ZW (USELESS FOR THE MOMENT, SET IT TO 0)
	19. Sexual effects of heterogametic sex (if equal to 1.5 in XY system, males have a 50 percent advantage to sire available ovules). Required but neglected if sexualSystem == 0 (USELESS FOR THE MOMENT, SET IT TO 1)
	20. Selfing rate of hermaphrodites, fixed over time (in [0-1])
	21. frequency of statistics calculation, all X generations (positive integer)
	22. Seed for the random generator (>0)
	23. initial situation is the phenotypic state of the metapopulation at generation 0 (0: freq is 1/3 of each morphes. 1: only ss mm. 2: only ss MM. 3: only SS mm. 4: only SS MM)
	24. FIRST newly introduced allele in deme 0 at generation generationNewAllele in a single individual (0 = S; 1 = s; 2 = M; 3 = m), this argument is not considered if argument 19 is equal to 0.
	25. generation at which ONE copy of the FIRST new allele is brought into the metapop. This argument is required but not considered if argument 20 is equal to 0.
	26. SECOND newly introduced allele in deme 0 at generation generationNewAllele2 in a single individual (0 = S; 1 = s; 2 = M; 3 = m), this argument is not considered if argument 19 is equal to 0.
	27. generation at which ONE copy of the SECOND new allele is brought into the metapop. This argument is required but not considered if argument 19 is equal to 0.
	28. Probability of homo-morphic pairing
	29. 0: all males can be in the pool of reproductive males from which fathers will be sampled (fitness = 1 for all individuals); 1: frequency dependent sampling of a subset of malesfrom which fathers will be sampled (fitness = 1 - morphe_frequency)
	30. 0: panmixia; if 1: trisStyle



Exemple of command:  
./triStyli_fluctu 100   200 20 0.1 0.1   2000 2001   20 0.00001 1 0.00001   2 1   0.5 0 1   1   0 1   0   1   123   0   0   2001   2   2001   0.2   0 1  

