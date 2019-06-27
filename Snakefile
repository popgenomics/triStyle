nDemes = 1 
K0 = [10, 100, 1000, 10000]
K2 = 10000
pK1_to_K2 = 0 
pK2_to_K1 =  0
nGenerations = 30
nGenerations_extinction = 31
nNtrlLoci = 1
mu_ntrl =  0
nQuantiLoci = 1 
mu_quanti = 0
nBabies_max = 100
prop_male_breeding = 1
I = 0
E = 0
k = 10
col_model = 0
sexualSystem = 0 
sexEffects = 1
selfingRate = 0
verbose = 1
nIterations = 10
initial_situation = 0 
first_introduced_allele = 0
time_first_introduction = 1
second_introduced_allele = 0
time_second_introduction = 999 
P_homo_pairing = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
male_choice = 0
panmiixa_or_trisStylie = 1

ITERATIONS = range(nIterations)

rule targets:
	input:
		expand('N_{K}/output_{iteration}.txt', K=K0, iteration=ITERATIONS)
	shell:
		"""
		head -n 1 N_10/output_1.txt > res.txt
		tail -n 1 -q {input} >>res.txt
		"""

rule simulations:
	output:
		'N_{K}/output_{iteration}.txt'
	shell:
		"""
		triStyli {nDemes} {wildcards.K} {K2} {pK1_to_K2} {pK2_to_K1} {nGenerations} {nGenerations_extinction} {nNtrlLoci} {mu_ntrl} {nQuantiLoci} {mu_quanti} {nBabies_max} {prop_male_breeding} {I} {E} {k} {col_model} {sexualSystem} {sexEffects} {selfingRate} {verbose} {wildcards.iteration} {initial_situation} {first_introduced_allele} {time_first_introduction} {second_introduced_allele} {time_second_introduction} {P_homo_pairing} {male_choice} {panmiixa_or_trisStylie}
		"""

