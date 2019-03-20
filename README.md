# triStyle
trisStyle is a metapopulation simulator to study the evolution of a heterostyle reproduction system and neutral unlinked loci.  
Are allowed:  
	+ Different models of migration are allowed (migrant __versus__ propagule pool).  
	+ Different initial situations (isoplethy, monomorphic for short, mid or long morph).  
	+ Different migration rates, extinction rates and number of recolonizers.  
	+ Different selfing rates. 
	+ Temporal fluctuations of carrying capacities by specifying the carrying capacities K1, K2 and the probabilities of transition from K1 to K2, and from K2 to K1.  


To compil the code:  
gcc triStyli.c -L/usr/local/lib -lgsl -lgslcblas -lm -Wall -Wextra -Wshadow -Werror -O3 -o triStyli  
  
Exemple of code:  
./triStyli_fluctu 100   200 20 0.1 0.1   2000 2001   20 0.00001 1 0.00001   2 1   0.5 0 1   1   0 1   0   1   123   0   0   2001   2   2001   0.2   0 1
