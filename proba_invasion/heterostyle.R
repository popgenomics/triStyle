# prob : get all the probabilities to get a possible genotype composition
P_x_to_y = function(prevGen=c(0,2,1), nextGen=genotypes_pop, N=N, P=0.1){
	# returns the probability for each genotypes contained in 'nextGen' to be reached in one generation given the frequencies at the previous generation 'prevGen'
	# prevGen contains c(n_SS, n_Ss, n_ss) qui are the number of SS, Ss and ss genotypes in the population respectively
	# example of the returned object:
	#> head(tmp)
	# next_generation SS Ss ss         prob
	#          0_0_40  0  0 40 2.492015e-20
	#          0_1_39  0  1 39 2.084231e-18
	#          0_2_38  0  2 38 8.497979e-17
	#          0_3_37  0  3 37 2.250677e-15
	#          0_4_36  0  4 36 4.353013e-14
	#          0_5_35  0  5 35 6.553264e-13

	f_SS = prevGen[1]/N # frequency of SS genotypes
	f_Ss = prevGen[2]/N # frequency of Ss genotypes
	f_ss = prevGen[3]/N # frequency of ss genotypes
	next_generation = c()
	SS = c()
	Ss = c()
	ss = c()
	prob = c()
	for(i in 1:nrow(nextGen)){
		if(f_ss==0){ # if no ss
			f_SS_next = f_SS + (1/4)*f_Ss**2 # frequency of SS at the next generation
			f_Ss_next = f_Ss*((1/2)*f_Ss + f_SS) # frequency of Ss at the next generation
			f_ss_next = (1/4) * f_Ss**2 # frequency of ss at the next generation
		}else{
			if(f_ss==1){ # if only ss
				f_SS_next = 0
				f_Ss_next = 0
				f_ss_next = 1
			}else{ # if at least one S and one s allele
				r_L = P*f_ss/((P - 1)*f_ss + 1) # probability to produce seeds via the Long x Long crosses
				r_S = (1-f_ss) * P / ( (1-f_ss)*P + f_ss ) # probability to produce seeds via the Short x Short crosses
			
				f_S = f_Ss + f_SS # frequency of the Short morphe in the population ( also equal to 1 - f_ss )
				
				f_SS_next = (f_SS*r_S*f_SS/f_S) + (f_SS*r_S*f_Ss/f_S*1/2) + (f_Ss*r_S*f_SS/f_S*1/2) + (f_Ss*r_S*f_Ss/f_S*1/4)
				f_ss_next = (f_ss * r_L) + ( f_ss * f_Ss/f_S * (1-r_L) * (1/2) ) + ( f_Ss * (1-r_S) * (1/2) ) + ( f_Ss * f_Ss/f_S * r_S * 1/4 )
				f_Ss_next = (f_ss * f_SS/f_S * (1-r_L) ) + ( f_ss * f_Ss/f_S * (1-r_L) * 1/2) + ( f_SS * (1-r_S)) + ( f_SS * f_Ss/f_S * r_S * 1/2 ) + ( f_Ss * (1-r_S) * 1/2 ) + ( f_Ss * f_Ss/f_S * r_S * 1/2 ) + ( f_Ss * f_SS/f_S * r_S * 1/2 )
			}
		}
		
		# records the tested genotypic composition (i.e, '1_4_5' for N=10 as instance)
		next_generation = c(next_generation, paste(nextGen[i,], collapse="_"))
		# records the probability of the tested genotypic composition (just above)
		prob = c(prob, dmultinom(nextGen[i,], N, c(f_SS_next, f_Ss_next, f_ss_next))) # dmultinom( c(n_SS, n_Ss, n_ss), N, c(f'SS, f'Ss, f'ss) ) returns the probability of trnasition from a state c(n_SS, n_Ss, n_ss) to a state c(f'SS, f'Ss, f'ss)*N
		SS = c(SS, nextGen[i, 1]) 
		Ss = c(Ss, nextGen[i, 2])
		ss = c(ss, nextGen[i, 3])
	}
	res = data.frame(next_generation=next_generation, SS=SS, Ss=Ss, ss=ss, prob=prob)
	return(res)
}



##############################################################################
# GET THE ANALYTICAL EXPECTATIONS OF MORPHE FREQUENCIES ACCORDING TO N AND P #
##############################################################################
# recorded results
p_long_fixed_new_S = c()
p_short_fixed_new_S = c()
p_polyM_new_S = c()
p_long_fixed_new_s = c()
p_short_fixed_new_s = c()
p_polyM_new_s = c()
vec_N = c()
vec_P = c()

for(N in c(10, 40)){
	# generates all the possible genotype compositions for a population of size N
	genotypes_pop = NULL
	for(i in 0:N){
		for(j in 0:(N-i)){
			SS = i
			Ss = j
			ss = N-i-j
			genotypes_pop = rbind(genotypes_pop, c(SS, Ss, ss))
		}
	}

	colnames(genotypes_pop) = c('SS', 'Ss', 'ss')


	#for(P in 0:10/10){
	for(P in 0:50/50){
		cat(paste('N=',N, 'P=', P, '\n', sep='\t'))
		
		# generates a transition matrix from i (row) to j (col) where i and j are population's genotypic composition, i.e, 0_5_3 for 0 SS individual, 5 Ss individuals and 3 ss.
		mat = NULL
		for(i in 1:nrow(genotypes_pop)){
			tmp = P_x_to_y(prevGen = genotypes_pop[i,], nextGen = genotypes_pop, N=N, P=P)
			mat = rbind(mat, tmp$prob)
		}
		colnames(mat) = tmp[,1]
		rownames(mat) = tmp[,1]

		# initial situation for 2 cases: introduction of a dominant S allele (pop_S) and a recessive s allele (pop_s)
		pop_S = mat[paste(0, 1, N-1, sep='_'), ] # introduction of a S allele in a ss population
		pop_s = mat[paste(N-1, 1, 0, sep='_'), ] # introduction of a s allele in a SS population
		for(generations in 2:(6*N)){ # loop over generations
			pop_S = pop_S%*%mat
			pop_s = pop_s%*%mat
		} # end of loop over generations
		
		vec_N = c(vec_N, N) # record the used 'N' value (population size)
		vec_P = c(vec_P, P) # record the used 'P' value (permasivity to homomorphic pairing)
		
		# initial condition: new S allele in the ancestral ss population c(0, 1, N-1)
		P1 = pop_S[1, paste(0, 0, N, sep='_')]
		P2 = pop_S[1, paste(N, 0, 0, sep='_')]
		p_long_fixed_new_S = c(p_long_fixed_new_S, P1) # probability to fix the long morphe 6.N generations after the introduction of a dominant S allele
		p_short_fixed_new_S = c(p_short_fixed_new_S, P2) # probability to fix the short morphe 6.N generations after the introduction of a dominant S allele
		p_polyM_new_S = c(p_polyM_new_S, 1-P1-P2) # probability to maintain short/long polymorphism 6.N generations after the introduction of a dominant S allele
		
		# initial condition: new s allele in the ancestral SS population c(N-1, 1, 0)
		P1 = pop_s[, paste(0, 0, N, sep='_')]
		P2 = pop_s[, paste(N, 0, 0, sep='_')]
		p_long_fixed_new_s = c(p_long_fixed_new_s, P1) # probability to fix the long morphe 6.N generations after the introduction of a recessive s allele
		p_short_fixed_new_s = c(p_short_fixed_new_s, P2) # probability to fix the short morphe 6.N generations after the introduction of a recessive s allele
		p_polyM_new_s = c(p_polyM_new_s, 1-P1-P2)# probability to maintain short/long polymorphism 6.N generations after the introduction of a recessive s allele
	}
}

expectations = data.frame(N=vec_N, P=vec_P, P_fixation_short_s=p_short_fixed_new_s, P_fixation_long_s=p_long_fixed_new_s, P_polymorphism_s=p_polyM_new_s, P_fixation_short_S=p_short_fixed_new_S, P_fixation_long_S=p_long_fixed_new_S, P_polymorphism_S=p_polyM_new_S)

### plot
library(viridis)
plot3var_v2(x=vec_P, y=vec_N, z=p_polyM_new_S , couleurs=viridis(10), xlab="Compatibility with same morph (P)", ylab="population size (N)", zlab="short/long polymorphism", FUN="mean", zlim=c(0,1))
plot3var_v2(x=vec_P, y=vec_N, z=p_long_fixed_new_S, couleurs=viridis(10), xlab="Compatibility with same morph (P)", ylab="population size (N)", zlab="fixation of the ancestral long morph", FUN="mean", zlim=c(0,1))
plot3var_v2(x=vec_P, y=vec_N, z=p_short_fixed_new_S, couleurs=viridis(10), xlab="Compatibility with same morph (P)", ylab="population size (N)", zlab="fixation of the new short morph", FUN="mean", zlim=c(0,1))


##########################################################################
# get the same quantities but obtained from individual based simulations #
##########################################################################
new_S = read.table("res_S_allele_in_ssmm.txt", h=T)
new_s = read.table("res_s_allele_in_SSmm.txt", h=T)

S_N10 = s_N10 = S_N40 = s_N40 = NULL

for(P in unique(new_S$proba_homomorphic_pairing)){
	# introduction of S
	## N10
	y = new_S[which(new_S$nIndMaxPerDeme_1==10 & new_S$proba_homomorphic_pairing==P), ]
	S_N10 = rbind(S_N10, c(P, mean(y$"nDemes_short_fixed"), mean(y$"nDemes_long_fixed"), mean(y$"nDemes_mid_lostOnly"))) # P; prop_fixation of short morph; prop_fixation of long morph; prop short/long polymorphism
	
	## N40
	y = new_S[which(new_S$nIndMaxPerDeme_1==40 & new_S$proba_homomorphic_pairing==P), ]
	S_N40 = rbind(S_N40, c(P, mean(y$"nDemes_short_fixed"), mean(y$"nDemes_long_fixed"), mean(y$"nDemes_mid_lostOnly"))) # P; prop_fixation of short morph; prop_fixation of long morph; prop short/long polymorphism
	
	# introduction of s
	## N10
	y = new_s[which(new_s$nIndMaxPerDeme_1==10 & new_s$proba_homomorphic_pairing==P), ]
	s_N10 = rbind(s_N10, c(P, mean(y$"nDemes_short_fixed"), mean(y$"nDemes_long_fixed"), mean(y$"nDemes_mid_lostOnly"))) # P; prop_fixation of short morph; prop_fixation of long morph; prop short/long polymorphism
	
	## N40
	y = new_s[which(new_s$nIndMaxPerDeme_1==40 & new_s$proba_homomorphic_pairing==P), ]
	s_N40 = rbind(s_N40, c(P, mean(y$"nDemes_short_fixed"), mean(y$"nDemes_long_fixed"), mean(y$"nDemes_mid_lostOnly"))) # P; prop_fixation of short morph; prop_fixation of long morph; prop short/long polymorphism
	
}
header = c('P', 'P_fixation_short', 'P_fixation_long', 'P_polymorphism')
colnames(S_N10) = header
colnames(s_N10) = header
colnames(S_N40) = header
colnames(s_N40) = header

#####################################
# plot expectations and simulations #
#####################################
par(mfrow=c(2,2))
par(las=1)
cex.lab=1.5
cex.axis=1.5
cex.main=1.7
# short/long polymorphism
	## expected values
	plot(expectations$P[which(expectations$N==10)], expectations$P_polymorphism_S[which(expectations$N==10)], ylim=c(0,1), type='l', col=viridis(2)[1], lwd=3, cex.lab=cex.lab, cex.axis=cex.axis, xlab = 'P', ylab = 'probability after 6.N generations', main = "short/long polymorphism after invasion", cex.main=cex.main)
	lines(expectations$P[which(expectations$N==10)], expectations$P_polymorphism_s[which(expectations$N==10)], col=viridis(2)[2], lwd=3)

	lines(expectations$P[which(expectations$N==40)], expectations$P_polymorphism_S[which(expectations$N==40)], col=viridis(2)[1], lwd=3, lty=4)
	lines(expectations$P[which(expectations$N==40)], expectations$P_polymorphism_s[which(expectations$N==40)], col=viridis(2)[2], lwd=3, lty=4)

	## simulated values
	points(S_N10[, 1], S_N10[,4], pch=16, cex=2.5, col=viridis(2)[1])
	points(s_N10[, 1], s_N10[,4], pch=16, cex=2.5, col=viridis(2)[2])

	points(S_N40[, 1], S_N40[,4], pch=18, cex=2.5, col=viridis(2)[1])
	points(s_N40[, 1], s_N40[,4], pch=18, cex=2.5, col=viridis(2)[2])

	legend('topright', lwd=c(3,3,3,3), lty=c(1, 1, 4, 4), col=c(viridis(2)[1], viridis(2)[2], viridis(2)[1], viridis(2)[2]), bty='n', legend=c('dominant S (N=10)', 'recessive s (N=10)','dominant S (N=40)', 'recessive s (N=40)'), cex=1.5)

par(las=1)
# fixation of the new morph
	## expected values
	plot(expectations$P[which(expectations$N==10)], expectations$P_fixation_short_S[which(expectations$N==10)], ylim=c(0,1), type='l', col=viridis(2)[1], lwd=3, cex.lab=cex.lab, cex.axis=cex.axis, xlab = 'P', ylab = 'probability after 6.N generations', main = "fixation of the new morph", cex.main=cex.main)
	lines(expectations$P[which(expectations$N==10)], expectations$P_fixation_long_s[which(expectations$N==10)], col=viridis(2)[2], lwd=3)

	lines(expectations$P[which(expectations$N==40)], expectations$P_fixation_short_S[which(expectations$N==40)], col=viridis(2)[1], lwd=3, lty=4)
	lines(expectations$P[which(expectations$N==40)], expectations$P_fixation_long_s[which(expectations$N==40)], col=viridis(2)[2], lwd=3, lty=4)

	## simulated values
	points(S_N10[, 1], S_N10[,2], pch=16, cex=2.5, col=viridis(2)[1])
	points(s_N10[, 1], s_N10[,3], pch=16, cex=2.5, col=viridis(2)[2])

	points(S_N40[, 1], S_N40[,2], pch=18, cex=2.5, col=viridis(2)[1])
	points(s_N40[, 1], s_N40[,3], pch=18, cex=2.5, col=viridis(2)[2])


par(las=1)
# fixation of the ancestral morph
	## expected values
	plot(expectations$P[which(expectations$N==10)], expectations$P_fixation_long_S[which(expectations$N==10)], ylim=c(0,1), type='l', col=viridis(2)[1], lwd=3, cex.lab=cex.lab, cex.axis=cex.axis, xlab = 'P', ylab = 'probability after 6.N generations', main = "elimination of the new morph", cex.main=cex.main)
	lines(expectations$P[which(expectations$N==10)], expectations$P_fixation_short_s[which(expectations$N==10)], col=viridis(2)[2], lwd=3)

	lines(expectations$P[which(expectations$N==40)], expectations$P_fixation_long_S[which(expectations$N==40)], col=viridis(2)[1], lwd=3, lty=4)
	lines(expectations$P[which(expectations$N==40)], expectations$P_fixation_short_s[which(expectations$N==40)], col=viridis(2)[2], lwd=3, lty=4)

	## simulated values
	points(S_N10[, 1], S_N10[,3], pch=16, cex=2.5, col=viridis(2)[1])
	points(s_N10[, 1], s_N10[,2], pch=16, cex=2.5, col=viridis(2)[2])

	points(S_N40[, 1], S_N40[,3], pch=18, cex=2.5, col=viridis(2)[1])
	points(s_N40[, 1], s_N40[,2], pch=18, cex=2.5, col=viridis(2)[2])

dev.print(pdf, "invasion.pdf", bg="white")
dev.off()

