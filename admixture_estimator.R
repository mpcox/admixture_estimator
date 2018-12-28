# admixture calculations

# Murray Cox
# 15 May 2009

# checked that still runs
# 16 October 2018, R v.3.5.1

# method of Chakraborty et al. 1992 [Am J Hum Genet 50:145-155]

# program outputs *ASIAN* admixture rates

# populations
# Southern Han, PNG Highlanders, Sumba/Rindi

# globals
sample_density <- 1000
iterations     <- 1000
rounding       <- 2

# parental frequencies
asian_n_raw  <- c(24, 28, 14, 16, 8, 18, 22, 12, 11, 10, 25, 25, 25, 21, 22, 24, 25, 24, 25, 25, 24, 25, 24, 25, 25, 25, 24, 50, 50, 50, 50, 44, 50, 50, 48, 48, 50, 50)
asian_f_raw  <- c(0.75, 1.00, 0.71, 0.75, 0.63, 0.89, 0.55, 1.00, 0.82, 0.70, 0.84, 0.76, 0.84, 1.00, 0.82, 1.00, 0.88, 0.92, 0.92, 0.88, 0.88, 0.92, 0.58, 0.92, 0.92, 0.68, 0.88, 0.88, 1.00, 0.74, 1.00, 0.98, 1.00, 0.80, 0.92, 0.63, 1.00, 0.84)
asian_f_raw  <- round(asian_f_raw + 0.000000001, digits=rounding)

papuan_n_raw <- c(16, 18, 18, 18, 9, 16, 18, 16, 9, 9, 6, 6, 6, 0, 6, 5, 6, 6, 5, 6, 6, 6, 6, 6, 6, 6, 6, 10, 12, 12, 12, 12, 12, 12, 10, 12, 12, 12)
papuan_f_raw <- c(0.19, 0.22, 0.11, 0.00, 1.00, 0.13, 0.17, 0.19, 0.00, 0.89, 0.17, 0.00, 0.00, 0, 0.00, 0.20, 0.00, 0.33, 0.20, 0.33, 0.00, 0.33, 1.00, 0.17, 0.00, 0.00, 0.00, 0.00, 0.08, 0.00, 0.08, 0.25, 0.50, 0.00, 0.00, 0.00, 0.08, 0.08)
papuan_f_raw <- round(papuan_f_raw + 0.000000001, digits=rounding)

# hybrid frequencies
hybrid_n_raw <- c(20, 20, 20, 20, 10, 20, 20, 20, 10, 10, 10, 10, 10, 10, 9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 18, 18, 20, 20, 20, 20, 20, 20, 18, 20, 20)
hybrid_f_raw <- c(0.8500, 0.8500, 0.6500, 0.6000, 0.5000, 0.8000, 0.5500, 0.7000, 0.7000, 0.7000, 0.8000, 0.8000, 0.7000, 0.8000, 0.4444, 0.9000, 0.4000, 0.8000, 0.5000, 0.8000, 0.8000, 1.0000, 0.5000, 0.8000, 1.0000, 0.6000, 0.7000, 0.8333, 0.8889, 0.5000, 0.8000, 0.6000, 0.7000, 0.7000, 0.6500, 0.7778, 0.7000, 0.6000)
hybrid_f_raw <- round(hybrid_f_raw + 0.000000001, digits=rounding)

# genomic compartment
# all
compartment <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# autosomes only
#compartment <- c(1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1)

# X chromosome only
#compartment <- c(0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)


### test case
#asian_n_raw <- c(10,10,10)
#asian_f_raw <- c(0,0,0)

#papuan_n_raw <- c(10,10,10)
#papuan_f_raw <- c(1,1,1)

#hybrid_n_raw <- c(10,10,10)
#hybrid_f_raw <- c(0.8, 0.8, 0.8)

#compartment <- c(1,1,1)
###


# generate working copies
asian_n_sampled <- asian_n_raw * compartment
asian_f_sampled <- asian_f_raw

papuan_n_sampled <- papuan_n_raw * compartment
papuan_f_sampled <- papuan_f_raw

hybrid_n_sampled <- hybrid_n_raw * compartment
hybrid_f_sampled <- hybrid_f_raw


# primary set of error checking
# remove SNPs with no individuals typed in any sample
if( length(which(asian_n_sampled == 0)) > 0 ){

	delete_a <- which(asian_n_sampled == 0)

	asian_n_sampled <- asian_n_sampled[-c(delete_a)]
	asian_f_sampled <- asian_f_sampled[-c(delete_a)]
	
	papuan_n_sampled <- papuan_n_sampled[-c(delete_a)]
	papuan_f_sampled <- papuan_f_sampled[-c(delete_a)]
	
	hybrid_n_sampled <- hybrid_n_sampled[-c(delete_a)]
	hybrid_f_sampled <- hybrid_f_sampled[-c(delete_a)]
}

if( length(which(papuan_n_sampled == 0)) > 0 ){

	delete_m <- which(papuan_n_sampled == 0)

	asian_n_sampled <- asian_n_sampled[-c(delete_m)]
	asian_f_sampled <- asian_f_sampled[-c(delete_m)]
	
	papuan_n_sampled <- papuan_n_sampled[-c(delete_m)]
	papuan_f_sampled <- papuan_f_sampled[-c(delete_m)]
	
	hybrid_n_sampled <- hybrid_n_sampled[-c(delete_m)]
	hybrid_f_sampled <- hybrid_f_sampled[-c(delete_m)]
}

if( length(which(hybrid_n_sampled == 0)) > 0 ){

	delete_h <- which(hybrid_n_sampled == 0)

	asian_n_sampled <- asian_n_sampled[-c(delete_h)]
	asian_f_sampled <- asian_f_sampled[-c(delete_h)]
	
	papuan_n_sampled <- papuan_n_sampled[-c(delete_h)]
	papuan_f_sampled <- papuan_f_sampled[-c(delete_h)]
	
	hybrid_n_sampled <- hybrid_n_sampled[-c(delete_h)]
	hybrid_f_sampled <- hybrid_f_sampled[-c(delete_h)]
}


# number of markers
marker_n <- length(hybrid_n_sampled)

# construct data frame
P <- data.frame(cbind(  1, 1, 1:marker_n))

# population data frame
P[,2] <- asian_f_sampled
P[,1] <- papuan_f_sampled
P[,3] <- hybrid_f_sampled

# error checking
# tweak SNPs fixed in hybrid
for(i in 1:marker_n){
	if(P[i,3]==0) P[i,3] <- 0.0000000001
}
for(i in 1:marker_n){
	if(P[i,3]==1) P[i,3] <- P[i,3] - 0.0000000001
}


# estimator of ancestry proportion, mu, across loci (equation 6)
# also gives correct answer for a single locus
mu <- function( freq.table ) {

	P <- freq.table
	
	upper <- 0
	for(l in 1:length(P[,1])){
		sum1  <- (P[l,2]-P[l,1])*(P[l,3]-P[l,1])/P[l,3]
		sum2  <- ((1-P[l,2])-(1-P[l,1]))*((1-P[l,3])-(1-P[l,1]))/(1-P[l,3])
		upper <- upper + sum1 + sum2
	}
	
	lower <- 0
	for(l in 1:length(P[,1])){
		sum1  <- ((P[l,2]-P[l,1])^2)/P[l,3]
		sum2  <- (((1-P[l,2])-(1-P[l,1]))^2)/(1-P[l,3])
		lower <- lower + sum1 + sum2
	}

	return(upper/lower)
}

# Asian admixture estimate at face value
estimated.admixture <- mu(P)


# iterations to obtain Monte Carlo confidence intervals

# re-generate working copies
asian_n <- asian_n_raw * compartment
asian_f <- asian_f_raw

papuan_n <- papuan_n_raw * compartment
papuan_f <- papuan_f_raw

hybrid_n <- hybrid_n_raw * compartment
hybrid_f <- hybrid_f_raw


# secondary set of error checking
# remove SNPs with no individuals typed in any sample
if( length(which(asian_n == 0)) > 0 ){

	delete_a <- which(asian_n == 0)

	asian_n <- asian_n[-c(delete_a)]
	asian_f <- asian_f[-c(delete_a)]
	
	papuan_n <- papuan_n[-c(delete_a)]
	papuan_f <- papuan_f[-c(delete_a)]
	
	hybrid_n <- hybrid_n[-c(delete_a)]
	hybrid_f <- hybrid_f[-c(delete_a)]
}

if( length(which(papuan_n == 0)) > 0 ){

	delete_m <- which(papuan_n == 0)

	asian_n <- asian_n[-c(delete_m)]
	asian_f <- asian_f[-c(delete_m)]
	
	papuan_n <- papuan_n[-c(delete_m)]
	papuan_f <- papuan_f[-c(delete_m)]
	
	hybrid_n <- hybrid_n[-c(delete_m)]
	hybrid_f <- hybrid_f[-c(delete_m)]
}

if( length(which(hybrid_n == 0)) > 0 ){

	delete_h <- which(hybrid_n == 0)

	asian_n <- asian_n[-c(delete_h)]
	asian_f <- asian_f[-c(delete_h)]
	
	papuan_n <- papuan_n[-c(delete_h)]
	papuan_f <- papuan_f[-c(delete_h)]
	
	hybrid_n <- hybrid_n[-c(delete_h)]
	hybrid_f <- hybrid_f[-c(delete_h)]
}


# generate allele frequency probability densities
# asian frequency sampling distribution
sam_f  <- asian_f
sam_n  <- asian_n
asian  <- vector()
range  <- seq(0, 1, by=0.01)
for(i in 1:length(sam_n)){

	this.snp <- vector()
	for(j in 1:length(range)){
		
		sample <- rbinom(sample_density, sam_n[i], range[j]) / sam_n[i]
		match  <- which(round(sample + 0.000000001, digits=rounding) == sam_f[i])
		# fixes IEC rounding protocol which is not implemented in Excel
		
		if( length(match) == 0 ) { next }
		
		this.snp <- c(this.snp, rep(range[j], length(match)) )
	}
	
	asian <- c(asian, list(this.snp))
}

# papuan frequency sampling distribution
sam_f  <- papuan_f
sam_n  <- papuan_n
papuan <- vector()
range  <- seq(0, 1, by=0.01)
for(i in 1:length(sam_n)){

	this.snp <- vector()
	for(j in 1:length(range)){
		
		sample <- rbinom(sample_density, sam_n[i], range[j]) / sam_n[i]
		match  <- which(round(sample + 0.000000001, digits=rounding) == sam_f[i])
		# fixes IEC rounding protocol which is not implemented in Excel
		
		if( length(match) == 0 ) { next }
		
		this.snp <- c(this.snp, rep(range[j], length(match)) )
	}
	
	papuan <- c(papuan, list(this.snp))
}

# hybrid frequency sampling distribution
sam_f  <- hybrid_f
sam_n  <- hybrid_n
hybrid  <- vector()
range  <- seq(0, 1, by=0.01)
for(i in 1:length(sam_n)){

	this.snp <- vector()
	for(j in 1:length(range)){
		
		sample <- rbinom(sample_density, sam_n[i], range[j]) / sam_n[i]
		match  <- which(round(sample + 0.000000001, digits=rounding) == sam_f[i])
		# fixes IEC rounding protocol which is not implemented in Excel
		
		if( length(match) == 0 ) { next }
		
		this.snp <- c(this.snp, rep(range[j], length(match)) )
	}
	
	hybrid <- c(hybrid, list(this.snp))
}


# run iterations
mu.values <- vector()
for(i in 1:iterations){
		
	# resample frequencies
	asn_f <- vector()
	pap_f <- vector()
	hyb_f <- vector()
		
	for(i in 1:marker_n){
		asn_f <- c(asn_f, sample(asian[[i]],  1, replace=T))
		pap_f <- c(pap_f, sample(papuan[[i]], 1, replace=T))
		hyb_f <- c(hyb_f, sample(hybrid[[i]], 1, replace=T))
	}
	
	# change fixed frequencies in hybrid sample
	for(j in 1:marker_n){
		if(hyb_f[j]==0) hyb_f[j] <- 0.0000000001
	}
	for(j in 1:marker_n){
		if(hyb_f[j]==1) hyb_f[j] <- hyb_f[j] - 0.0000000001
	}
	
	# create resampled data frame
	d <- data.frame(pap_f, asn_f, hyb_f)

	mu.values <- c(mu.values, mu(d))
}


# calculate summary statistics
mean.estimate <- mean(mu.values)
median.estimate <- quantile(mu.values, 0.5)
quantile.estimate <- quantile(mu.values, c(0.025,0.975))

# print summary of results
marker_n
mean.estimate
median.estimate
quantile.estimate

