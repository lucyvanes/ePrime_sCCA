# sparse CCA analysis on ePrime MRI and 4-year outcome data
# June 2020
# Lucy Vanes


setwd("C:/Users/k1456336/Documents/KCL BMEIS/ePrime/data")

library(PMA)
library(caret)
library(dplyr)	

# Make some decisions
#=====================================================
which_psych <- "item" # "composite" or "item"
control_vars <- c("sex","age","gest_weeks") # which covariates to regress out before CCA
Nkeep <- 100		# number of modes (not more than number of features in a single set)
Nperm <- 10000		# number of permutations for significance testing
Nboot <- 500		# number of bootstraps to assess feature significance
thresh <- 0.1		# threshold of loading coefficient for feature to be meaningful

output_dir <- ("C:/Users/k1456336/Documents/KCL BMEIS/ePrime/data/CCA output/CCA_vars_weights/bootstrapped_features/item_psych_500_permute_500_bootstrap_thresh0p5")

# Prepare data
#=====================================================
dat <- read.csv("ePrime_merged_data2.csv", header=T)
dat$id <- factor(dat$id)
dat$sex <- factor(dat$sex)

for (i in 3:length(dat))
{
  dat[[i]] <- as.numeric(dat[[i]])
}


composite_psych_vars <- c("ecbq_neg_affect", "ecbq_surgency",  "ecbq_effort_control","emque_emo_contagion","emque_attent_others","emque_prosocial", "sdq_emo", "sdq_conduct", "sdq_adhd", "sdq_peer", "sdq_prosocial", "adhd_inattention", "adhd_hyper", "srs_soc_awareness", "srs_soc_cog", "srs_soc_comm", "srs_soc_mot", "srs_rrb")

item_vars <- names(dat)[25:188]

brain_vars <- names(dat)[193:328]

if (which_psych=="composite"){
	cca_dat_all <- dat[,c("id","sex","age","gest_weeks",composite_psych_vars, brain_vars)]
} else if (which_psych=="item"){
	cca_dat_all <- dat[,c("id","sex","age","gest_weeks",item_vars, brain_vars)]
}
cca_dat <- cca_dat_all[complete.cases(cca_dat_all),]


# regress out control variables
#=======================================
if (which_psych=="composite"){
	vars <- c(composite_psych_vars, brain_vars)
} else if (which_psych=="item"){
	vars <- c(item_vars, brain_vars)
}


for (v in vars){
	lm1 <- lm(paste(v, "~",paste(control_vars, collapse = "+")), cca_dat)
	cca_dat[v] <- lm1$resid
}

# Get data ready for CCA
#========================================
if (which_psych=="composite"){
	mat1 <- as.matrix(cca_dat[,composite_psych_vars])
} else if (which_psych=="item"){
	mat1 <- as.matrix(cca_dat[,item_vars])
}
mat2 <- as.matrix(cca_dat[,brain_vars])			  # 170*136 / 130*136


#===================================================
# 				run original CCA
#===================================================

# Find best parameters for CCA
#==========================================
perm.out <- CCA.permute(mat1,mat2,typex="standard",typez="standard")
print(perm.out)
plot(perm.out)

# run CCA using parameters from perm.out
#==========================================
out <- CCA(mat1,mat2,typex="standard",typez="standard",K=Nkeep, penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)


#=====================================================
# 	Permutation testing of significance of each mode
#=====================================================

# create permutations of mat2 (brain data)
mat2_perms <- vector(mode = "list", length = Nperm)
for (p in 1:Nperm){
	mat2_perms[[p]] <- mat2[sample(nrow(mat2)),]
}

# Permutation loop
#=========================================
# loop over Nperm, run CCA on mat1 and each mat2 permutation in mat2_perms
# save corr for each permutation
cca_perms <- vector(mode = "list", length = Nperm)
perms_mode1_cors <- numeric()
	
for (p in 1:Nperm){
	print("=====================================")
	print(p)
	print("=====================================")
	cca_perms[[p]] <-  CCA(mat1,mat2_perms[[p]],typex="standard",typez="standard",
						   K=Nkeep, penaltyx=perm.out$bestpenaltyx,
						   penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)

	perms_mode1_cors[p] <- cca_perms[[p]]$cors[1]
}


# Test significance of each mode
#=============================================
# N = number of correlations from the first mode of each (permuted) CCA that are greater than correlation of mode d from original CCA

pvals <- numeric()
for (d in 1:Nkeep){
	cca_cor <- out$cors[d]
	N <- sum(perms_mode1_cors > cca_cor)
	pvals[d] <- (1+N)/Nperm 
}


Ncca <- sum(pvals<0.05) 				# number of significant CCA modes
cca_sig_modes <- which(pvals < 0.05) 	# indices for significant modes

# Save interesting data
#=========================
U <- as.data.frame(out$u)
V <- as.data.frame(out$v)

U_sig <- U[,cca_sig_modes]				# significant U variates
V_sig <- V[,cca_sig_modes]				# significant V variates


U_vars <- vector(mode = "list", length = Ncca)
V_vars <- vector(mode = "list", length = Ncca)
for (n in 1:Ncca){
	U_vars[[n]] <- which(U[cca_sig_modes[n]] > 0 | U[cca_sig_modes[n]] < 0)
	V_vars[[n]] <- which(V[cca_sig_modes[n]] > 0 | V[cca_sig_modes[n]] < 0)
}

# print variables for each component
for (n in 1:Ncca){
	print(paste("Significant mode", n))
	print("==========================")
	print(names(as.data.frame(mat1))[U_vars[[n]]])
	print(names(as.data.frame(mat2))[V_vars[[n]]])
	print("", quote=FALSE)
}


#=========================================================================
# Resampling procedure to assess feature significance
#=========================================================================

# create bootstrap samples
#==========================
samples_twothirds <- createDataPartition(as.numeric(cca_dat$id), p=.667, times=Nboot)
samples_onethird <- vector(mode = "list", length = Nboot)
for (s in 1:Nboot){
	samples_onethird[[s]] <- sample(samples_twothirds[[s]], length(cca_dat$id)-length(samples_twothirds[[s]]), replace = TRUE, prob = NULL)
}

boot_samples_mat1 <- vector(mode = "list", length = Nboot)
boot_samples_mat2 <- vector(mode = "list", length = Nboot)
for (s in 1:Nboot){
	boot_samples_mat1[[s]] <- mat1[c(samples_twothirds[[s]], samples_onethird[[s]]), ] 
	boot_samples_mat2[[s]] <- mat2[c(samples_twothirds[[s]], samples_onethird[[s]]), ] 
}

# Bootstrap loop
#============================
boot_cca <- vector(mode = "list", length = Nboot)
for (s in 1:Nboot){
	print(paste("BOOTSTRAP", s, "RUNNING"))
	print("==========================================")
	boot_cca[[s]] <- CCA(boot_samples_mat1[[s]],boot_samples_mat2[[s]],typex="standard",typez="standard",K=Nkeep, penaltyx=perm.out$bestpenaltyx,penaltyz=perm.out$bestpenaltyz, v=perm.out$v.init)
}


# Matching procedure
#===============================

## prepare variables 
U_boot_sig<- vector(mode = "list", length = Nboot)
V_boot_sig<- vector(mode = "list", length = Nboot)
boot_distr_U <- vector(mode = "list", length = Nboot)

df_U1 <- as.data.frame(matrix(nrow=Nboot, ncol=dim(U_sig[1])))
df_V1 <- as.data.frame(matrix(nrow=Nboot, ncol=dim(V_sig[1])))

df_U <- vector(mode = "list", length = Ncca)
df_V <- vector(mode = "list", length = Ncca)

for (n in 1:Ncca){
	df_U[[n]] <- df_U1
	df_V[[n]] <- df_V1
}

# Loop to match variates
#========================
for (s in 1:Nboot){
	# Find variate that most closely matches significant variates from original CCA
	#==============================================================================
	closest.match.u <-apply(abs(cor(U_sig,boot_cca[[s]]$u)),1,which.max) 
	closest.match.u.cor <-apply(abs(cor(U_sig,boot_cca[[s]]$u)),1,max) 

	closest.match.v <-apply(abs(cor(V_sig,boot_cca[[s]]$v)),1,which.max)
	closest.match.v.cor <-apply(abs(cor(V_sig,boot_cca[[s]]$v)),1,max) 

	U_boot_sig[[s]] <- as.data.frame(matrix(nrow=dim(U_sig[1]), ncol=Ncca))
	V_boot_sig[[s]] <- as.data.frame(matrix(nrow=dim(V_sig[1]), ncol=Ncca))

	for (n in 1:Ncca){
		if (closest.match.u[n]==closest.match.v[n]){
			U_boot_sig[[s]][n] <- boot_cca[[s]]$u[,closest.match.u[n]]
			V_boot_sig[[s]][n] <- boot_cca[[s]]$v[,closest.match.v[n]]

		} else {
			print(paste("Bootstrap",s,": could not find consistently matching variate for canonical variate",n))
			U_boot_sig[[s]][n] <- NA
			V_boot_sig[[s]][n] <- NA
		}

	}

	# fix sign if it's flipped
	#===================================
	for (c in 1:Ncca){
		corr <- cor(U_sig[,c], U_boot_sig[[s]][,c])
		if (is.na(corr)){
		} else if (corr < 0){
			U_boot_sig[[s]][c] <- U_boot_sig[[s]][c] * -1
			V_boot_sig[[s]][c] <- V_boot_sig[[s]][c] * -1
		} else {
		}
	}

	# save bootstrap distribution for each feature, each variate in U
	#=================================================================
	for (v in 1:Ncca){
		for (f in 1:dim(U_sig)[1]){
			df_U[[v]][s,f] <- U_boot_sig[[s]][f,v]
		}
	}
	# save bootstrap distribution for each feature, each variate in V
	#=================================================================
	for (v in 1:Ncca){
		for (f in 1:dim(V_sig)[1]){
			df_V[[v]][s,f] <- V_boot_sig[[s]][f,v]
		}
	}
}



# Determine bootstrap 95% CI for each feature, each variate
#============================================================

# Psych
#==============
feature_loadings_U <- vector(mode = "list", length = Ncca)
for (v in 1:Ncca){
	feature_loadings_U[[v]] <- as.data.frame(matrix(nrow=0, ncol=5))
	names(feature_loadings_U[[v]]) <- c("sig_feature_name", "sig_feature_mean", "sig_feature_std", "sig_feature_left", "sig_feature_right")
	for (f in 1:dim(U_sig)[1]){
		mean <- mean(df_U[[v]][,f], na.rm=T)
		std <- sd(df_U[[v]][,f], na.rm=T)
		n <- sum(!is.na(df_U[[v]][,f])) # count only bootstraps for which a match was found
		
		error <- qnorm(0.975)*std/sqrt(n)
		left <- mean-error
		right <- mean+error
		if (abs(mean) > thresh){
			if (left > 0 | right < 0){
				print(paste("Canonical variate", v ,"FEATURE NUMBER", f, "95% CI does not include 0"))
				# save mean loadings
				sig_feature_name <- colnames(mat1)[f]
				sig_feature_mean <- mean
				sig_feature_std <- std
				sig_feature_left <- left
				sig_feature_right <- right
				df_f <- data.frame(sig_feature_name, sig_feature_mean, sig_feature_std, sig_feature_left, sig_feature_right)
				feature_loadings_U[[v]] <- rbind(feature_loadings_U[[v]], df_f)
			}
		}
	}
}


# brain
#=============
feature_loadings_V <- vector(mode = "list", length = Ncca)

for (v in 1:Ncca){
	feature_loadings_V[[v]] <- as.data.frame(matrix(nrow=0, ncol=5))
	names(feature_loadings_V[[v]]) <- c("sig_feature_name", "sig_feature_mean", "sig_feature_std", "sig_feature_left", "sig_feature_right")
	for (f in 1:dim(V_sig)[1]){
		mean <- mean(df_V[[v]][,f], na.rm=T)
		std <- sd(df_V[[v]][,f], na.rm=T)
		n <- sum(!is.na(df_V[[v]][,f])) # count only bootstraps for which a match was found

		error <- qnorm(0.975)*std/sqrt(n)
		left <- mean-error
		right <- mean+error
				if (abs(mean) > thresh){
					if (left > 0 | right < 0){
				print(paste("Canonical variate", v ,"FEATURE NUMBER", f, "95% CI does not include 0"))
				# save mean loadings
				sig_feature_name <- colnames(mat2)[f]
				sig_feature_mean <- mean
				sig_feature_std <- std
				sig_feature_left <- left
				sig_feature_right <- right
				df_f <- data.frame(sig_feature_name, sig_feature_mean, sig_feature_std, sig_feature_left, sig_feature_right)
				feature_loadings_V[[v]] <- rbind(feature_loadings_V[[v]], df_f)
					}
				}
	}
}

#==============================================
# Save data on significant features
#==============================================

setwd(output_dir)

cca_U_vars_weights <- vector(mode = "list", length = Ncca)
cca_V_vars_weights <- vector(mode = "list", length = Ncca)

for (n in 1:Ncca){
	cca_U_vars_weights[[n]] <- data.frame("Unames"=feature_loadings_U[[n]][,1], "Uweights"=feature_loadings_U[[n]][,2])
	U_n_filename <- paste("CCA_U_mode_", n, "_vars_weights.csv", sep="")
	write.csv(cca_U_vars_weights[[n]], U_n_filename, row.names=F, quote=F)

	cca_V_vars_weights[[n]] <- data.frame("Vnames"=feature_loadings_V[[n]][,1], "Vweights"=feature_loadings_V[[n]][,2])
	V_n_filename <- paste("CCA_V_mode_", n, "_vars_weights.csv", sep="")
	write.csv(cca_V_vars_weights[[n]], V_n_filename, row.names=F, quote=F)
}






