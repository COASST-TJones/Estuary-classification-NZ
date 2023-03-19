# Estuary Classification Analyses

################################################################################
# Libraries and data
library(vegan)
library(mclust)
library(clValid)
library(randomForest)
library(MASS)
library(rpart)
est.data <- read.csv(here::here("data","EEC_data.csv"))
est.data$est.class[est.data$est.class=="C"] <- "B"

################################################################################

################################################################################
# Data exploration
# identify original classification
# original factors
# P/V = tidal prism/total estuary volume
# R12/V = River inflow per tidal cycle/total volume
# SC = shoreline complexity
# CI = closure index
# EE = estuary elongation

names(est.data)

# "no"			- Estuary number         
# "est_name"   	- Estuary name
# "council"		- Council
# "EWA"			- Estuary water area at HW spring tide (EWA) (m2)         
# "CLA"			- Catchment land area (CLA) (km2)
# "p_IA"		- Percent intertidal area (%)
# "mean_depth"	- Mean depth (m)
# "STP"			- Spring tidal prism (STP) (m3)
# "R12"			- River input over tidal cycle 12.4 hr (R12) (m3)
# "TEV"    		- Total estuary volume at HW (TEV) (m3)   
# "MCI"			- Mouth closure index (MCI)
# "SCI"			- Shoreline complexity index (SCI)
# "R12_TEV"		- Ratio of River inflow to total volume
# "STP_TEV"		- Ratio of spring tidal prism to total volume
# "CLA_EWA" 	- Ratio of catchment land area to estuary water area at high tide
# "EE"			- Elongation Index

# Plot all of the initial variables against each other to examine parameter correlations
# Initial factors are 
# STP_TEV
# R12_TEV
# SCI
# MCI
# EE

# re-calculate R12_TEV and STP_TEV values as the datasheet contains zero's for these values
est.data$R12_TEV  <- est.data$R12/est.data$TEV
est.data$STP_TEV  <- est.data$STP/est.data$TEV

# Plot correlations
par(mfrow=c(4,4), mar=c(4,4,0,0), oma=c(1,2,1,2))
pairs(est.data[,c("STP_TEV","R12_TEV","SCI","MCI")])

# Examine histograms of variables under different transformations 
par(mfrow=c(2,2))
hist(est.data$STP_TEV,100, xlab="STP/TEV", main="No Transform")
hist(sqrt(est.data$STP_TEV),100, xlab="sqrt(STP/TEV)", 
     main="Square-root Transform")
hist(sqrt(sqrt(est.data$STP_TEV)),100, xlab="(STP/TEV)^(1/4)", 
     main="Fourth-root Transform")
hist(log(est.data$STP_TEV+1),100, xlab="log(STP/TEV + 1)", 
     main="log(X+1) Transform")

# Biplots of transformed variables
par(mfrow=c(4,4), mar=c(4,4,0,0), oma=c(1,2,1,2))

plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="(STP/TEV)^(1/4)", xlab="", main="")
text(x=0, y=0, label="(STP/TEV)^(1/4)")
plot(0, type="n", ylab="", xlab="", axes=F, main="")
plot(0, type="n", ylab="", xlab="", axes=F, main="")
plot(0, type="n", ylab="", xlab="", axes=F, main="")

plot(sqrt(sqrt(R12_TEV)) ~ sqrt(sqrt(STP_TEV)), data=est.data, col=1, pch=16, ylab="(R12/TEV)^(1/4)", xlab="")
plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="", main="")
text(x=0, y=0, label="(R12/TEV)^(1/4)")
plot(0, type="n", ylab="", xlab="", axes=F, main="")
plot(0, type="n", ylab="", xlab="", axes=F, main="")

plot(SCI ~ sqrt(sqrt(STP_TEV)), data=est.data, col=1, pch=16, ylab="SCI", xlab="")
plot(SCI ~ sqrt(sqrt(R12_TEV)), data=est.data, col=1, pch=16, ylab="", xlab="")
plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="", main="")
text(x=0, y=0, label="SCI")
plot(0, type="n", ylab="", xlab="", axes=F, main="")

plot(MCI ~ sqrt(sqrt(STP_TEV)), data=est.data, col=1, pch=16, ylab="MCI", xlab="(STP/TEV)^(1/4)")
plot(MCI ~ sqrt(sqrt(R12_TEV)), data=est.data, col=1, pch=16, ylab="", xlab="(R12_TEV)^(1/4)")
plot(MCI ~ SCI, data=est.data, col=1, pch=16, ylab="", xlab="SCI")
plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="MCI", main="")
text(x=0, y=0, label="MCI")

###############################################################################

###############################################################################
# PCA
# Run a principal components analysis to examine multivariate dispersion
names(est.data)

# create data copy with transforms
est.data.a1 <- est.data[,c("no","MCI","SCI","R12_TEV","STP_TEV")]
est.data.a1$R12_TEV <- sqrt(sqrt(est.data.a1$R12_TEV))
est.data.a1$STP_TEV <- sqrt(sqrt(est.data.a1$STP_TEV))

pc1 <- princomp(est.data.a1)

pc1.scores <- pc1$scores

par(mfrow=c(2,2), mar=c(4,4,1,1))

plot(Comp.1 ~ Comp.2, data=pc1.scores, ylab="PCA1", xlab="PCA2", pch=16)
plot(Comp.1 ~ Comp.3, data=pc1.scores, ylab="PCA1", xlab="PCA3", pch=16)
plot(Comp.2 ~ Comp.3, data=pc1.scores, ylab="PCA2", xlab="PCA3", pch=16)

###############################################################################

###############################################################################
# Clustering Type 1 - Hierarchical Cluster Analyses

# Create distance matrix based on Euclidean distances
dist.est.data.a1 <- dist(as.matrix(est.data.a1[,-1]), method="euclidean")
hc.est.data.a1 <- hclust(dist.est.data.a1, method="ward")

# Plot dendrogram
par(mfrow=c(1,1))
plot(hc.est.data.a1, xlab="Estuary_no")

# Examine the number of clusters as a function of height
height.vec <- (1:8000)/100
ngroups <- 0

for ( i in 1:8000) {
	ngroups[i] <- max(cutree(hc.est.data.a1, h=height.vec[i]))
}

plot(ngroups ~ height.vec, pch=16, ylab="Number of groups", xlab="Height")

# examine a few cluster scenarios
cut.k10 <- cutree(hc.est.data.a1, k=4)

# Function for cross-tabulation of categorisations based on original vs new
table.meth <- function(clust.dat, act.dat) {
	
	ncols <- length(unique(clust.dat))
	nrows <- length(unique(act.dat))
	out.tab <- array(NA, dim=c(nrows, ncols))
	
	for ( i in 1:ncols) {
		for(j in 1:nrows){
			out.tab[j,i] <- length(clust.dat[as.factor(clust.dat)==levels(as.factor(clust.dat))[i] & as.factor(act.dat)==levels(as.factor(act.dat))[j]])
		}	
	}
	rownames(out.tab) <- levels(as.factor(act.dat))
	colnames(out.tab) <- levels(as.factor(clust.dat))
	return(out.tab)
}

table.meth(clust.dat=cut.k10, act.dat=as.character(est.data$est.class))

# Plot the characteristics color coded by grouping scheme
par(mfrow=c(4,4), mar=c(4,4,0,0), oma=c(1,2,1,2))

plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="(STP/TEV)^(1/4)", xlab="", main="")
text(x=0, y=0, label="(STP/TEV)^(1/4)")
plot(0, type="n", ylab="", xlab="", axes=F, main="")
plot(0, type="n", ylab="", xlab="", axes=F, main="")
plot(0, type="n", ylab="", xlab="", axes=F, main="")

plot(sqrt(sqrt(R12_TEV)) ~ sqrt(sqrt(STP_TEV)), data=est.data, col=cut.k10, pch=16, ylab="(R12/TEV)^(1/4)", xlab="")
plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="", main="")
text(x=0, y=0, label="(R12/TEV)^(1/4)")
plot(0, type="n", ylab="", xlab="", axes=F, main="")
plot(0, type="n", ylab="", xlab="", axes=F, main="")

plot(SCI ~ sqrt(sqrt(STP_TEV)), data=est.data, col=cut.k10, pch=16, ylab="SCI", xlab="")
plot(SCI ~ sqrt(sqrt(R12_TEV)), data=est.data, col=cut.k10, pch=16, ylab="", xlab="")
plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="", main="")
text(x=0, y=0, label="SCI")
plot(0, type="n", ylab="", xlab="", axes=F, main="")

plot(MCI ~ sqrt(sqrt(STP_TEV)), data=est.data, col=cut.k10, pch=16, ylab="MCI", xlab="(STP/TEV)^(1/4)")
plot(MCI ~ sqrt(sqrt(R12_TEV)), data=est.data, col=cut.k10, pch=16, ylab="", xlab="(R12_TEV)^(1/4)")
plot(MCI ~ SCI, data=est.data, col=cut.k10, pch=16, ylab="", xlab="SCI")
plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="MCI", main="")
text(x=0, y=0, label="MCI")

################################################################################

################################################################################
# finding the optimal number of groups
# Can be assessed by using a leave-one out or v-fold cross validation scheme
valid.clust <- function (inp.data, ngroups=3, n.left=5, 
                         lim=1000, verb=F, test.type="DIST") {
	
	len <- dim(inp.data)[1]													# length of original data
	perms <- rbind(c(1:ngroups), allPerms(ngroups, control=how(maxperm=1e6)))
	nperms <- dim(perms)[1]
	
	perc.match <- array(NA, dim=c(lim))
	diff.struct <- array(NA, dim=c(lim))
	
	# perform whole dataset clustering
	dist.orig <- dist(as.matrix(inp.data[,-1]), method="euclidean")			# calculation of distance matrix
	hc.orig <- hclust(dist.orig, method="ward")								# cluster analysis
	cut.orig <- cutree(hc.orig, k=ngroups)									# Find original grouping labels
	
	if(verb==TRUE) {														# The verb==TRUE option enables the reporting of detailed information on individual misclassifications
		est.levels <- levels(as.factor(inp.data$no))						# defines a list of estuary labels
		err.frame <- array(0, dim=c(length(est.levels), 2*ngroups+1))		# creates an array for storing the number of individual classifications
		err.frame[,1] <- cut.orig											# place full dataset classifications into the first column
	}
	
	# **
	# Now loop to examine the cross-validated error
	for ( e in 1:lim) {
		# Split data into training and testing dataset
		id.left.out <- floor(runif(n.left,min=1,max=len+0.9999999))				# identify left out testing data
		left.out <- inp.data[id.left.out,]										# test dataset
		train <- inp.data[-id.left.out,]										# Training dataset
		# **
	
		# Build tree based on training dataset
		dist.train <- dist(as.matrix(train[,-1]), method="euclidean")			# calculation of distance matrix
		hc.train <- hclust(dist.train, method="ward")							# Perform cluster analyses
		cut.train <- cutree(hc.train, k=ngroups)								# make cut at k=ngroups
		# **
	
		# There is potential for the original labellings and new labellings to be mismatched
		# The following code tests to see whether the two can be aligned 
		confu <- array(NA, dim=c(ngroups, ngroups))								# Definition of confusion matrix between labellings
		
		for (i in 1:ngroups) {													# Loop over groups
			orig.names <- inp.data[cut.orig==i,1]								# Find list of estuary id's for the original labelling of group i
			for (j in 1:ngroups) {						
				new.names <- train[cut.train==j,1]								# Find list of estuary id's for the training data labelling of group i
				confu[i,j] <- sum(new.names %in% orig.names)					# Find the number of matching entities
			}
		}
		
		diag.test <- array(NA, dim=c(nperms))									# creation of vector for storing test statistics
		
		for(i in 1:nperms) {													# loop through all possible arrangements of columns corresponding to re-labelling the training data cluster id's
			test.mat <- array(NA, dim=c(ngroups,ngroups))						# the best matching arrangement is determined as the one that achieves the highest matrix diagonal
			for ( j in 1:ngroups){												# Loop over groups to create the test matrix
				test.mat[,j]<- confu[,perms[i,j]]
			}
			diag.test[i] <- sum(diag(test.mat))									# calculate the diagonal test statistic
		}
		new.labs <- perms[which(diag.test==max(diag.test))[1],]					# list of new labels
		
		cut.train.adj <- array(NA, dim=c(length(cut.train)))					# assign new labels to data
		
		for(i in 1:ngroups) {													# Loop over groups
			cut.train.adj[cut.train==i] <- new.labs[i]							# Create new labelling based on match between train and original labelling schemes
		}
		
		# **
		# Now we need to identify which groups the training data fall into 
		
		if(test.type=="VAR") {
			# This first method minimises the total within group variance
			var.orig <- 0
			for (i in 1:ngroups) {
			# evaluate within group variance without the test data
				var.orig[i] <- sum(var(train[cut.train.adj==i,-1]))
			}
	
			test.label <- 0
			var.check <- 0	
			for(j in 1:n.left) {													# Loop over left out entities
				for(i in 1:ngroups) {												# Loop over possible groups for this entity
					t.data <- rbind(train[cut.train.adj==i,-1], left.out[j,-1])		# place left.out into each group and evaluate sum of within group variances
					var.check[i] <- sum(var(t.data)) + sum(var.orig[-i])
				}
				test.label[j] <- which(var.check==min(var.check))					# Identify the group label as the one which minimises the sum of within group variances
			}
			
		} else {
		# Identify groups based on minimum distance to group centroids
		
			MCI.centro <- tapply(train$MCI, cut.train.adj, FUN=mean)
			SCI.centro <- tapply(train$SCI, cut.train.adj, FUN=mean)
			STP_TEV.centro <- tapply(train$STP_TEV, cut.train.adj, FUN=mean)
			R12_TEV.centro <- tapply(train$R12_TEV, cut.train.adj, FUN=mean)

			test.label <- array(NA, dim=c(n.left))
			for(j in 1:n.left) {													# Loop over left out entities
				dist.check <- array(NA, dim=c(ngroups))
				for(i in 1:ngroups) {												# Loop over possible groups for this entity
					dist.check[i] <- (MCI.centro[i] - left.out$MCI[j])^2 + (SCI.centro[i] - left.out$SCI[j])^2 + (STP_TEV.centro[i] - left.out$STP_TEV[j])^2 + (R12_TEV.centro[i] - left.out$R12_TEV[j])^2
				}

				test.label[j] <- which(dist.check==min(dist.check))					# Identify the group label as the one which minimises the sum of within group variances
			}
		}
		
		
		match.stat <- array(NA, dim=c(n.left))
		for (i in 1:n.left) {
			match.stat[i] <- if(test.label[i] == cut.orig[id.left.out[i]]) 1 else {0}		# Test to see if group labellings math between training and test labellings
		}	
	
		perc.match[e] <- 100*mean(match.stat)									# calculation of % classification accuracy
		diff.struct[e] <- 100*(sum(confu)-max(diag.test))/(sum(confu))			# calculation of the difference in classification structure
		
		if(verb==TRUE) {
		# record which points are being classified incorrectly and in different classes

			for(i in 1:length(est.levels)){															# loop over all of the labels
																									# training errors and classificaton errors are treated seperately
				cl.group <- cut.train.adj[train$no==est.levels[i]]									# determine the cluster that this point was placed in
				if(length(cl.group)!=0) {																# (only when this point was in the training set)
					err.frame[i,cl.group+1] <- err.frame[i,cl.group+1] + 1							# identifies the column of the error array and adds one to the respective column
				}
				ne.group <- test.label[left.out$no==est.levels[i]]									# determine the cluster that the point was placed into
				if(length(ne.group)!=0) {																# only when this point was in the test dataset
					err.frame[i,ne.group+ngroups+1] <- err.frame[i,ne.group+ngroups+1] + 1			# identifies the column of the error array and adds one to the respective column
				}
			}		
		}
		
	}
	
	par(mfrow=c(2,1))
	hist(perc.match, 50)
	hist(diff.struct, 50)
	out.dat <- data.frame(perc=mean(perc.match, na.rm=T), dif=mean(diff.struct, na.rm=T))			# combine % mismatch error and % difference in cluster structure into output dataframe
	if(verb==TRUE) {
	
		out.dat <- data.frame(inp.data, err.frame)									# combine input data with error array
		names.add <- "P"
		
		for ( i in 1:ngroups) {														# create list of names for data frame
			names.add[i] <- paste("Trainclass_", i, sep="")
			names.add[i+ngroups] <- paste("Testclass_", i, sep="")
		}
		names(out.dat) <- c(names(inp.data), "Fullclass", names.add)
		
		# calculate overall test and training errors for each point
		
		Train.err <- array(NA, dim=c(dim(inp.data)[1]))								# Define arrays to fill for the overall train and test errors
		Test.err <- array(NA, dim=c(dim(inp.data)[1]))									
		
		sum.train <- rowSums(err.frame[,2:(ngroups+1)]) 								# calculate sum of times each point was in training and test sets
		sum.test <- rowSums(err.frame[,(ngroups+2):((2*ngroups)+1)]) 
		
		for( i in 1:ngroups) {
			for(j in 1:length(est.levels)){
				Train.err[out.dat$Fullclass==i] <- err.frame[out.dat$Fullclass==i,i+1]/sum.train[out.dat$Fullclass==i]
				Test.err[out.dat$Fullclass==i] <- err.frame[out.dat$Fullclass==i,ngroups+i+1]/sum.test[out.dat$Fullclass==i]
			}
		}
		
		out.dat <- data.frame(out.dat, Train.err, Test.err)
		
		return(out.dat)
	}
	return(out.dat)
}

# Evaluate this for different amount of omitted data
fold.vec <- c(10,20,30,40)
group.vec <- c(2,3,4,5,6,7,8,9)

error.vec <- 0
diff.vec <- 0
fold.out <- 0
group.out <- 0

K <- 1
for ( i in 1:4) {
	for( j in 1:8) {
		vc <- valid.clust(inp.data=est.data.a1, ngroups=group.vec[j],n.left=fold.vec[i], lim=500)
		error.vec[K] <- vc$perc
		diff.vec[K] <- vc$dif
		fold.out[K] <- fold.vec[i]
		group.out[K] <- group.vec[j]
		K <- K+1
	}
}

cv.dataframe <- data.frame(fold.out, group.out, error.vec, diff.vec)

# Plot allocation success as a function of the number of groups 
# acording to different amounts of omitted data
par(mfrow=c(2,1), mar=c(5,4,1,2))
plot(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==10,],
     pch=16, col=1, ylim=c(0,100), ylab="% classification error", 
     xlab="", type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], 
      pch=16, col=2, type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==30,], 
      pch=16, col=3, type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], 
      pch=16, col=4, type="b")

X <- c(2,3,4,5,6,7,8,9)
Y <- 100/X

lines(Y ~ X, col=1, lty=2, lwd=2)

# Plot 
plot(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==10,], pch=16, col=1, ylim=c(0,50), ylab="% change in structure", xlab="Number of Groups", type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=2, type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==30,], pch=16, col=3, type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=4, type="b")

legend("topleft", col=c(1,2,3,4,1), lty=c(1,1,1,1,2), legend=c("10 left out", "20 left out", "30 left out", "40 left out", "Expected by Chance"), bty="n")

# calculate cohens kappa

coh.kappa <- (cv.dataframe$error.vec - (100/cv.dataframe$group.out))/(100-(100/cv.dataframe$group.out))
cv.dataframe <- data.frame(cv.dataframe, coh.kappa)

plot(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==10,], pch=16, col=1, ylim=c(0,1), ylab="Cohens Kappa", xlab="", type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=2, type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==30,], pch=16, col=3, type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=4, type="b")

abline(h=0.2, col=1, lty=2, lwd=1)
abline(h=0.4, col=2, lty=2, lwd=1)
abline(h=0.6, col=3, lty=2, lwd=1)
abline(h=0.8, col=4, lty=2, lwd=1)

legend("bottomleft", col=c(1,2,3,4), lty=c(1,1,1,1), legend=c("10 left out", "20 left out", "30 left out", "40 left out"), bty="n")
legend("bottomright", col=c(1,2,3,4), lty=c(2,2,2,2), legend=c("Slight", "Fair", "Substantial", "Good"), bty="n")


# 2 - 5 groups is a reasonable number of clusters
plotfunc <- function(colscheme, leglab, colo) {

	par(mfrow=c(4,4), mar=c(4,4,0,0), oma=c(1,2,1,2), xpd=TRUE)

	plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="(STP/TEV)^(1/4)", xlab="", main="")
	text(x=0, y=0, label="(STP/TEV)^(1/4)")
	plot(0, type="n", ylab="", xlab="", axes=F, main="")
	plot(0, type="n", ylab="", xlab="", axes=F, main="")
	plot(0, type="n", ylab="", xlab="", axes=F, main="")

	plot(sqrt(sqrt(R12_TEV)) ~ sqrt(sqrt(STP_TEV)), data=est.data, col=colscheme, pch=16, ylab="(R12/TEV)^(1/4)", xlab="")
	plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="", main="")
	text(x=0, y=0, label="(R12/TEV)^(1/4)")
	plot(0, type="n", ylab="", xlab="", axes=F, main="")

	legend("center", col=colo, legend=leglab, pch=16, bty="n", ncol=2)

	plot(0, type="n", ylab="", xlab="", axes=F, main="")

	plot(SCI ~ sqrt(sqrt(STP_TEV)), data=est.data, col=colscheme, pch=16, ylab="SCI", xlab="")
	plot(SCI ~ sqrt(sqrt(R12_TEV)), data=est.data, col=colscheme, pch=16, ylab="", xlab="")
	plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="", main="")
	text(x=0, y=0, label="SCI")
	plot(0, type="n", ylab="", xlab="", axes=F, main="")

	plot(MCI ~ sqrt(sqrt(STP_TEV)), data=est.data, col=colscheme, pch=16, ylab="MCI", xlab="(STP/TEV)^(1/4)")
	plot(MCI ~ sqrt(sqrt(R12_TEV)), data=est.data, col=colscheme, pch=16, ylab="", xlab="(R12_TEV)^(1/4)")
	plot(MCI ~ SCI, data=est.data, col=colscheme, pch=16, ylab="", xlab="SCI")
	plot(5, ylim=c(-1,1), xlim=c(-1,1), ylab="", xlab="MCI", main="")
	text(x=0, y=0, label="MCI")
}

plotdiff <- function (group.labs) {

	t1 <- tapply(est.data.a1$STP_TEV,group.labs, FUN=mean)
	t2 <- tapply(est.data.a1$R12_TEV, group.labs, FUN=mean)
	t3 <- tapply(est.data.a1$SCI, group.labs, FUN=mean)
	t4 <- tapply(est.data.a1$MCI, group.labs, FUN=mean)

	s1 <- tapply(est.data.a1$STP_TEV,group.labs, FUN=sd)
	s2 <- tapply(est.data.a1$R12_TEV, group.labs, FUN=sd)
	s3 <- tapply(est.data.a1$SCI, group.labs, FUN=sd)
	s4 <- tapply(est.data.a1$MCI, group.labs, FUN=sd)
		
	par(mfrow=c(2,2), mar=c(3,5,2,1))
	
	plot(t1, pch=16, ylab="STP/TEV", xlab="", ylim=c(min(t1-s1), max(t1+s1)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t1[i]-s1[i], x1=i, y1=t1[i]+s1[i], col=1, lwd=2)
	}
	plot(t2, pch=16, ylab="R12/TEV", xlab="", ylim=c(min(t2-s2), max(t2+s2)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t2[i]-s2[i], x1=i, y1=t2[i]+s2[i], col=1, lwd=2)
	}
	plot(t3, pch=16, ylab="SCI", xlab="", ylim=c(min(t3-s3), max(t3+s3)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t3[i]-s3[i], x1=i, y1=t3[i]+s3[i], col=1, lwd=2)
	}
	plot(t4, pch=16, ylab="MCI", xlab="", ylim=c(min(t4-s4), max(t4+s4)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t4[i]-s4[i], x1=i, y1=t4[i]+s4[i], col=1, lwd=2)
	}
}

# 5 clusters

hc.cut.k5 <- cutree(hc.est.data.a1, k=5)

plotfunc(colscheme=hc.cut.k5, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"), colo=c(1,2,3,4,5))

plotdiff(hc.cut.k5)

table.meth(clust.dat=hc.cut.k5, act.dat=as.character(est.data$est.class))

vc <- valid.clust(inp.data=est.data.a1, ngroups=5,n.left=40, lim=1000, verb=TRUE)

Train.col <- 1 - vc$Train.err
Test.col <- 1 - vc$Test.err
Train.col <- topo.colors(10)[round(Train.col*9)+1]
Test.col <- topo.colors(10)[round(Test.col*9)+1]

plotfunc(colscheme=Train.col, leglab=c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60","60-70", "70-80", "80-90", "90-100"), colo=topo.colors(10))
plotfunc(colscheme=Test.col, leglab=c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60","60-70", "70-80", "80-90", "90-100"), colo=topo.colors(10))

# 6 clusters

hc.cut.k6 <- cutree(hc.est.data.a1, k=6)

plotfunc(colscheme=hc.cut.k6, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6"), colo=c(1,2,3,4,5,6))

plotdiff(hc.cut.k6)

table.meth(clust.dat=hc.cut.k6, act.dat=as.character(est.data$est.class))

###############################################################################

###############################################################################
# Apply the same but to scaled variables
# scale all four axes 
est.data.a2 <- est.data.a1
est.data.a2$R12_TEV <- scale(est.data.a2$R12_TEV)
est.data.a2$STP_TEV <- scale(est.data.a2$STP_TEV)
est.data.a2$SCI <- scale(est.data.a2$SCI)
est.data.a2$MCI <- scale(est.data.a2$MCI)

head(est.data.a2)

dist.est.data.a2 <- dist(as.matrix(est.data.a2[,-1]), method="euclidean")
hc.est.data.a2 <- hclust(dist.est.data.a2, method="ward")
plot(hc.est.data.a2, xlab="Estuary_no")

# identify the number of clusters based on allocation success
fold.vec <- c(10,20,30,40)
group.vec <- c(2,3,4,5,6,7,8)

error.vec <- 0
diff.vec <- 0
fold.out <- 0
group.out <- 0

K <- 1
for ( i in 1:4) {
	for( j in 1:7) {
		vc <- valid.clust(inp.data=est.data.a2, ngroups=group.vec[j],n.left=fold.vec[i], lim=300, test.type="DIST")
		error.vec[K] <- vc$perc
		diff.vec[K] <- vc$dif
		fold.out[K] <- fold.vec[i]
		group.out[K] <- group.vec[j]
		K <- K+1
	}
}

cv.dataframe <- data.frame(fold.out, group.out, error.vec, diff.vec)

par(mfrow=c(2,1), mar=c(5,4,1,2))
plot(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==10,], pch=16, col=1, ylim=c(0,100), ylab="% classification error", xlab="", type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=2, type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==30,], pch=16, col=3, type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=4, type="b")

X <- c(2,3,4,5,6,7,8,9)
Y <- 100/X

lines(Y ~ X, col=1, lty=2, lwd=2)

plot(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==10,], pch=16, col=1, ylim=c(0,50), ylab="% change in structure", xlab="Number of Groups", type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=2, type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==30,], pch=16, col=3, type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=4, type="b")

legend("topleft", col=c(1,2,3,4,1), lty=c(1,1,1,1,2), legend=c("10 left out", "20 left out", "30 left out", "40 left out", "Expected by Chance"), bty="n")

# calculate cohens kappa

coh.kappa <- (cv.dataframe$error.vec - (100/cv.dataframe$group.out))/(100-(100/cv.dataframe$group.out))
cv.dataframe <- data.frame(cv.dataframe, coh.kappa)

plot(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==10,], pch=16, col=1, ylim=c(0,1), ylab="Cohens Kappa", xlab="", type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=2, type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==30,], pch=16, col=3, type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=4, type="b")

abline(h=0.2, col=1, lty=2, lwd=1)
abline(h=0.4, col=2, lty=2, lwd=1)
abline(h=0.6, col=3, lty=2, lwd=1)
abline(h=0.8, col=4, lty=2, lwd=1)

legend("bottomleft", col=c(1,2,3,4), lty=c(1,1,1,1), legend=c("10 left out", "20 left out", "30 left out", "40 left out"), bty="n")
legend("bottomright", col=c(1,2,3,4), lty=c(2,2,2,2), legend=c("Slight", "Fair", "Substantial", "Good"), bty="n")

# 4 clusters

hc.cut.k4 <- cutree(hc.est.data.a2, k=4)

plotfunc(colscheme=hc.cut.k4, leglab=c("Group 1", "Group 2", "Group 3", "Group 4"), colo=c(1,2,3,4))

plotdiff(hc.cut.k4)

table.meth(clust.dat=hc.cut.k4, act.dat=as.character(est.data$est.class))

vc <- valid.clust(inp.data=est.data.a2, ngroups=4,n.left=40, lim=1000, verb=TRUE)

Train.col <- 1 - vc$Train.err
Test.col <- 1 - vc$Test.err
Train.col <- topo.colors(10)[round(Train.col*9)+1]
Test.col <- topo.colors(10)[round(Test.col*9)+1]

plotfunc(colscheme=Train.col, leglab=c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60","60-70", "70-80", "80-90", "90-100"), colo=topo.colors(10))
plotfunc(colscheme=Test.col, leglab=c("0-10", "10-20", "20-30", "30-40", "40-50", "50-60","60-70", "70-80", "80-90", "90-100"), colo=topo.colors(10))

###########################################################################################################################################
# Random Forest Clustering
###########################################################################################################################################

# random forests determines a distance matrix among points based on the 
# proportion of times those two data points are assigned to the same group
rf.data.a1 <- randomForest(est.data.a1[, -1], 
                           type="unsupervised", proximity=TRUE)

dist.rf.a1 <- as.dist(sqrt(1-rf.data.a1$proximity))

# this distance can then be placed into a hierarchical cluster analysis 
# scenario

hc.rf.a1 <- hclust(dist.rf.a1, method="ward")

# how to decide on the most appropriate number of clusters?
# For a given classification scheme (i.e. the results of selecting k=2 above) 
# we have a set of candidate labels
# Previously these were placed into groups based on euclidean distance to 
# group centroids
# IN this case the groups are not defined through euclidean distances, but 
# really upon random Forest dissimilarities
# One potential method is to use the random Forest method in a classification
# analysis using the candidate labels as the classes
# We can then record the out of bag error rate for this candidate set as a 
# measure of group strength

rf.clust.ident <- function(data.obj,tree.obj, kclust=c(2:10), nrepeats=50) {

	nclusts <- length(kclust)
	oob.error <- array(NA, dim=c(nclusts))
	oob.min <- array(NA, dim=c(nclusts))
	oob.max <- array(NA, dim=c(nclusts))
	
	for (i in 1:nclusts) {
		
		# make cut in the tree to identify the candidate group labellings
		cut.tree <- as.factor(cutree(tree.obj, k=kclust[i]))
		rf.dat.test <- data.frame(data.obj, cut.tree)
		oob.rep <- array(NA, dim=c(nrepeats))
		for (j in 1:nrepeats) {
			# create random forest classification tree
			rf.cut.tree <- randomForest(cut.tree ~ SCI + MCI + STP_TEV + R12_TEV, data=rf.dat.test, sampsize=221)
			conf.mat <- rf.cut.tree$confusion[,1:kclust[i]]
			oob.rep[j] <- 100 * (1- (sum(diag(conf.mat))/sum(conf.mat)))
		}
		oob.error[i] <- mean(oob.rep)
		oob.min[i] <- min(oob.rep)
		oob.max[i] <- max(oob.rep)
	}
	
	plot(oob.error ~ kclust, xlab="Number of Groups", ylab="OOB error (%)", type="l", lwd=2, ylim=c(0, max(oob.max)))
	lines(oob.min ~ kclust, lty=2, col=2, lwd=1.5)
	lines(oob.max ~ kclust, lty=2, col=2, lwd=1.5)
	return(oob.error)
}

rf.clust.ident(data.obj=est.data.a1,tree.obj=hc.rf.a1, kclust=c(2:50), nrepeats=10) 

# 2 clusters

hc.cut.k2 <- cutree(hc.rf.a1, k=2)

plotfunc(colscheme=hc.cut.k2, leglab=c("Group 1", "Group 2"), colo=c(1,2))

plotdiff(hc.cut.k2)

table.meth(clust.dat=hc.cut.k2, act.dat=as.character(est.data$est.class))

# 3 clusters

hc.cut.k3 <- cutree(hc.rf.a1, k=3)

plotfunc(colscheme=hc.cut.k3, leglab=c("Group 1", "Group 2", "Group 3"), colo=c(1,2,3))

plotdiff(hc.cut.k3)

table.meth(clust.dat=hc.cut.k3, act.dat=as.character(est.data$est.class))

# 4 clusters

hc.cut.k4 <- cutree(hc.rf.a1, k=4)

plotfunc(colscheme=hc.cut.k4, leglab=c("Group 1", "Group 2", "Group 3", "Group 4"), colo=c(1,2,3,4))

plotdiff(hc.cut.k4)

table.meth(clust.dat=hc.cut.k4, act.dat=as.character(est.data$est.class))

# 5 clusters

hc.cut.k5 <- cutree(hc.rf.a1, k=5)

plotfunc(colscheme=hc.cut.k5, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"), colo=c(1,2,3,4,5))

plotdiff(hc.cut.k5)

table.meth(clust.dat=hc.cut.k5, act.dat=as.character(est.data$est.class))

# 6 clusters

hc.cut.k6 <- cutree(hc.rf.a1, k=6)

plotfunc(colscheme=hc.cut.k6, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6"), colo=c(1,2,3,4,5,6))

plotdiff(hc.cut.k6)

table.meth(clust.dat=hc.cut.k6, act.dat=as.character(est.data$est.class))

# 7 clusters

hc.cut.k7 <- cutree(hc.rf.a1, k=7)

plotfunc(colscheme=hc.cut.k7, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7"), colo=c(1,2,3,4,5,6,7))

plotdiff(hc.cut.k7)

table.meth(clust.dat=hc.cut.k7, act.dat=as.character(est.data$est.class))


# assessment of random forests invariance to transformations

rf.test.raw <- est.data[,c(1,11,12,13,14)]
rf.test.trans <- est.data.a2

# plot MDS representations of the two schemes

raw.dist <- dist(rf.test.raw[,-1])
mds.raw.dist <- isoMDS(raw.dist, k=2)
mds.raw.scores <- mds.raw.dist$points

trans.dist <- dist(rf.test.trans[,-1])
mds.trans.dist <- isoMDS(trans.dist, k=2)
mds.trans.scores <- mds.trans.dist$points

par(mfrow=c(2,1), mar=c(4,4,2,1))
plot(mds.raw.scores[,1] ~ mds.raw.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=est.data$est.class, cex=0.75, main="Raw Data")
legend("topleft", col=c(1,2,3,4,5,6,7), legend=c("A", "B",'D','E','F','G','H'), pch=16, bty="n")
plot(mds.trans.scores[,1] ~ mds.trans.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=est.data$est.class, cex=0.75, main="Transformed Data")

# Now examine the random forest distances through MDS

rf.raw <- randomForest(rf.test.raw[, -1], type="unsupervised", proximity=TRUE)
rf.raw.prox <- rf.raw$proximity
rf.raw.prox[415,436] <- 0.9999999
rf.raw.prox[436,415] <- 0.9999999
dist.rf.raw <- as.dist(sqrt(1-rf.raw.prox))
mds.raw.rf <- isoMDS(dist.rf.raw, k=2)
mds.raw.rf.scores <- mds.raw.rf$points

rf.trans <- randomForest(rf.test.trans[, -1], type="unsupervised", proximity=TRUE)
rf.trans.prox <- rf.trans$proximity
rf.trans.prox[415,436] <- 0.9999999
rf.trans.prox[436,415] <- 0.9999999
dist.rf.trans <- as.dist(sqrt(1-rf.trans.prox))
mds.trans.rf <- isoMDS(dist.rf.trans, k=2)
mds.trans.rf.scores <- mds.trans.rf$points

mds.trans.rf.scores[,2] <- -mds.trans.rf.scores[,2]

par(mfrow=c(2,1), mar=c(4,4,2,1))
plot(mds.raw.rf.scores[,1] ~ mds.raw.rf.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=est.data$est.class, cex=0.75, main="Raw Data")
plot(mds.trans.rf.scores[,1] ~ mds.trans.rf.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=est.data$est.class, cex=0.75, main="Transformed Data")

# Now examine whether the groupings that result are the same
rf.raw <- randomForest(rf.test.raw[, -1], type="unsupervised", proximity=TRUE)
rf.raw.prox <- rf.raw$proximity
dist.rf.raw <- as.dist(sqrt(1-rf.raw.prox))
rf.trans <- randomForest(rf.test.trans[, -1], type="unsupervised", proximity=TRUE)
rf.trans.prox <- rf.trans$proximity
dist.rf.trans <- as.dist(sqrt(1-rf.trans.prox))

hc.rf.raw <- hclust(dist.rf.raw, method="ward")
hc.rf.trans <- hclust(dist.rf.trans, method="ward")
 
ct.rf.raw <- cutree(hc.rf.raw, k=2)
ct.rf.trans <- cutree(hc.rf.trans, k=2)

table.meth(clust.dat=ct.rf.raw, act.dat=ct.rf.trans)

# Stochasticity and changes to structure



valid.forest <- function (inp.data, ngroups, lim=100, ntrees=c(100,200,300,400)) {
	
	error.tet <- 0
	error.min <- 0
	error.max <- 0
	
	rec.class <- array(NA, dim=c(dim(inp.data)[1], 2*lim))
		
	for(r in 1:length(ntrees)) {
	
		error.vec <- array(NA, dim=c(lim))
		perms <- rbind(c(1:ngroups), allPerms(ngroups, control=how(maxperm=1e6)))
		nperms <- dim(perms)[1]
	
		for(e in 1:lim) {
		
			rf.test <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntree=ntrees[r])
			rf.prox <- rf.test$proximity
			rf.dist <- as.dist(sqrt(1-rf.prox))
			hc.rf <- hclust(rf.dist, method="ward")
			class.scheme.1 <- cutree(hc.rf, k=ngroups)
		
			rf.test <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntree=ntrees[r])
			rf.prox <- rf.test$proximity
			rf.dist <- as.dist(sqrt(1-rf.prox))
			hc.rf <- hclust(rf.dist, method="ward")
			class.scheme.2 <- cutree(hc.rf, k=ngroups)
		
			# There is potential for the original labellings and new labellings to be mismatched
			# The following code tests to see whether the two can be aligned 
			confu <- array(NA, dim=c(ngroups, ngroups))								# Definition of confusion matrix between labellings
	
			for (i in 1:ngroups) {													# Loop over groups
				orig.names <- inp.data[class.scheme.1==i,1]								# Find list of estuary id's for the original labelling of group i
				for (j in 1:ngroups) {						
					new.names <- inp.data[class.scheme.2==j,1]								# Find list of estuary id's for the training data labelling of group i
					confu[i,j] <- sum(new.names %in% orig.names)					# Find the number of matching entities
				}
			}
		
			diag.test <- array(NA, dim=c(nperms))									# creation of vector for storing test statistics
		
			for(i in 1:nperms) {													# loop through all possible arrangements of columns corresponding to re-labelling the training data cluster id's
				test.mat <- array(NA, dim=c(ngroups,ngroups))						# the best matching arrangement is determined as the one that achieves the highest matrix diagonal
				for ( j in 1:ngroups){												# Loop over groups to create the test matrix
					test.mat[,j]<- confu[,perms[i,j]]
				}
				diag.test[i] <- sum(diag(test.mat))									# calculate the diagonal test statistic
			}
			new.labs <- perms[which(diag.test==max(diag.test))[1],]					# list of new labels
			class.scheme.1.adj <- array(NA, dim=c(length(class.scheme.1)))					# assign new labels to data
		
			for(i in 1:ngroups) {													# Loop over groups
				class.scheme.1.adj[class.scheme.1==i] <- new.labs[i]							# Create new labelling based on match between train and original labelling schemes
			}
			
			test.matrix <- table.meth(clust.dat=class.scheme.1.adj, act.dat=class.scheme.2)
			error.vec[e] <- 100*sum(diag(as.matrix(test.matrix)))/sum(test.matrix)
		}
	
		error.tet[r] <- mean(error.vec)
		error.min[r] <- error.vec[rank(error.vec, ties.method="random")== ceiling(lim*0.05)]
		error.max[r] <- error.vec[rank(error.vec, ties.method="random")== ceiling(lim*0.95)]
	}
	
	plot(error.tet ~ ntrees, ylab="% Similarity among classifications", xlab="Number of Trees", ylim=c(0,100), type="l", lty=1, col=1, lwd=2)
	lines(error.min ~ ntrees, lty=2, type="l",lwd=1.5, col=2) 
	lines(error.max ~ ntrees, lty=2, type="l",lwd=1.5, col=2) 
	
	error.dat <- data.frame(ntrees, error.tet, error.min, error.max)
	return(error.dat)
}	
	
X <- 50*c(1:10)
Y <- 1000*c(1:10)
	
treenum <- c(X,Y)
vf.k4 <- valid.forest(est.data.a1, ngroups=4, lim=100, ntrees=treenum)



valid.forest.check <- function (inp.data, ngroups, lim=100, ntrees=2000) {

	rec.class <- array(NA, dim=c(dim(inp.data)[1], 2*lim))

	perms <- rbind(c(1:ngroups), allPerms(ngroups, control=how(maxperm=1e6)))
	nperms <- dim(perms)[1]

	for(e in 1:lim) {
	
		rf.test <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntree=ntrees)
		rf.prox <- rf.test$proximity
		rf.dist <- as.dist(sqrt(1-rf.prox))
		hc.rf <- hclust(rf.dist, method="ward")
		class.scheme.1 <- cutree(hc.rf, k=ngroups)
		
		rf.test <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntree=ntrees)
		rf.prox <- rf.test$proximity
		rf.dist <- as.dist(sqrt(1-rf.prox))
		hc.rf <- hclust(rf.dist, method="ward")
		class.scheme.2 <- cutree(hc.rf, k=ngroups)
		
		# There is potential for the original labellings and new labellings to be mismatched
		# The following code tests to see whether the two can be aligned 
		confu <- array(NA, dim=c(ngroups, ngroups))								# Definition of confusion matrix between labellings
	
		for (i in 1:ngroups) {													# Loop over groups
			orig.names <- inp.data[class.scheme.1==i,1]								# Find list of estuary id's for the original labelling of group i
			for (j in 1:ngroups) {						
				new.names <- inp.data[class.scheme.2==j,1]								# Find list of estuary id's for the training data labelling of group i
				confu[i,j] <- sum(new.names %in% orig.names)					# Find the number of matching entities
			}
		}
	
		diag.test <- array(NA, dim=c(nperms))									# creation of vector for storing test statistics
	
		for(i in 1:nperms) {													# loop through all possible arrangements of columns corresponding to re-labelling the training data cluster id's
			test.mat <- array(NA, dim=c(ngroups,ngroups))						# the best matching arrangement is determined as the one that achieves the highest matrix diagonal
			for ( j in 1:ngroups){												# Loop over groups to create the test matrix
				test.mat[,j]<- confu[,perms[i,j]]
			}
			diag.test[i] <- sum(diag(test.mat))									# calculate the diagonal test statistic
		}
		new.labs <- perms[which(diag.test==max(diag.test))[1],]					# list of new labels
		class.scheme.1.adj <- array(NA, dim=c(length(class.scheme.1)))					# assign new labels to data
		
		for(i in 1:ngroups) {													# Loop over groups
			class.scheme.1.adj[class.scheme.1==i] <- new.labs[i]							# Create new labelling based on match between train and original labelling schemes
		}
		
		rec.class[,((2*e)-1)] <- class.scheme.1.adj
		rec.class[,(2*e)] <- class.scheme.2
	}
	
	return(rec.class)
	
}
	
	
vf.k4 <- valid.forest.check(est.data.a1,ngroups=4, lim=5, ntrees=2000)

vf.k4.adj <- array(NA, dim=dim(vf.k4))
perms <- rbind(c(1:4), allPerms(4, control=how(maxperm=1e6)))
nperms <- dim(perms)[1]
confu <- array(NA, dim=c(4, 4))								
for (e in 2:200) {
	for (i in 1:4) {													
		orig.names <- est.data.a1[vf.k4[,e]==i,1]								
		for (j in 1:4) {						
			new.names <- est.data.a1[vf.k4[,1]==j,1]								
			confu[i,j] <- sum(new.names %in% orig.names)					
		}
	}
	
	diag.test <- array(NA, dim=c(nperms))	
								
	for(i in 1:nperms) {
												
		test.mat <- array(NA, dim=c(4,4))						
		for ( j in 1:4){												
			test.mat[,j]<- confu[,perms[i,j]]
		}
		diag.test[i] <- sum(diag(test.mat))									
	}
	new.labs <- perms[which(diag.test==max(diag.test))[1],]					
				
		
	for(i in 1:4) {													
		vf.k4.adj[vf.k4[,e]==i,e] <- new.labs[i]							
	}
}

nf <- layout(pl.mat,heights=c(2,1))
nf
plot(SCI ~ R12_TEV, data=est.data.a1, col=vf.k4.adj[,1], pch=16, cex=0.6)
lines(SCI ~ R12_TEV, data=est.data.a1[c(28,247,19,67,350,394),], col=vf.k4.adj[c(28,247,19,67,350,394),1], pch=16, cex=2, type="p")
text(x=est.data.a1$R12_TEV[c(28,247,19,67,350,394)], y=est.data.a1$SCI[c(28,247,19,67,350,394)], labels=c("A","B","C","D","E","F"), pos=2, cex=2)

plot(SCI ~ STP_TEV, data=est.data.a1, col=vf.k4.adj[,1], pch=16, cex=0.6)
lines(SCI ~ STP_TEV, data=est.data.a1[c(28,247,19,67,350,394),], col=vf.k4.adj[c(28,247,19,67,350,394),1], pch=16, cex=2, type="p")
text(x=est.data.a1$STP_TEV[c(28,247,19,67,350,394)], y=est.data.a1$SCI[c(28,247,19,67,350,394)], labels=c("A","B","C","D","E","F"), pos=2, cex=2)


thist1 <- vf.k4.adj[28,]
nom <- paste("A - ", est.data$est_name[28])
hist(thist1, breaks=c(0.5,1.5,2.5,3.5,4.5), col=c(1,2,3,4), main=nom, xlab="Group")

thist1 <- vf.k4.adj[247,]
nom <- paste("B - ", est.data$est_name[247])
hist(thist1, breaks=c(0.5,1.5,2.5,3.5,4.5), col=c(1,2,3,4), main=nom, xlab="Group")

thist1 <- vf.k4.adj[19,]
nom <- paste("C - ", est.data$est_name[19])
hist(thist1, breaks=c(0.5,1.5,2.5,3.5,4.5), col=c(1,2,3,4), main=nom, xlab="Group")

thist1 <- vf.k4.adj[67,]
nom <- paste("D - ", est.data$est_name[67])
hist(thist1, breaks=c(0.5,1.5,2.5,3.5,4.5), col=c(1,2,3,4), main=nom, xlab="Group")

thist1 <- vf.k4.adj[350,]
nom <- paste("E - ", est.data$est_name[350])
hist(thist1, breaks=c(0.5,1.5,2.5,3.5,4.5), col=c(1,2,3,4), main=nom, xlab="Group")

thist1 <- vf.k4.adj[394,]
nom <- paste("F - ", est.data$est_name[394])
hist(thist1, breaks=c(0.5,1.5,2.5,3.5,4.5), col=c(1,2,3,4), main=nom, xlab="Group")

##############################################################################################################################################
# Model Based Clustering
##############################################################################################################################################

# perform initial model-based cluster analysis

Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:12)

summary(Mclust.data.a1)
summary(Mclust.data.a1, parameters = TRUE)
plot(Mclust.data.a1)
Mclust.class <- Mclust.data.a1$classification

MclustBIC.data.a1 <- mclustBIC(est.data.a1[,-1], G=2:20)

# EEV 12 clusters

Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:12)
Mclust.class <- Mclust.data.a1$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9", "Group 10", "Group 11", "Group 12"), colo=c(1,2,3,4,5,6,7,8,9,10,11,12))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))

# VEV 9 clusters

Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:9)
Mclust.class <- Mclust.data.a1$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9"), colo=c(1,2,3,4,5,6,7,8,9))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))

# VVI 6 clusters

Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:6)
Mclust.class <- Mclust.data.a1$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6"), colo=c(1,2,3,4,5,6))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))

# Repeat the analyses but using scaled variables

Mclust.data.a2 <- Mclust(est.data.a2[,-1], G=2:20)
MclustBIC.data.a2 <- mclustBIC(est.data.a2[,-1], G=2:20)

summary(Mclust.data.a2)
summary(Mclust.data.a2, parameters = TRUE)
plot(Mclust.data.a2)
Mclust.class <- Mclust.data.a2$classification

# VEV 6 clusters

Mclust.data.a2 <- Mclust(est.data.a2[,-1], G=2:6)
Mclust.class <- Mclust.data.a2$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6"), colo=c(1,2,3,4,5,6))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))


# VEV 7 clusters

Mclust.data.a2 <- Mclust(est.data.a2[,-1], G=2:7)
Mclust.class <- Mclust.data.a2$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7"), colo=c(1,2,3,4,5,6,7))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))


# Generate decision rule trees for each classification

Tree.Construct <- function(class.dat, type="restricted") {

	len.vector <- tapply(class.dat, class.dat, FUN=length)
	buck <- round(min(len.vector)*0.5)
	mc.dat <- data.frame(est.data, class.dat)
	if(type=="restricted") {
		tree <- rpart(class.dat ~ R12_TEV + STP_TEV + SCI + MCI, method="class", data=mc.dat, control=rpart.control(minbucket=buck, xval=20))
	} else {
		tree <- rpart(class.dat ~ R12_TEV + STP_TEV + SCI + MCI, method="class", data=mc.dat, control=rpart.control(xval=20))
	}
	frum <- tree$frame[,c(1:5)]
	node.vec <- rownames(frum)[frum$var=="<leaf>"]
	node.data <- frum$yval[frum$var=="<leaf>"]
	attach(mc.dat, warn.conflicts=F)
	rule.vec <- list()
	for (i in 1:length(node.vec)) {
	
		rty <- path.rpart(tree, node=node.vec[i], print.it=FALSE)
		rule.vec[[i]] <- rty[[1]]
	
	}
	
	node.list <- array(NA, dim=c(dim(est.data)[1]))
	for(i in 1:length(node.vec)) {
		rule.list <- "P"
		tdata <- mc.dat
		for (j in 2:length(rule.vec[[i]])){
			attach(tdata, warn.conflicts=F)
			tdata <- tdata[eval(parse(text=rule.vec[[i]][j])),]
		}
		rownumba <- which(as.character(est.data$est_name) %in% as.character(tdata$est_name))
		node.list[rownumba] <- i
	}
	outdat <- data.frame(class.dat, node.list)
	
	break.vec <- c(1:(length(unique(class.dat))+1))-0.5
	
	# nclass should be the number of terminal nodes
	nclass <- length(node.vec)
	nrows.plot <- ceiling(nclass/4)
	plot.mat <- matrix(c(rep(c(1,1,1,1), 7-nrows.plot),2:(nclass+1), rep(0,(28-(4*(7-nrows.plot))-length(2:(nclass+1))))),nrow=4,ncol=7, byrow=F)
	nf <- layout(plot.mat,heights=c(1,1))
	nf
	
	plot(tree)
	text(tree)
	
	for (i in 1:length(node.vec)) {
		plot.tit <- paste("Node ", i, sep="")
		hist(outdat$class.dat[outdat$node.list==i], breaks=break.vec, main=plot.tit, ylab="", xlab="Group", col=c(1:length(unique(class.dat))))	
	}
	
	tree.class0 <- 0
	bin.test <- 0
	for(i in 1:length(node.list)){
		tree.class0[i] <- node.data[node.list[i]]
		bin.test[i] <- if(tree.class0[i]==class.dat[i])1 else {0}
	}
	
	tm1 <- table.meth(class.dat, tree.class0)
	ind.error <- 0
	for(i in 1:length(unique(class.dat))) {
		ind.error[i] <- 100*diag(tm1)[i]/length(class.dat[class.dat==i])
	}
	
	perc.diff <- 100*mean(bin.test)
	Group <- c(as.character(unique(class.dat)),"Total")
	Error <- c(ind.error, perc.diff)
	Number <- c(len.vector, sum(len.vector))
	ret.dat <- data.frame(Group, Error, Number)
	return(ret.dat)
}

Tree.Construct(Mclust.class)

# Euclidean distance measure – ward grouping, unscaled variables
# 5 groups
dist.est.data.a1 <- dist(as.matrix(est.data.a1[,-1]), method="euclidean")
hc.est.data.a1 <- hclust(dist.est.data.a1, method="ward")
hc.cut.k5 <- cutree(hc.est.data.a1, k=5)
plotfunc(colscheme=hc.cut.k5, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"), colo=c(1,2,3,4,5))

Tree.Construct(hc.cut.k5, type="restricted")
Tree.Construct(hc.cut.k5, type="non")

# Euclidean distance measure – ward grouping, scaled variables
# 4 groups
dist.est.data.a2 <- dist(as.matrix(est.data.a2[,-1]), method="euclidean")
hc.est.data.a2 <- hclust(dist.est.data.a2, method="ward")
hc.cut.k4 <- cutree(hc.est.data.a2, k=4)
plotfunc(colscheme=hc.cut.k4, leglab=c("Group 1", "Group 2", "Group 3", "Group 4"), colo=c(1,2,3,4))

Tree.Construct(hc.cut.k4, type="restricted")
Tree.Construct(hc.cut.k4, type="non")

# Random Forests distance - ward grouping
# 7 groups

rf.data.a1 <- randomForest(est.data.a1[, -1], type="unsupervised", proximity=TRUE, ntree=50000)
dist.rf.a1 <- as.dist(sqrt(1-rf.data.a1$proximity))
hc.rf.a1 <- hclust(dist.rf.a1, method="ward")
rf.cut.k7 <- cutree(hc.rf.a1, k=7)

plotfunc(colscheme=rf.cut.k7, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7"), colo=c(1,2,3,4,5,6,7))

Tree.Construct(rf.cut.k7, type="restricted")
Tree.Construct(rf.cut.k7, type="non")

# Model based 7 clusters - scaled

Mclust.data.a2 <- Mclust(est.data.a2[,-1], G=2:7)
Mclust.class <- Mclust.data.a2$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7"), colo=c(1,2,3,4,5,6,7))

Tree.Construct(Mclust.class, type="restricted")
Tree.Construct(Mclust.class, type="non")



# Random Forests with average clustering method
# determine the number of groups that are supported

rf.data.a1 <- randomForest(est.data.a1[, -1], type="unsupervised", proximity=TRUE, ntree=10000)
rf.prox <- rf.data.a1$proximity
rf.prox[415,436] <- 0.9999999
rf.prox[436,415] <- 0.9999999
dist.rf.a1 <- as.dist(sqrt(1-rf.prox))

mds.trans.rf <- isoMDS(dist.rf.a1, k=2)
mds.trans.rf.scores <- mds.trans.rf$points

hc.rf.a1 <- hclust(dist.rf.a1, method="average")
rf.clust.ave <- cutree(hc.rf.a1, k=7)
table.meth(clust.dat=rf.clust.ave, act.dat=as.character(est.data$est.class))
Tree.Construct(rf.clust.ave, type="restricted")

hc.rf.a1 <- hclust(dist.rf.a1, method="ward")
rf.clust.ward <- cutree(hc.rf.a1, k=7)
table.meth(clust.dat=rf.clust.ward, act.dat=as.character(est.data$est.class))
Tree.Construct(rf.clust.ward, type="restricted")

hc.rf.a1 <- hclust(dist.rf.a1, method="single")
rf.clust.single <- cutree(hc.rf.a1, k=7)
table.meth(clust.dat=rf.clust.single, act.dat=as.character(est.data$est.class))
Tree.Construct(rf.clust.single, type="restricted")

hc.rf.a1 <- hclust(dist.rf.a1, method="complete")
rf.clust.complete <- cutree(hc.rf.a1, k=7)
table.meth(clust.dat=rf.clust.complete, act.dat=as.character(est.data$est.class))
Tree.Construct(rf.clust.complete, type="restricted")

pc1 <- princomp(est.data.a2[,-1])
pc1.scores <- pc1$scores

par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(Comp.1 ~ Comp.2, data=pc1.scores, ylab="PCA1", xlab="PCA2", pch=16, col=rf.clust.ward, main="Ward")
plot(Comp.1 ~ Comp.2, data=pc1.scores, ylab="PCA1", xlab="PCA2", pch=16, col=Mclust.class, main="Single")
plot(Comp.1 ~ Comp.2, data=pc1.scores, ylab="PCA1", xlab="PCA2", pch=16, col=rf.clust.complete, main="Complete")
plot(Comp.1 ~ Comp.2, data=pc1.scores, ylab="PCA1", xlab="PCA2", pch=16, col=rf.clust.ave, main="Average")

par(mfrow=c(2,2), mar=c(4,4,2,1))
plot(mds.trans.rf.scores[,1] ~ mds.trans.rf.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=rf.clust.ward, cex=1, main="Ward")
plot(mds.trans.rf.scores[,1] ~ mds.trans.rf.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=Mclust.class, cex=1, main="Single")
plot(mds.trans.rf.scores[,1] ~ mds.trans.rf.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=rf.clust.complete, cex=1, main="Complete")
plot(mds.trans.rf.scores[,1] ~ mds.trans.rf.scores[,2], ylab="Axis 1", xlab="Axis 2", pch=16, col=rf.clust.ave, cex=1, main="Average")


table.meth(rf.clust.ward, rf.clust.complete)

table.meth(rf.clust.ward, rf.clust.ave)

table.meth(rf.clust.complete, rf.clust.ave)


