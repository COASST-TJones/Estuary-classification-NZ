setwd("C:/Users/Timothy/Desktop/NZ estuarine classification")
library(vegan)
library(mclust)
library(clValid)
library(randomForest)
library(MASS)
library(rpart)
est.data <- read.csv(file="EEC_data.csv")
est.data$est.class[est.data$est.class=="C"] <- "B"

names(est.data)

est.data$R12_TEV  <- est.data$R12/est.data$TEV
est.data$STP_TEV  <- est.data$STP/est.data$TEV

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

# Examine variable distributions to determine transformations

par(mfrow=c(2,2))
hist(est.data$CLA,100, xlab="CLA", main="No Transform")
hist(sqrt(est.data$CLA),100, xlab="sqrt(CLA)", main="Square-root Transform")
hist(sqrt(sqrt(est.data$CLA)),100, xlab="(CLA)^(1/4)", main="Fourth-root Transform")
hist(log(est.data$CLA+1),100, xlab="log(CLA+1)", main="log(X+1) Transform")

plot(ecdf(est.data$CLA), xlab="CLA", main="No Transform", ylab="CDF(CLA)")
plot(ecdf(sqrt(est.data$CLA)), xlab="sqrt(CLA)", main="Square-root Transform", ylab="CDF(sqrt(CLA))")
plot(ecdf(sqrt(sqrt(est.data$CLA))), xlab="(CLA)^(1/4)", main="Fourth-root Transform",ylab="CDF((CLA)^(1/4))" )
plot(ecdf(log(est.data$CLA+1)), xlab="log(CLA+1)", main="log (X+1) Transform", ylab="CDF(log(CLA+1))" )

# prepare dataset

names(est.data)

set2.dat <- est.data[,c(1,4:12)]

set2.dat$STP <- sqrt(sqrt(set2.dat$STP))

set2.dat$R12 <- sqrt(sqrt(set2.dat$R12))

set2.dat$mean_depth <- sqrt(sqrt(set2.dat$mean_depth))

set2.dat$p_IA[set2.dat$p_IA==0] <- 0.05
set2.dat$p_IA[set2.dat$p_IA==100] <- 99.95
set2.dat$p_IA <- log(set2.dat$p_IA/(100-set2.dat$p_IA))

set2.dat$CLA <- log(set2.dat$CLA + 1)

set2.dat$TEV <- log(set2.dat$TEV)

set2.dat$EWA <- log(set2.dat$EWA)

set2.scale <- set2.dat

set2.scale$STP <- scale(set2.scale$STP)
set2.scale$R12 <- scale(set2.scale$R12)
set2.scale$mean_depth <- scale(set2.scale$mean_depth)
set2.scale$p_IA <- scale(set2.scale$p_IA)
set2.scale$CLA <- scale(set2.scale$CLA)
set2.scale$SCI <- scale(set2.scale$SCI)
set2.scale$MCI <- scale(set2.scale$MCI)
set2.scale$TEV <- scale(set2.scale$TEV)
set2.scale$EWA <- scale(set2.scale$EWA)

head(set2.scale)

# remove CLA and EWA 
# visualising the data

pc.set2 <- princomp(set2.scale[,-1])
pc.set2

plot.PCA <- function (PCA, colscheme) {
	
	advect <- function(loading,colx, coly, cirsize,cut.off){
	
		factnames <- rownames(loading)
		for(i in 1:length(factnames)){
			ynum <- loading[i,coly]
			xnum <- loads[i,colx]
			vect.len <- sqrt((ynum^2) + (xnum^2))
			if(vect.len > cut.off) {
			
				segments(x0=0, y0=0, x1=cirsize*xnum, y1=cirsize*ynum, lwd=2, col=2, lty=1)
				text(x=cirsize*xnum, y=cirsize*ynum, pos=3, labels=factnames[i], cex=1.4, col="grey50")
			
			}
		}
	
	}
	scores <- PCA$scores
	loads <- PCA$loadings
	par(mfrow=c(2,3), mar=c(4,4,1,1))
	
	plot(Comp.1 ~ Comp.2, data=scores, col=colscheme, pch=16, ylab="PCA1", xlab="PCA2")
	advect(loading=loads, colx=2, coly=1,cirsize=6,cut.off=0.5)
	plot(Comp.1 ~ Comp.3, data=scores, col=colscheme, pch=16, ylab="PCA1", xlab="PCA3")
	advect(loading=loads, colx=3, coly=1,cirsize=2.5,cut.off=0.5)
	plot(Comp.1 ~ Comp.4, data=scores, col=colscheme, pch=16, ylab="PCA1", xlab="PCA4")
	advect(loading=loads, colx=4, coly=1,cirsize=2.5,cut.off=0.5)

	plot(Comp.2 ~ Comp.3, data=scores, col=colscheme, pch=16, ylab="PCA2", xlab="PCA3")
	advect(loading=loads, colx=3, coly=2,cirsize=2.5,cut.off=0.5)
	plot(Comp.2 ~ Comp.4, data=scores, col=colscheme, pch=16, ylab="PCA2", xlab="PCA4")
	advect(loading=loads, colx=4, coly=2,cirsize=2.5,cut.off=0.5)
	
	plot(Comp.3 ~ Comp.4, data=scores, col=colscheme, pch=16, ylab="PCA3", xlab="PCA4")
	advect(loading=loads, colx=4, coly=3,cirsize=2.5,cut.off=0.5)
	
		
}


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

plotdiff <- function (group.labs) {

	t1 <- tapply(set2.scale$STP,group.labs, FUN=mean)
	t2 <- tapply(set2.scale$R12, group.labs, FUN=mean)
	t3 <- tapply(set2.scale$SCI, group.labs, FUN=mean)
	t4 <- tapply(set2.scale$MCI, group.labs, FUN=mean)
	t5 <- tapply(set2.scale$mean_depth, group.labs, FUN=mean)
	t6 <- tapply(set2.scale$p_IA, group.labs, FUN=mean)
	t7 <- tapply(set2.scale$TEV, group.labs, FUN=mean)
	
	s1 <- tapply(set2.scale$STP,group.labs, FUN=sd)
	s2 <- tapply(set2.scale$R12, group.labs, FUN=sd)
	s3 <- tapply(set2.scale$SCI, group.labs, FUN=sd)
	s4 <- tapply(set2.scale$MCI, group.labs, FUN=sd)
	s5 <- tapply(set2.scale$mean_depth, group.labs, FUN=sd)
	s6 <- tapply(set2.scale$p_IA, group.labs, FUN=sd)
	s7 <- tapply(set2.scale$TEV, group.labs, FUN=sd)
	
	par(mfrow=c(2,4), mar=c(3,5,2,1))
	
	plot(t1, pch=16, ylab="STP", xlab="", ylim=c(min(t1-s1), max(t1+s1)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t1[i]-s1[i], x1=i, y1=t1[i]+s1[i], col=1, lwd=2)
	}
	plot(t2, pch=16, ylab="R12", xlab="", ylim=c(min(t2-s2), max(t2+s2)), cex=1.5)
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
	plot(t5, pch=16, ylab="Mean depth", xlab="", ylim=c(min(t5-s5), max(t5+s5)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t5[i]-s5[i], x1=i, y1=t5[i]+s5[i], col=1, lwd=2)
	}
	plot(t6, pch=16, ylab="% Intertidal Area", xlab="", ylim=c(min(t6-s6), max(t6+s6)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t6[i]-s6[i], x1=i, y1=t6[i]+s6[i], col=1, lwd=2)
	}
	plot(t7, pch=16, ylab="TEV", xlab="", ylim=c(min(t7-s7), max(t7+s7)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t7[i]-s7[i], x1=i, y1=t7[i]+s7[i], col=1, lwd=2)
	}
}



# identifying an adequate number of groups

valid.clust <- function (inp.data, ngroups=3, n.left=5, lim=1000, verb=F, test.type="DIST") {
	
	len <- dim(inp.data)[1]													# length of original data
	perms <- rbind(c(1:ngroups), allPerms(ngroups, control=how(maxperm=1e7)))
	nperms <- dim(perms)[1]
	
	perc.match <- array(NA, dim=c(lim))
	diff.struct <- array(NA, dim=c(lim))
	
	# perform whole dataset clustering
	dist.orig <- dist(as.matrix(inp.data[,-1]), method="euclidean")			# calculation of distance matrix
	hc.orig <- hclust(dist.orig, method="ward")								# cluster analysis
	cut.orig <- cutree(hc.orig, k=ngroups)									# Find original grouping labels
	#print(cut.orig)
	if(verb==TRUE) {														# The verb==TRUE option enables the reporting of detailed information on individual misclassifications
		est.levels <- levels(as.factor(inp.data$no))						# defines a list of estuary labels
		err.frame <- array(0, dim=c(length(est.levels), 2*ngroups+1))		# creates an array for storing the number of individual classifications
		err.frame[,1] <- cut.orig											# place full dataset classifications into the first column
	}
	
	# **
	# Now loop to examine the cross-validated error
	for ( e in 1:lim) {
		# Split data into training and testing dataset
		id.left.out <- array(NA, dim=c(n.left))
		j.test <- 1
		while (j.test <= n.left) {
			ty <- floor(runif(1,min=1,max=len+0.9999999))# identify left out testing data
			if(ty %in% id.left.out==FALSE) {
				id.left.out[j.test] <- ty
				j.test <- j.test + 1
			}
		}
		left.out <- inp.data[id.left.out,]										# test dataset
		train <- inp.data[-id.left.out,]										# Training dataset
		# **
		#print(left.out)
		# Build tree based on training dataset
		dist.train <- dist(as.matrix(train[,-1]), method="euclidean")			# calculation of distance matrix
		hc.train <- hclust(dist.train, method="ward")							# Perform cluster analyses
		cut.train <- cutree(hc.train, k=ngroups)								# make cut at k=ngroups
		# **
		#print(data.frame(train,cut.train))
	
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
		#print(confu)
		
		diag.test <- array(NA, dim=c(nperms))									# creation of vector for storing test statistics
		
		for(i in 1:nperms) {													# loop through all possible arrangements of columns corresponding to re-labelling the training data cluster id's
			test.mat <- array(NA, dim=c(ngroups,ngroups))						# the best matching arrangement is determined as the one that achieves the highest matrix diagonal
			for ( j in 1:ngroups){												# Loop over groups to create the test matrix
				test.mat[,j]<- confu[,perms[i,j]]
			}
			diag.test[i] <- sum(diag(test.mat))									# calculate the diagonal test statistic
		}
		new.labs <- perms[which(diag.test==max(diag.test))[1],]					# list of new labels
		#print(new.labs)
		cut.train.adj <- array(NA, dim=c(length(cut.train)))					# assign new labels to data
		
		for(i in 1:ngroups) {													# Loop over groups
			cut.train.adj[cut.train==new.labs[i]] <- i							# Create new labelling based on match between train and original labelling schemes
		}
		#print(data.frame(train, cut.train.adj))
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
			STP.centro <- tapply(train$STP, cut.train.adj, FUN=mean)
			R12.centro <- tapply(train$R12, cut.train.adj, FUN=mean)
			p_IA.centro <- tapply(train$p_IA, cut.train.adj, FUN=mean)
			mean_depth.centro <- tapply(train$mean_depth, cut.train.adj, FUN=mean)
			TEV.centro <- tapply(train$TEV, cut.train.adj, FUN=mean)
			
			test.label <- array(NA, dim=c(n.left))
			for(j in 1:n.left) {													# Loop over left out entities
				dist.check <- array(NA, dim=c(ngroups))
				for(i in 1:ngroups) {												# Loop over possible groups for this entity
					dist.check[i] <- (MCI.centro[i] - left.out$MCI[j])^2 + (SCI.centro[i] - left.out$SCI[j])^2 + (STP.centro[i] - left.out$STP[j])^2 + (R12.centro[i] - left.out$R12[j])^2 + (p_IA.centro[i] - left.out$p_IA[j])^2 + (mean_depth.centro[i] - left.out$mean_depth[j])^2 + (TEV.centro[i] - left.out$TEV[j])^2
				}

				test.label[j] <- which(dist.check==min(dist.check))					# Identify the group label as the one which minimises the sum of within group variances
			}
		}
		#print(test.label)
		
		
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
	
	par(mfrow=c(2,2))
	hist(perc.match, 50)
	hist(diff.struct, 50)
	plot(perc.match ~ diff.struct, pch=16)
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

fold.vec <- c(20,40)
group.vec <- c(2,3,4,5,6,7,8,9,10)
lim.vec <- c(300,300,300,300,300,300,300,100,100)

error.vec <- 0
diff.vec <- 0
fold.out <- 0
group.out <- 0

K <- 1
for ( i in 1:2) {
	for( j in 1:9) {
		vc <- valid.clust(inp.data=set2.scale, ngroups=group.vec[j],n.left=fold.vec[i], lim=lim.vec[j])
		error.vec[K] <- vc$perc
		diff.vec[K] <- vc$dif
		fold.out[K] <- fold.vec[i]
		group.out[K] <- group.vec[j]
		K <- K+1
	}
}

cv.dataframe <- data.frame(fold.out, group.out, error.vec, diff.vec)

par(mfrow=c(2,2), mar=c(5,4,1,2))
plot(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=1, ylim=c(0,100), ylab="% classification error", xlab="", type="b")
lines(error.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=2, type="b")
abline(h=80, lty=2, col=3)
abline(h=90, lty=2, col=3)
text(x=10, y=86.5, labels="80%")
text(x=10, y=96.5, labels="90%")

X <- c(2,3,4,5,6,7,8,9,10)
Y <- 100/X

lines(Y ~ X, col=1, lty=2, lwd=2)

plot(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=1, ylim=c(0,50), ylab="% change in structure", xlab="Number of Groups", type="b")
lines(diff.vec ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=2, type="b")

legend("topleft", col=c(1,2,1), lty=c(1,1,2), legend=c("20 left out","40 left out", "Expected by Chance"), bty="n")

# calculate cohens kappa

coh.kappa <- (cv.dataframe$error.vec - (100/cv.dataframe$group.out))/(100-(100/cv.dataframe$group.out))
cv.dataframe <- data.frame(cv.dataframe, coh.kappa)

plot(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==20,], pch=16, col=1, ylim=c(0,1), ylab="Cohens Kappa", xlab="", type="b")
lines(coh.kappa ~ group.out, data=cv.dataframe[cv.dataframe$fold.out==40,], pch=16, col=2, type="b")

abline(h=0.2, col=1, lty=2, lwd=1)
abline(h=0.4, col=2, lty=2, lwd=1)
abline(h=0.6, col=3, lty=2, lwd=1)
abline(h=0.8, col=4, lty=2, lwd=1)

legend("bottomleft", col=c(1,2), lty=c(1,1), legend=c("20 left out", "40 left out"), bty="n")
legend("bottomright", col=c(1,2,3,4), lty=c(2,2,2,2), legend=c("Slight", "Fair", "Substantial", "Good"), bty="n")

# 5 groups

dist.orig <- dist(as.matrix(set2.scale[,-1]), method="euclidean")			# calculation of distance matrix
hc.eucl.set2 <- hclust(dist.orig, method="ward")

k5.eucl.set2 <- cutree(hc.eucl.set2, k=5)
plot.PCA(PCA=pc.set2, colscheme=k5.eucl.set2)
plotdiff(k5.eucl.set2)
table.meth(k5.eucl.set2, act.dat=as.character(est.data$est.class))

est.data[k5.eucl.set2==1,1:2]

# 6 groups

k6.eucl.set2 <- cutree(hc.eucl.set2, k=6)
plot.PCA(PCA=pc.set2, colscheme=k6.eucl.set2)
plotdiff(k6.eucl.set2)
table.meth(k6.eucl.set2, act.dat=as.character(est.data$est.class))

est.data[k6.eucl.set2==1,1:2]

###################################################################################################################################################################
# Random Forests Clustering
###################################################################################################################################################################

set2.rf <- est.data[,c(1,6:12)]
rf.prox.set2 <- randomForest(set2.rf[, -1], type="unsupervised", proximity=TRUE, ntrees=500000)
dist.rf.set2 <- as.dist(sqrt(1-rf.prox.set2$proximity))
hc.rf.set2 <- hclust(dist.rf.set2, method="ward")

# Determine the number of groups

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
			rf.cut.tree <- randomForest(cut.tree ~ SCI + MCI + STP + R12 + p_IA + mean_depth + TEV, data=rf.dat.test, sampsize=221,mtry=2, ntrees=500)
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

rf.clust.ident(data.obj=set2.rf,tree.obj=hc.rf.set2, kclust=c(2:30), nrepeats=1000) 

# 5 groups

k5.rf.set2 <- cutree(hc.rf.set2, k=5)
plot.PCA(PCA=pc.set2, colscheme=k5.rf.set2)
plotdiff(k5.rf.set2)
table.meth(k5.rf.set2, act.dat=as.character(est.data$est.class))

est.data[k5.rf.set2==1,1:2]

# 7 groups

k7.rf.set2 <- cutree(hc.rf.set2, k=7)
plot.PCA(PCA=pc.set2, colscheme=k7.rf.set2)
plotdiff(k7.rf.set2)
table.meth(k7.rf.set2, act.dat=as.character(est.data$est.class))

est.data[k7.rf.set2==1,1:2]

# 9 groups

k9.rf.set2 <- cutree(hc.rf.set2, k=9)
plot.PCA(PCA=pc.set2, colscheme=k9.rf.set2)
plotdiff(k9.rf.set2)
table.meth(k9.rf.set2, act.dat=as.character(est.data$est.class))

est.data[k9.rf.set2==1,1:2]

###################################################################################################################################################################
# Model Based Clustering
###################################################################################################################################################################


MclustBIC.set2 <- mclustBIC(set2.scale[,-1], G=2:30)

# VEV 7 groups

Mclust.data.a1 <- Mclust(set2.scale[,-1], G=2:7)
Mclust.class <- Mclust.data.a1$classification
plot.PCA(PCA=pc.set2, colscheme=Mclust.class)
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))

est.data[Mclust.class==1,1:2]
est.data[Mclust.class==2,1:2]
est.data[Mclust.class==3,1:2]
est.data[Mclust.class==4,1:2]
est.data[Mclust.class==5,1:2]
est.data[Mclust.class==6,1:2]
est.data[Mclust.class==7,1:2]


# VEV 9 groups

Mclust.data.a1 <- Mclust(set2.scale[,-1], G=9:9)
Mclust.class <- Mclust.data.a1$classification
plot.PCA(PCA=pc.set2, colscheme=Mclust.class)
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))

est.data[Mclust.class==1,1:2]
est.data[Mclust.class==2,1:2]
est.data[Mclust.class==3,1:2]
est.data[Mclust.class==4,1:2]
est.data[Mclust.class==5,1:2]
est.data[Mclust.class==6,1:2]
est.data[Mclust.class==7,1:2]


####
# additional bits

setwd("C:/Users/Timothy/Desktop/NZ estuarine classification")
library(vegan)
library(mclust)
library(clValid)
library(randomForest)
library(MASS)
library(rpart)
est.data <- read.csv(file="EEC_data.csv")
est.data$est.class[est.data$est.class=="C"] <- "B"
names(est.data)

est.data$R12_TEV  <- est.data$R12/est.data$TEV
est.data$STP_TEV  <- est.data$STP/est.data$TEV

# est.data.a1 = dataset 1 transformed, unscaled

est.data.a1 <- est.data[,c(1,11,12,13,14)]
est.data.a1$R12_TEV <- sqrt(sqrt(est.data.a1$R12_TEV))
est.data.a1$STP_TEV <- sqrt(sqrt(est.data.a1$STP_TEV))

# est.data.a2 = dataset 1 transformed, scaled

est.data.a2 <- est.data.a1
est.data.a2$R12_TEV <- scale(est.data.a2$R12_TEV)
est.data.a2$STP_TEV <- scale(est.data.a2$STP_TEV)
est.data.a2$SCI <- scale(est.data.a2$SCI)
est.data.a2$MCI <- scale(est.data.a2$MCI)

# set1.scale = dataset 2, transformed and scaled

set1.dat <- est.data[,c(1,6,7,11,12,13,14,15)]
set1.dat$STP_TEV <- sqrt(sqrt(set1.dat$STP_TEV))
set1.dat$R12_TEV <- sqrt(sqrt(set1.dat$R12_TEV))
set1.dat$mean_depth <- sqrt(sqrt(set1.dat$mean_depth))
set1.dat$p_IA[set1.dat$p_IA==0] <- 0.05
set1.dat$p_IA[set1.dat$p_IA==100] <- 99.95
set1.dat$p_IA <- log(set1.dat$p_IA/(100-set1.dat$p_IA))
set1.dat$CLA_EWA <- log(set1.dat$CLA_EWA + 1)
set1.scale <- set1.dat
set1.scale$STP_TEV <- scale(set1.scale$STP_TEV)
set1.scale$R12_TEV <- scale(set1.scale$R12_TEV)
set1.scale$mean_depth <- scale(set1.scale$mean_depth)
set1.scale$p_IA <- scale(set1.scale$p_IA)
set1.scale$CLA_EWA <- scale(set1.scale$CLA_EWA)
set1.scale$SCI <- scale(set1.scale$SCI)
set1.scale$MCI <- scale(set1.scale$MCI)

# set2.scale = dataset 3, transformed and scaled

set2.dat <- est.data[,c(1,6:12)]
set2.dat$STP <- sqrt(sqrt(set2.dat$STP))
set2.dat$R12 <- sqrt(sqrt(set2.dat$R12))
set2.dat$mean_depth <- sqrt(sqrt(set2.dat$mean_depth))
set2.dat$p_IA[set2.dat$p_IA==0] <- 0.05
set2.dat$p_IA[set2.dat$p_IA==100] <- 99.95
set2.dat$p_IA <- log(set2.dat$p_IA/(100-set2.dat$p_IA))
set2.dat$TEV <- log(set2.dat$TEV)
set2.scale <- set2.dat
set2.scale$STP <- scale(set2.scale$STP)
set2.scale$R12 <- scale(set2.scale$R12)
set2.scale$mean_depth <- scale(set2.scale$mean_depth)
set2.scale$p_IA <- scale(set2.scale$p_IA)
set2.scale$SCI <- scale(set2.scale$SCI)
set2.scale$MCI <- scale(set2.scale$MCI)
set2.scale$TEV <- scale(set2.scale$TEV)

###############
## Functions ##
###############

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

Tree.Construct <- function(class.dat, type="restricted") {

	len.vector <- tapply(class.dat, class.dat, FUN=length)
	buck <- round(min(len.vector)*0.5)
	mc.dat <- data.frame(est.data, class.dat)
	if(type=="restricted") {
		tree <- rpart(class.dat ~ R12_TEV + STP_TEV + SCI + MCI + p_IA + mean_depth + CLA_EWA, method="class", data=mc.dat, control=rpart.control(minbucket=buck, xval=20))
	} else {
		tree <- rpart(class.dat ~ R12_TEV + STP_TEV + SCI + MCI + p_IA + mean_depth + CLA_EWA, method="class", data=mc.dat, control=rpart.control(xval=20))
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

RForestCode <- function (inp.data, ngroups, lim, set) {
	
	class.lab <- array(NA, dim=c(dim(inp.data)[1],lim))
	rf.obj <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntrees=100000)
	rf.dist <- as.dist(sqrt(1-rf.obj$proximity))
	hc.rf <- hclust(rf.dist, method="ward")
	class.lab[,1] <- cutree(hc.rf, k=ngroups)
	
	perms <- rbind(c(1:ngroups), allPerms(ngroups, control=how(maxperm=1e7)))
	nperms <- dim(perms)[1]
	
	for (e in 2:lim){
		
		rf.obj <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntrees=100000)
		rf.dist <- as.dist(sqrt(1-rf.obj$proximity))
		hc.rf <- hclust(rf.dist, method="ward")
		
		test.lab <- cutree(hc.rf, k=ngroups)
		orig.lab <- class.lab[,1]
		
		confu <- array(NA, dim=c(ngroups, ngroups))								# Definition of confusion matrix between labellings
		
		for (i in 1:ngroups) {													# Loop over groups
			orig.names <- inp.data[orig.lab==i,1]								# Find list of estuary id's for the original labelling of group i
			for (j in 1:ngroups) {						
				new.names <- inp.data[test.lab==j,1]								# Find list of estuary id's for the training data labelling of group i
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
		
		for(i in 1:ngroups) {													# Loop over groups
			class.lab[test.lab==new.labs[i],e] <- i							# Create new labelling based on match between train and original labelling schemes
		}
	}
	
	class.mat <- array(NA, dim=c(dim(inp.data)[1],ngroups))
	Fin.class <- array(NA, dim=c(dim(inp.data)[1]))
	
	for (i in 1:dim(inp.data)[1]) {
		for(j in 1:ngroups){
			tlist <- class.lab[i,]
			class.mat[i,j] <-  100*length(tlist[tlist==j])/lim
		}
		Fin.class[i] <- which(class.mat[i,]==max(class.mat[i,]))
	}
	
	filename <- paste("RF_class_",set,"_", ngroups,".csv", sep="")
	outdat <- data.frame(est.data[,1:2], Fin.class, class.mat)
	namecol <- "P"
	for(i in 1:ngroups){
		namecol[i] <- paste("% Class ", i, sep="")
	}
	names(outdat) <- c("Number", "Name", "Modal Class", namecol)
	write.csv(outdat, file=filename)
	return(Fin.class)
}

ForestGroup <- function (inp.data, ngroups, lim, set) {
	
	class.lab <- array(NA, dim=c(dim(inp.data)[1],lim))
	rf.obj <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntrees=100000)
	rf.dist <- as.dist(sqrt(1-rf.obj$proximity))
	hc.rf <- hclust(rf.dist, method="ward")
	class.lab[,1] <- cutree(hc.rf, k=ngroups)
	
	perms <- rbind(c(1:ngroups), allPerms(ngroups, control=how(maxperm=1e7)))
	nperms <- dim(perms)[1]
	
	for (e in 2:lim){
		
		rf.obj <- randomForest(inp.data[, -1], type="unsupervised", proximity=TRUE, ntrees=100000)
		rf.dist <- as.dist(sqrt(1-rf.obj$proximity))
		hc.rf <- hclust(rf.dist, method="ward")
		
		test.lab <- cutree(hc.rf, k=ngroups)
		orig.lab <- class.lab[,1]
		
		confu <- array(NA, dim=c(ngroups, ngroups))								# Definition of confusion matrix between labellings
		
		for (i in 1:ngroups) {													# Loop over groups
			orig.names <- inp.data[orig.lab==i,1]								# Find list of estuary id's for the original labelling of group i
			for (j in 1:ngroups) {						
				new.names <- inp.data[test.lab==j,1]								# Find list of estuary id's for the training data labelling of group i
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
		
		for(i in 1:ngroups) {													# Loop over groups
			class.lab[test.lab==new.labs[i],e] <- i							# Create new labelling based on match between train and original labelling schemes
		}
	}
	
	class.mat <- array(NA, dim=c(dim(inp.data)[1],ngroups))
	Fin.class <- array(NA, dim=c(dim(inp.data)[1]))
	max.class.mem <- array(NA, dim=c(dim(inp.data)[1]))
	
	for (i in 1:dim(inp.data)[1]) {
		for(j in 1:ngroups){
			tlist <- class.lab[i,]
			class.mat[i,j] <-  100*length(tlist[tlist==j])/lim
		}
		Fin.class[i] <- which(class.mat[i,]==max(class.mat[i,]))
		max.class.mem[i] <- max(class.mat[i,])
	}
	
	filename <- paste("RF_class_",set,"_", ngroups,".csv", sep="")
	outdat <- data.frame(est.data[,1:2], Fin.class, class.mat)
	namecol <- "P"
	for(i in 1:ngroups){
		namecol[i] <- paste("% Class ", i, sep="")
	}
	names(outdat) <- c("Number", "Name", "Modal Class", namecol)
	write.csv(outdat, file=filename)
	
	mean.class.mem <- mean(max.class.mem)
	lo.class.mem <- max.class.mem[rank(max.class.mem, ties.method="random")==round(0.05*443)] 
	hi.class.mem <- max.class.mem[rank(max.class.mem, ties.method="random")==round(0.95*443)] 
	
	class.mem.out <- c(mean.class.mem, lo.class.mem, hi.class.mem)
	return(class.mem.out)
}


# data set 1

groupvec <- c(2,3,4,5,6,7,8,9,10)
limvec <- c(500,500,500,500,500,400,300,200,100)
err.vec.rf <- array(NA, dim=c(9,3))
for(i in 1:9){
	err.vec.rf[i,] <- ForestGroup(inp.data=est.data[,c(1,11:14)], ngroups=groupvec[i], lim=limvec[i], set="S1")
}

# data set 2

err.vec.rf2 <- array(NA, dim=c(9,3))
for(i in 1:9){
	err.vec.rf2[i,] <- ForestGroup(inp.data=est.data[,c(1,6,7,11,12,13,14,15)], ngroups=groupvec[i], lim=limvec[i], set="S2")
}

# data set 3

err.vec.rf3 <- array(NA, dim=c(9,3))
for(i in 1:9){
	err.vec.rf3[i,] <- ForestGroup(inp.data=est.data[,c(1,6:12)], ngroups=groupvec[i], lim=limvec[i], set="S3")
}

set1.rf.dec <- data.frame(Groups, err.vec.rf)
names(set1.rf.dec) <- c("Groups", "mean.err", "low.err", "high.err")
set2.rf.dec <- data.frame(Groups, err.vec.rf2)
names(set2.rf.dec) <- c("Groups", "mean.err", "low.err", "high.err")
set3.rf.dec <- data.frame(Groups, err.vec.rf3)
names(set3.rf.dec) <- c("Groups", "mean.err", "low.err", "high.err")

#
par(mfrow=c(2,2), oma=c(2,2,0,0), mar=c(3,3,2,2))
plot(mean.err ~ Groups,data=set1.rf.dec, type="l", lwd=2, col=1, ylim=c(0,100), xlab="", ylab="", main="Dataset 1")
lines(low.err ~ Groups, data=set1.rf.dec, type="l", lwd=2, col=2,lty=2)
lines(high.err ~ Groups, data=set1.rf.dec, type="l", lwd=2, col=2, lty=2)
abline(h=80, lty=2, col="grey50")

plot(mean.err ~ Groups,data=set2.rf.dec, type="l", lwd=2, col=1, ylim=c(0,100), xlab="", ylab="", main="Dataset 2")
lines(low.err ~ Groups, data=set2.rf.dec, type="l", lwd=2, col=2,lty=2)
lines(high.err ~ Groups, data=set2.rf.dec, type="l", lwd=2, col=2, lty=2)
abline(h=80, lty=2, col="grey50")

plot(mean.err ~ Groups,data=set3.rf.dec, type="l", lwd=2, col=1, ylim=c(0,100), xlab="", ylab="", main="Dataset 3")
lines(low.err ~ Groups, data=set3.rf.dec, type="l", lwd=2, col=2,lty=2)
lines(high.err ~ Groups, data=set3.rf.dec, type="l", lwd=2, col=2, lty=2)
abline(h=80, lty=2, col="grey50")

mtext(side=2, outer=T, text="Mean maximum group membership (%)")
mtext(side=1, outer=T, text="Number of Groups")

#########
#########

# list of cluster scenarios to evaluate

# Set 1 (NS) Eucl Ward k=5 		s1.NS.EW.k5
dist.est.data.a1 <- dist(as.matrix(est.data.a1[,-1]), method="euclidean")
hc.est.data.a1 <- hclust(dist.est.data.a1, method="ward")
s1.NS.EW.k5 <- cutree(hc.est.data.a1, k=5)

# Set 1 (NS) Eucl Ward k=6	 	s1.NS.EW.k6
s1.NS.EW.k6 <- cutree(hc.est.data.a1, k=6)

# Set 1 (S) Eucl Ward k=4	 	s1.S.EW.k4
dist.est.data.a2 <- dist(as.matrix(est.data.a2[,-1]), method="euclidean")
hc.est.data.a2 <- hclust(dist.est.data.a2, method="ward")
s1.S.EW.k4 <- cutree(hc.est.data.a2, k=4)

# Set 1 (S) Eucl Ward k=5	 	s1.S.EW.k5
s1.S.EW.k5 <- cutree(hc.est.data.a2, k=5)

# Set 1 RF Ward k=6 			s1.RFW.k6
s1.RFW.k6 <- RForestCode(inp.data=est.data[,c(1,11:14)], ngroups=6, lim=500, set="S1")

# Set 1 RF Ward k=7 			s1.RFW.k7
s1.RFW.k7 <- RForestCode(inp.data=est.data[,c(1,11:14)], ngroups=7, lim=500, set="S1")

# Set 1 (NS) Model EEV k=12  	s1.NS.EEV.k12
Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:12)
s1.NS.EEV.k12 <- Mclust.data.a1$classification

# Set 1 (NS) Model VEV k=9		s1.NS.VEV.k9
Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:9)
s1.NS.VEV.k9 <- Mclust.data.a1$classification

# Set 1 (NS) Model VVI k=6		s1.NS.VVI.k6
Mclust.data.a1 <- Mclust(est.data.a1[,-1], G=2:6)
s1.NS.VVI.k6 <- Mclust.data.a1$classification

# Set 1 (S) Model VEV k=6		s1.S.VEV.k6
Mclust.data.a1 <- Mclust(est.data.a2[,-1], G=2:6)
s1.S.VEV.k6 <- Mclust.data.a1$classification

# Set 1 (S) Model VEV k=7		s1.S.VEV.k7
Mclust.data.a1 <- Mclust(est.data.a2[,-1], G=7:7)
s1.S.VEV.k7 <- Mclust.data.a1$classification

# Set 2 Eucl Ward k=8			s2.EW.k8
dist.est.data.a1 <- dist(as.matrix(set1.scale[,-1]), method="euclidean")
hc.est.data.a1 <- hclust(dist.est.data.a1, method="ward")
s2.EW.k8 <- cutree(hc.est.data.a1, k=8)

# Set 2 Eucl Ward k=9			s2.EW.k9
s2.EW.k9 <- cutree(hc.est.data.a1, k=9)

# Set 2 RF Ward k=5				s2.RFW.k5
s2.RFW.k5 <- RForestCode(inp.data=est.data[,c(1,6,7,11,12,13,14,15)], ngroups=5, lim=500, set="S2")

# Set 2 RF Ward k=7				s2.RFW.k7
s2.RFW.k7 <- RForestCode(inp.data=est.data[,c(1,6,7,11,12,13,14,15)], ngroups=7, lim=500, set="S2")

# Set 2 RF Ward k=9				s2.RFW.k9
s2.RFW.k9 <- RForestCode(inp.data=est.data[,c(1,6,7,11,12,13,14,15)], ngroups=9, lim=300, set="S2")

# Set 2 Model VEV k=9			s2.VEV.k9
Mclust.data.a1 <- Mclust(set1.scale[,-1], G=2:9)
s2.VEV.k9 <- Mclust.data.a1$classification

# Set 2 Model VEV k=10			s2.VEV.k10
Mclust.data.a1 <- Mclust(set1.scale[,-1], G=2:10)
s2.VEV.k10 <- Mclust.data.a1$classification

# Set 3 Eucl Ward k=5			s3.EW.k5
dist.est.data.a1 <- dist(as.matrix(set2.scale[,-1]), method="euclidean")
hc.est.data.a1 <- hclust(dist.est.data.a1, method="ward")
s3.EW.k5 <- cutree(hc.est.data.a1, k=5)

# Set 3 Eucl Ward k=6			s3.EW.k6
s3.EW.k6 <- cutree(hc.est.data.a1, k=6)

# Set 3 RF Ward k=5				s3.RFW.k5
s3.RFW.k5 <- RForestCode(inp.data=est.data[,c(1,6:12)], ngroups=5, lim=500, set="S3")

# Set 3 RF Ward k=7				s3.RFW.k7
s3.RFW.k7 <- RForestCode(inp.data=est.data[,c(1,6:12)], ngroups=7, lim=500, set="S3")

# Set 3 RF Ward k=9				s3.RFW.k9
s3.RFW.k9 <- RForestCode(inp.data=est.data[,c(1,6:12)], ngroups=9, lim=300, set="S3")

# Set 3 Model VEV k=7			s3.VEV.k7
Mclust.data.a1 <- Mclust(set2.scale[,-1], G=2:7)
s3.VEV.k7 <- Mclust.data.a1$classification

# Set 3 Model VEV k=9			s3.VEV.k9
Mclust.data.a1 <- Mclust(set2.scale[,-1], G=9:9)
s3.VEV.k9 <- Mclust.data.a1$classification

# Evaluate rule sets

Tree.Construct <- function(inp.dat=est.data, class.dat, type="restricted", contl=1) {

	mc.dat <- data.frame(inp.dat, class.dat)
	
	len.vector <- tapply(class.dat, class.dat, FUN=length)
	buck <- round(min(len.vector)*contl)
		
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
	
	node.list <- array(NA, dim=c(dim(mc.dat)[1]))
	for(i in 1:length(node.vec)) {
		rule.list <- "P"
		tdata <- mc.dat
		for (j in 2:length(rule.vec[[i]])){
			attach(tdata, warn.conflicts=F)
			tdata <- tdata[eval(parse(text=rule.vec[[i]][j])),]
		}
		rownumba <- which(as.character(mc.dat$est_name) %in% as.character(tdata$est_name))
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

rangefin <- function(class.scheme, set=1, ngroups) {
	if(set==1) {
		range.vec <- array(NA, dim=c(4,2*ngroups))
		Tile <- array("P", dim=c(2*ngroups))
		for(i in 1:ngroups){
			range.vec[1,((2*i)-1)] <- min(est.data$R12_TEV[class.scheme==i])
			range.vec[1,(2*i)] <- max(est.data$R12_TEV[class.scheme==i])
			range.vec[2,((2*i)-1)] <- min(est.data$STP_TEV[class.scheme==i])
			range.vec[2,(2*i)] <- max(est.data$STP_TEV[class.scheme==i])
			range.vec[3,((2*i)-1)] <- min(est.data$SCI[class.scheme==i])
			range.vec[3,(2*i)] <- max(est.data$SCI[class.scheme==i])
			range.vec[4,((2*i)-1)] <- min(est.data$MCI[class.scheme==i])
			range.vec[4,(2*i)] <- max(est.data$MCI[class.scheme==i])
			
			Tile[(2*i)-1] <- paste("Min ", i, sep="")
			Tile[(2*i)] <- paste("Max ", i, sep="")
		}
		Var.Name <- c("R12/TEV", "STP/TEV", "SCI", "MCI")
		outdata <- data.frame(Var.Name, range.vec)
		names(outdata) <- c("Var.Name", Tile)
	} else {
		if(set==2) {
			range.vec <- array(NA, dim=c(7,2*ngroups))
			Tile <- array("P", dim=c(2*ngroups))
			for(i in 1:ngroups){
				range.vec[1,((2*i)-1)] <- min(est.data$R12_TEV[class.scheme==i])
				range.vec[1,(2*i)] <- max(est.data$R12_TEV[class.scheme==i])
				range.vec[2,((2*i)-1)] <- min(est.data$STP_TEV[class.scheme==i])
				range.vec[2,(2*i)] <- max(est.data$STP_TEV[class.scheme==i])
				range.vec[3,((2*i)-1)] <- min(est.data$SCI[class.scheme==i])
				range.vec[3,(2*i)] <- max(est.data$SCI[class.scheme==i])
				range.vec[4,((2*i)-1)] <- min(est.data$MCI[class.scheme==i])
				range.vec[4,(2*i)] <- max(est.data$MCI[class.scheme==i])
				range.vec[5,((2*i)-1)] <- min(est.data$mean_depth[class.scheme==i])
				range.vec[5,(2*i)] <- max(est.data$mean_depth[class.scheme==i])
				range.vec[6,((2*i)-1)] <- min(est.data$p_IA[class.scheme==i])
				range.vec[6,(2*i)] <- max(est.data$p_IA[class.scheme==i])
				range.vec[7,((2*i)-1)] <- min(est.data$CLA_EWA[class.scheme==i])
				range.vec[7,(2*i)] <- max(est.data$CLA_EWA[class.scheme==i])
				Tile[(2*i)-1] <- paste("Min ", i, sep="")
				Tile[(2*i)] <- paste("Max ", i, sep="")
			}
			
			
			Var.Name <- c("R12/TEV", "STP/TEV", "SCI", "MCI", "Mean Depth", "% IA", "CLA/EWA")
			outdata <- data.frame(Var.Name, range.vec)
			names(outdata) <- c("Var.Name", Tile)
		} else {
			range.vec <- array(NA, dim=c(7,2*ngroups))
			Tile <- array("P", dim=c(2*ngroups))
			for(i in 1:ngroups){
				range.vec[1,((2*i)-1)] <- min(est.data$R12[class.scheme==i])
				range.vec[1,(2*i)] <- max(est.data$R12[class.scheme==i])
				range.vec[2,((2*i)-1)] <- min(est.data$STP[class.scheme==i])
				range.vec[2,(2*i)] <- max(est.data$STP[class.scheme==i])
				range.vec[3,((2*i)-1)] <- min(est.data$SCI[class.scheme==i])
				range.vec[3,(2*i)] <- max(est.data$SCI[class.scheme==i])
				range.vec[4,((2*i)-1)] <- min(est.data$MCI[class.scheme==i])
				range.vec[4,(2*i)] <- max(est.data$MCI[class.scheme==i])
				range.vec[5,((2*i)-1)] <- min(est.data$mean_depth[class.scheme==i])
				range.vec[5,(2*i)] <- max(est.data$mean_depth[class.scheme==i])
				range.vec[6,((2*i)-1)] <- min(est.data$p_IA[class.scheme==i])
				range.vec[6,(2*i)] <- max(est.data$p_IA[class.scheme==i])
				range.vec[7,((2*i)-1)] <- min(est.data$TEV[class.scheme==i])
				range.vec[7,(2*i)] <- max(est.data$TEV[class.scheme==i])
				Tile[(2*i)-1] <- paste("Min ", i, sep="")
				Tile[(2*i)] <- paste("Max ", i, sep="")
			}
			Var.Name <- c("R12", "STP", "SCI", "MCI", "Mean Depth", "% IA", "TEV")
			outdata <- data.frame(Var.Name, range.vec)
			names(outdata) <- c("Var.Name", Tile)
		}
	}
	return(outdata)
}

# s1.NS.EW.k5

plotdiff(s1.NS.EW.k5)
rangefin(class.scheme=s1.NS.EW.k5, set=1, ngroups=5)
Tree.Construct(class.dat=s1.NS.EW.k5, type="restricted", contl=1)

# s1.NS.EW.k6

plotdiff(s1.NS.EW.k6)
rangefin(class.scheme=s1.NS.EW.k6, set=1, ngroups=6)
Tree.Construct(class.dat=s1.NS.EW.k6, type="restricted", contl=1)

# s1.S.EW.k4

plotdiff(s1.S.EW.k4)
rangefin(class.scheme=s1.S.EW.k4, set=1, ngroups=4)

ne.class <- s1.S.EW.k4[s1.S.EW.k4 != 2]
ne.data <- est.data[s1.S.EW.k4 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

# s1.S.EW.k5

plotdiff(s1.S.EW.k5)
rangefin(class.scheme=s1.S.EW.k5, set=1, ngroups=5)

ne.class <- s1.S.EW.k5[s1.S.EW.k5 != 2]
ne.data <- est.data[s1.S.EW.k5 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

# Set 1 RF Ward k=6 			s1.RFW.k6
fgh <- read.csv(file="Class_s1_RFW_k6.csv")
s1.RFW.k6 <- fgh$class
plotdiff(s1.RFW.k6)
rangefin(class.scheme=s1.RFW.k6, set=1, ngroups=6)
Tree.Construct(class.dat=s1.RFW.k6, type="restricted", contl=1)

plot(sqrt(sqrt(est.data$R12_TEV)) ~ s1.RFW.k6, pch=16, col=rgb(0,0,0,0.4))

# Set 1 RF Ward k=7 			s1.RFW.k7
plotdiff(s1.RFW.k7)
rangefin(class.scheme=s1.RFW.k7, set=1, ngroups=7)
Tree.Construct(class.dat=s1.RFW.k7, type="restricted", contl=1)

plot(sqrt(sqrt(est.data$R12_TEV)) ~ s1.RFW.k7, pch=16, col=rgb(0,0,0,0.4))


# Set 1 (NS) Model EEV k=12  	s1.NS.EEV.k12
plotdiff(s1.NS.EEV.k12)
rangefin(class.scheme=s1.NS.EEV.k12, set=1, ngroups=12)
Tree.Construct(class.dat=s1.NS.EEV.k12, type="restricted", contl=1)

# Set 1 (NS) Model VEV k=9		s1.NS.VEV.k9

plotdiff(s1.NS.VEV.k9)
rangefin(class.scheme=s1.NS.VEV.k9, set=1, ngroups=9)

Tree.Construct(class.dat=s1.NS.VEV.k9, type="restricted", contl=1)

# Set 1 (NS) Model VVI k=6		s1.NS.VVI.k6
plotdiff(s1.NS.VVI.k6)
rangefin(class.scheme=s1.NS.VVI.k6, set=1, ngroups=6)

Tree.Construct(class.dat=s1.NS.VVI.k6, type="restricted", contl=1)

# Set 1 (S) Model VEV k=6		s1.S.VEV.k6

plotdiff(s1.S.VEV.k6)
rangefin(class.scheme=s1.S.VEV.k6, set=1, ngroups=6)

ne.class <- s1.S.VEV.k6[s1.S.VEV.k6 != 2]
ne.data <- est.data[s1.S.VEV.k6 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

Tree.Construct(class.dat=s1.S.VEV.k6, type="restricted", contl=1)

# Set 1 (S) Model VEV k=7		s1.S.VEV.k7

plotdiff(s1.S.VEV.k7)
rangefin(class.scheme=s1.S.VEV.k7, set=1, ngroups=7)

ne.class <- s1.S.VEV.k7[s1.S.VEV.k7 != 2]
ne.data <- est.data[s1.S.VEV.k7 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6
Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

Tree.Construct(class.dat=s1.S.VEV.k7, type="restricted", contl=1)

# Set 2 Eucl Ward k=8			s2.EW.k8
s2.EW.k8
plotdiff(s2.EW.k8)
rangefin(class.scheme=s2.EW.k8, set=2, ngroups=8)

ne.class <- s2.EW.k8[s2.EW.k8 != 2 & s2.EW.k8 != 8]
ne.data <- est.data[s2.EW.k8 != 2 & s2.EW.k8 != 8,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)
Tree.Construct(class.dat=s2.EW.k8, type="restricted", contl=1)

# Set 2 Eucl Ward k=9			s2.EW.k9
s2.EW.k9
plotdiff(s2.EW.k9)
rangefin(class.scheme=s2.EW.k9, set=2, ngroups=9)

ne.class <- s2.EW.k9[s2.EW.k9 != 2 & s2.EW.k9 != 9]
ne.data <- est.data[s2.EW.k9 != 2 & s2.EW.k9 != 9,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6
ne.class[ne.class==8] <- 7

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

# Set 2 RF Ward k=5				s2.RFW.k5
s2.RFW.k5
plotdiff(s2.RFW.k5)
rangefin(class.scheme=s2.RFW.k5, set=2, ngroups=5)

ne.class <- s2.RFW.k5[s2.RFW.k5 != 2]
ne.data <- est.data[s2.RFW.k5 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

# Set 2 RF Ward k=7				s2.RFW.k7
s2.RFW.k7
plotdiff(s2.RFW.k7)
rangefin(class.scheme=s2.RFW.k7, set=2, ngroups=7)

ne.class <- s2.RFW.k7[s2.RFW.k7 != 2]
ne.data <- est.data[s2.RFW.k7 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6
Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)


# Set 2 RF Ward k=9				s2.RFW.k9
s2.RFW.k9
plotdiff(s2.RFW.k9)
rangefin(class.scheme=s2.RFW.k9, set=2, ngroups=9)

ne.class <- s2.RFW.k9[s2.RFW.k9 != 2]
ne.data <- est.data[s2.RFW.k9 != 2,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6
ne.class[ne.class==8] <- 7
ne.class[ne.class==9] <- 8
Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)


# Set 2 Model VEV k=9			s2.VEV.k9
s2.VEV.k9
plotdiff(s2.VEV.k9)
rangefin(class.scheme=s2.VEV.k9, set=2, ngroups=9)

ne.class <- s2.VEV.k9[s2.VEV.k9 != 2 & s2.VEV.k9 != 9]
ne.data <- est.data[s2.VEV.k9 != 2 & s2.VEV.k9 != 9,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6
ne.class[ne.class==8] <- 7

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)


# Set 2 Model VEV k=10			s2.VEV.k10
s2.VEV.k10
plotdiff(s2.VEV.k10)
rangefin(class.scheme=s2.VEV.k10, set=2, ngroups=10)

ne.class <- s2.VEV.k10[s2.VEV.k10 != 2 & s2.VEV.k10 != 10]
ne.data <- est.data[s2.VEV.k10 != 2 & s2.VEV.k10 != 10,]
ne.class[ne.class==3] <- 2
ne.class[ne.class==4] <- 3
ne.class[ne.class==5] <- 4
ne.class[ne.class==6] <- 5
ne.class[ne.class==7] <- 6
ne.class[ne.class==8] <- 7
ne.class[ne.class==9] <- 8
Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

# Set 3 Eucl Ward k=5			s3.EW.k5
s3.EW.k5
plotdiff(s3.EW.k5)
rangefin(class.scheme=s3.EW.k5, set=3, ngroups=5)
Tree.Construct(inp.dat=est.data,class.dat=s3.EW.k5, type="restricted", contl=1)
# Set 3 Eucl Ward k=6			s3.EW.k6
s3.EW.k6
plotdiff(s3.EW.k6)
rangefin(class.scheme=s3.EW.k6, set=3, ngroups=6)

ne.class <- s3.EW.k6[s3.EW.k6 != 6]
ne.data <- est.data[s3.EW.k6 != 6,]

Tree.Construct(inp.dat=ne.data,class.dat=ne.class, type="restricted", contl=1)

# Set 3 RF Ward k=5				s3.RFW.k5

s3.RFW.k5
plotdiff(s3.RFW.k5)
rangefin(class.scheme=s3.RFW.k5, set=3, ngroups=5)
Tree.Construct(inp.dat=est.data,class.dat=s3.RFW.k5, type="restricted", contl=1)

# Set 3 RF Ward k=7				s3.RFW.k7
s3.RFW.k7
plotdiff(s3.RFW.k7)
rangefin(class.scheme=s3.RFW.k7, set=3, ngroups=7)
Tree.Construct(inp.dat=est.data,class.dat=s3.RFW.k7, type="restricted", contl=1)


# Set 3 RF Ward k=9				s3.RFW.k9
s3.RFW.k9
plotdiff(s3.RFW.k9)
rangefin(class.scheme=s3.RFW.k9, set=3, ngroups=9)
Tree.Construct(inp.dat=est.data,class.dat=s3.RFW.k9, type="restricted", contl=1)



# Set 3 Model VEV k=7			s3.VEV.k7
s3.VEV.k7
plotdiff(s3.VEV.k7)
rangefin(class.scheme=s3.VEV.k7, set=3, ngroups=7)
Tree.Construct(inp.dat=est.data,class.dat=s3.VEV.k7, type="restricted", contl=1)


# Set 3 Model VEV k=9			s3.VEV.k9
s3.VEV.k9
plotdiff(s3.VEV.k9)
rangefin(class.scheme=s3.VEV.k9, set=3, ngroups=9)
Tree.Construct(inp.dat=est.data,class.dat=s3.VEV.k9, type="restricted", contl=1)


