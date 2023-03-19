###########################################################################################################################################
# Section 2
###########################################################################################################################################

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
hist(est.data$CLA_EWA,100, xlab="CLA/EWA", main="No Transform")
hist(sqrt(est.data$CLA_EWA),100, xlab="sqrt(CLA/EWA)", main="Square-root Transform")
hist(sqrt(sqrt(est.data$CLA_EWA)),100, xlab="(CLA/EWA)^(1/4)", main="Fourth-root Transform")
hist(log(est.data$CLA_EWA+1),100, xlab="log(CLA/EWA + 1)", main="log(X+1) Transform")

plot(ecdf(est.data$CLA_EWA), xlab="CLA/EWA", main="No Transform", ylab="CDF(CLA/EWA)")
plot(ecdf(sqrt(est.data$CLA_EWA)), xlab="sqrt(CLA/EWA)", main="Square-root Transform", ylab="CDF(sqrt(CLA/EWA))")
plot(ecdf(sqrt(sqrt(est.data$CLA_EWA))), xlab="(CLA/EWA)^(1/4)", main="Fourth-root Transform",ylab="CDF((CLA/EWA)^(1/4))" )
plot(ecdf(log(est.data$CLA_EWA+1)), xlab="log(CLA/EWA + 1)", main="log(X+1) Transform", ylab="CDF(log(CLA/EWA+1))" )


# Set 1
# •	SCI
# •	MCI
# •	STP_TEV (use 4th root)
# •	R12_TEV (use 4th root)
# •	mean_depth (4th root)
# •	p_IA (logit)
# •	CLA_EWA (log(X+1))

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

head(set1.scale)

# plot factors against each other


plotfunc <- function(colscheme, leglab, colo) {

	par(mfrow=c(4,6), mar=c(3,3,0.5,0.5), xpd=TRUE)

	plot(R12_TEV ~ STP_TEV, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="STP/TEV", side=1, line=1,cex=0.8)
	mtext(text="R12/TEV", side=2, line=1,cex=0.8)
	plot(R12_TEV ~ SCI, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="SCI", side=1, line=1,cex=0.8)
	mtext(text="R12/TEV", side=2, line=1,cex=0.8)
	plot(R12_TEV ~ MCI, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="MCI", side=1, line=1,cex=0.8)
	mtext(text="R12/TEV", side=2, line=1,cex=0.8)
	plot(R12_TEV ~ p_IA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="% Int Area", side=1, line=1,cex=0.8)
	mtext(text="R12/TEV", side=2, line=1,cex=0.8)
	plot(R12_TEV ~ mean_depth, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="Mean depth", side=1, line=1,cex=0.8)
	mtext(text="R12/TEV", side=2, line=1,cex=0.8)
	plot(R12_TEV ~ CLA_EWA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="CLA/EWA", side=1, line=1,cex=0.8)
	mtext(text="R12/TEV", side=2, line=1,cex=0.8)
	
	plot(STP_TEV ~ SCI, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="SCI", side=1, line=1,cex=0.8)
	mtext(text="STP/TEV", side=2, line=1,cex=0.8)
	plot(STP_TEV ~ MCI, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="MCI", side=1, line=1,cex=0.8)
	mtext(text="STP/TEV", side=2, line=1,cex=0.8)
	plot(STP_TEV ~ p_IA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="% Int Area", side=1, line=1,cex=0.8)
	mtext(text="STP/TEV", side=2, line=1,cex=0.8)
	plot(STP_TEV ~ mean_depth, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="Mean depth", side=1, line=1,cex=0.8)
	mtext(text="STP/TEV", side=2, line=1,cex=0.8)
	plot(STP_TEV ~ CLA_EWA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="CLA/EWA", side=1, line=1,cex=0.8)
	mtext(text="STP/TEV", side=2, line=1,cex=0.8)
	
	plot(SCI ~ MCI, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="MCI", side=1, line=1,cex=0.8)
	mtext(text="SCI", side=2, line=1,cex=0.8)
	plot(SCI ~ p_IA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="% Int Area", side=1, line=1,cex=0.8)
	mtext(text="SCI", side=2, line=1,cex=0.8)
	plot(SCI ~ mean_depth, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="Mean depth", side=1, line=1,cex=0.8)
	mtext(text="SCI", side=2, line=1,cex=0.8)
	plot(SCI ~ CLA_EWA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="CLA/EWA", side=1, line=1,cex=0.8)
	mtext(text="SCI", side=2, line=1,cex=0.8)
	
	plot(MCI ~ p_IA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="% Int Area", side=1, line=1,cex=0.8)
	mtext(text="MCI", side=2, line=1,cex=0.8)
	plot(MCI ~ mean_depth, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="Mean depth", side=1, line=1,cex=0.8)
	mtext(text="MCI", side=2, line=1,cex=0.8)
	plot(MCI ~ CLA_EWA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="CLA/EWA", side=1, line=1,cex=0.8)
	mtext(text="MCI", side=2, line=1,cex=0.8)

	plot(p_IA ~ mean_depth, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="Mean depth", side=1, line=1,cex=0.8)
	mtext(text="% Int Area", side=2, line=1,cex=0.8)
	plot(p_IA ~ CLA_EWA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="CLA/EWA", side=1, line=1,cex=0.8)
	mtext(text="% Int Area", side=2, line=1,cex=0.8)
	
	plot(mean_depth ~ CLA_EWA, data=set1.dat, col=colscheme, pch=16, ylab="", xlab="",xaxt="n",yaxt="n")
	mtext(text="CLA/EWA", side=1, line=1,cex=0.8)
	mtext(text="Mean depth", side=2, line=1,cex=0.8)
	
	plot(0, type="n", ylab="", xlab="", axes=F, main="")
	legend("center", col=colo, legend=leglab, pch=16, bty="n", ncol=2)

}

plotfunc(colscheme=1, leglab=c("All"),colo=1)

plotdiff <- function (group.labs) {

	t1 <- tapply(set1.dat$STP_TEV,group.labs, FUN=mean)
	t2 <- tapply(set1.dat$R12_TEV, group.labs, FUN=mean)
	t3 <- tapply(set1.dat$SCI, group.labs, FUN=mean)
	t4 <- tapply(set1.dat$MCI, group.labs, FUN=mean)
	t5 <- tapply(set1.dat$mean_depth, group.labs, FUN=mean)
	t6 <- tapply(set1.dat$p_IA, group.labs, FUN=mean)
	t7 <- tapply(set1.dat$CLA_EWA, group.labs, FUN=mean)
	
	s1 <- tapply(set1.dat$STP_TEV,group.labs, FUN=sd)
	s2 <- tapply(set1.dat$R12_TEV, group.labs, FUN=sd)
	s3 <- tapply(set1.dat$SCI, group.labs, FUN=sd)
	s4 <- tapply(set1.dat$MCI, group.labs, FUN=sd)
	s5 <- tapply(set1.dat$mean_depth, group.labs, FUN=sd)
	s6 <- tapply(set1.dat$p_IA, group.labs, FUN=sd)
	s7 <- tapply(set1.dat$CLA_EWA, group.labs, FUN=sd)
	
	par(mfrow=c(2,4), mar=c(3,5,2,1))
	
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
	plot(t5, pch=16, ylab="Mean depth", xlab="", ylim=c(min(t5-s5), max(t5+s5)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t5[i]-s5[i], x1=i, y1=t5[i]+s5[i], col=1, lwd=2)
	}
	plot(t6, pch=16, ylab="% Intertidal Area", xlab="", ylim=c(min(t6-s6), max(t6+s6)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t6[i]-s6[i], x1=i, y1=t6[i]+s6[i], col=1, lwd=2)
	}
	plot(t7, pch=16, ylab="CLA/EWA", xlab="", ylim=c(min(t7-s7), max(t7+s7)), cex=1.5)
	for (i in 1:max(group.labs)) {
		segments(x0=i, y0=t7[i]-s7[i], x1=i, y1=t7[i]+s7[i], col=1, lwd=2)
	}
}


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
			STP_TEV.centro <- tapply(train$STP_TEV, cut.train.adj, FUN=mean)
			R12_TEV.centro <- tapply(train$R12_TEV, cut.train.adj, FUN=mean)
			p_IA.centro <- tapply(train$p_IA, cut.train.adj, FUN=mean)
			mean_depth.centro <- tapply(train$mean_depth, cut.train.adj, FUN=mean)
			CLA_EWA.centro <- tapply(train$CLA_EWA, cut.train.adj, FUN=mean)
			
			test.label <- array(NA, dim=c(n.left))
			for(j in 1:n.left) {													# Loop over left out entities
				dist.check <- array(NA, dim=c(ngroups))
				for(i in 1:ngroups) {												# Loop over possible groups for this entity
					dist.check[i] <- (MCI.centro[i] - left.out$MCI[j])^2 + (SCI.centro[i] - left.out$SCI[j])^2 + (STP_TEV.centro[i] - left.out$STP_TEV[j])^2 + (R12_TEV.centro[i] - left.out$R12_TEV[j])^2 + (p_IA.centro[i] - left.out$p_IA[j])^2 + (mean_depth.centro[i] - left.out$mean_depth[j])^2 + (CLA_EWA.centro[i] - left.out$CLA_EWA[j])^2
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

vc <- valid.clust(inp.data=set1.scale, ngroups=4,n.left=20, lim=1)


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
###########################################################################################################################
# Cluster Method 1 - Hierarchical clustering based on Euclidean distances using wards clustering algorithm
###########################################################################################################################

dist.eucl.set1 <- dist(as.matrix(set1.scale[,-1]), method="euclidean")
hc.eucl.set1 <- hclust(dist.eucl.set1, method="ward")
plot(hc.eucl.set1, xlab="Estuary_no")

# determining the number of groups

fold.vec <- c(20,40)
group.vec <- c(2,3,4,5,6,7,8,9,10)
lim.vec <- c(300,300,300,300,300,300,200,100,100)

error.vec <- 0
diff.vec <- 0
fold.out <- 0
group.out <- 0

K <- 1
for ( i in 1:2) {
	for( j in 1:9) {
		vc <- valid.clust(inp.data=set1.scale, ngroups=group.vec[j],n.left=fold.vec[i], lim=lim.vec[j])
		error.vec[K] <- vc$perc
		diff.vec[K] <- vc$dif
		fold.out[K] <- fold.vec[i]
		group.out[K] <- group.vec[j]
		K <- K+1
	}
}

cv.dataframe <- data.frame(fold.out, group.out, error.vec, diff.vec)

par(mfrow=c(2,1), mar=c(5,4,1,2))
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

# 8 and 9 groups are supported

# 8 groups
k8.eucl.set1 <- cutree(hc.eucl.set1, k=8)
plotfunc(colscheme=k8.eucl.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8"), colo=c(1,2,3,4,5,6,7,8))
plotdiff(k8.eucl.set1)
table.meth(k8.eucl.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k8.eucl.set1)

# 9 groups
k9.eucl.set1 <- cutree(hc.eucl.set1, k=9)
plotfunc(colscheme=k9.eucl.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9"), colo=c(1,2,3,4,5,6,7,8,9))
plotdiff(k9.eucl.set1)
table.meth(k9.eucl.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k9.eucl.set1)

###########################################################################################################################################################
# Random Forests clustering - random forest distance coupled with wards clustering algorithm
###########################################################################################################################################################

set1.rf <- est.data[,c(1,6,7,11,12,13,14,15)]
rf.prox.set1 <- randomForest(set1.rf[, -1], type="unsupervised", proximity=TRUE, ntrees=100000)
dist.rf.set1 <- as.dist(sqrt(1-rf.prox.set1$proximity))
hc.rf.set1 <- hclust(dist.rf.set1, method="ward")

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
			rf.cut.tree <- randomForest(cut.tree ~ SCI + MCI + STP_TEV + R12_TEV + p_IA + mean_depth + CLA_EWA, data=rf.dat.test, sampsize=221,mtry=2)
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

rf.clust.ident(data.obj=set1.rf,tree.obj=hc.rf.set1, kclust=c(2:30), nrepeats=1000) 

# 5 groups

k5.rf.set1 <- cutree(hc.rf.set1, k=5)
plotfunc(colscheme=k5.rf.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5"), colo=c(1,2,3,4,5))
plotdiff(k5.rf.set1)
table.meth(k5.rf.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k5.rf.set1)

est.data[k5.rf.set1==1,1:2]
est.data[k5.rf.set1==2,1:2]
est.data[k5.rf.set1==3,1:2]
est.data[k5.rf.set1==4,1:2]
est.data[k5.rf.set1==5,1:2]


# 7 groups

k7.rf.set1 <- cutree(hc.rf.set1, k=7)
plotfunc(colscheme=k7.rf.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7"), colo=c(1,2,3,4,5,6,7))
plotdiff(k7.rf.set1)
table.meth(k7.rf.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k7.rf.set1)

est.data[k7.rf.set1==1,1:2]
est.data[k7.rf.set1==2,1:2]
est.data[k7.rf.set1==3,1:2]
est.data[k7.rf.set1==4,1:2]
est.data[k7.rf.set1==5,1:2]
est.data[k7.rf.set1==6,1:2]
est.data[k7.rf.set1==7,1:2]

# 8 groups

k8.rf.set1 <- cutree(hc.rf.set1, k=8)
plotfunc(colscheme=k8.rf.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8"), colo=c(1,2,3,4,5,6,7,8))
plotdiff(k8.rf.set1)
table.meth(k8.rf.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k8.rf.set1)

est.data[k8.rf.set1==1,1:2]
est.data[k8.rf.set1==2,1:2]
est.data[k8.rf.set1==3,1:2]
est.data[k8.rf.set1==4,1:2]
est.data[k8.rf.set1==5,1:2]
est.data[k8.rf.set1==6,1:2]
est.data[k8.rf.set1==7,1:2]

# 9 groups

k9.rf.set1 <- cutree(hc.rf.set1, k=9)
plotfunc(colscheme=k9.rf.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9"), colo=c(1,2,3,4,5,6,7,8,9))
plotdiff(k9.rf.set1)
table.meth(k9.rf.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k9.rf.set1)

est.data[k9.rf.set1==1,1:2]
est.data[k9.rf.set1==2,1:2]
est.data[k9.rf.set1==3,1:2]
est.data[k9.rf.set1==4,1:2]
est.data[k9.rf.set1==5,1:2]
est.data[k9.rf.set1==6,1:2]
est.data[k9.rf.set1==7,1:2]

# 10 groups

k10.rf.set1 <- cutree(hc.rf.set1, k=10)
plotfunc(colscheme=k10.rf.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9", "Group 10"), colo=c(1,2,3,4,5,6,7,8,9,10))
plotdiff(k10.rf.set1)
table.meth(k10.rf.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k10.rf.set1)

est.data[k10.rf.set1==1,1:2]
est.data[k10.rf.set1==2,1:2]
est.data[k10.rf.set1==3,1:2]
est.data[k10.rf.set1==4,1:2]
est.data[k10.rf.set1==5,1:2]
est.data[k10.rf.set1==6,1:2]
est.data[k10.rf.set1==7,1:2]

# 11 groups

k11.rf.set1 <- cutree(hc.rf.set1, k=11)
plotfunc(colscheme=k11.rf.set1, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9", "Group 10", "Group 11"), colo=c(1,2,3,4,5,6,7,8,9,10,11))
plotdiff(k11.rf.set1)
table.meth(k11.rf.set1, act.dat=as.character(est.data$est.class))
Tree.Construct(k11.rf.set1)

est.data[k11.rf.set1==1,1:2]
est.data[k11.rf.set1==2,1:2]
est.data[k11.rf.set1==3,1:2]
est.data[k11.rf.set1==4,1:2]
est.data[k11.rf.set1==5,1:2]
est.data[k11.rf.set1==6,1:2]
est.data[k11.rf.set1==7,1:2]

############################################################################################################################################################
# Model based clustering
############################################################################################################################################################

MclustBIC.set1 <- mclustBIC(set1.scale[,-1], G=2:30)

# VEV 9 groups

Mclust.data.a1 <- Mclust(set1.scale[,-1], G=2:9)
Mclust.class <- Mclust.data.a1$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9"), colo=c(1,2,3,4,5,6,7,8,9))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))
Tree.Construct(k11.rf.set1)


est.data[Mclust.class==1,1:2]
est.data[Mclust.class==2,1:2]
est.data[Mclust.class==3,1:2]
est.data[Mclust.class==4,1:2]
est.data[Mclust.class==5,1:2]
est.data[Mclust.class==6,1:2]
est.data[Mclust.class==7,1:2]
est.data[Mclust.class==8,1:2]
est.data[Mclust.class==9,1:2]

# VEV 10 groups

Mclust.data.a1 <- Mclust(set1.scale[,-1], G=2:10)
Mclust.class <- Mclust.data.a1$classification
plotfunc(colscheme=Mclust.class, leglab=c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5", "Group 6", "Group 7", "Group 8", "Group 9", "Group 10"), colo=c(1,2,3,4,5,6,7,8,9,10))
plotdiff(Mclust.class)
table.meth(clust.dat=Mclust.class, act.dat=as.character(est.data$est.class))
Tree.Construct(Mclust.class)


est.data[Mclust.class==1,1:2]
est.data[Mclust.class==2,1:2]
est.data[Mclust.class==3,1:2]
est.data[Mclust.class==4,1:2]
est.data[Mclust.class==5,1:2]
est.data[Mclust.class==6,1:2]
est.data[Mclust.class==7,1:2]
est.data[Mclust.class==8,1:2]
est.data[Mclust.class==9,1:2]
est.data[Mclust.class==10,1:2]


