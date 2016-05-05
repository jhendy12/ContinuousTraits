#You can use code you wrote for the correlation exercise here.
setwd("~/GitHub/ContinuousTraits")
source("ContinuousFunctions.R")
library(ape)
library(geiger)
library(phytools)
library(OUwie)

library(phylolm)
library(corHMM)
library(rotl)
library(phylobase)

tree<-pbtree(n=20,scale=1,tip.label=(c("T.ochroleucum","T.palmeri","T.calocephalum","T.nanum","T.pallescens","T.repens","T.amabile","T.decorum","T.wigginsii","T.virginicum","T.vesiculosum","T.vernum","T.uniflorum","T.tomentosa","T.thalii","T.sylvatica","T.alpestre","T.argentinense","T.pratense","T.hybridum")))

con.1<-c(8,  6, 14, 16, 22, 10, 23, 19,  7,  6, 12, 14, 18, 24, 20,  9, 28,  6, 22, 17)
con.2<-c(1, 2, 3, 4, 1, 5, 5, 1, 6, 3, 2, 4, 4, 2, 2, 5, 1, 4, 1, 5)
continuous.data<-cbind.data.frame(row.names=tree$tip.label,con.1,con.2)

#same as last exercise
# using a simulated tree since my data has polytomies
#tree <- read.tree("____PATH_TO_TREE_OR_SOME_OTHER_WAY_OF_GETTING_A_TREE____")
#discrete.data <- read.csv(file="____PATH_TO_DATA_OR_SOME_OTHER_WAY_OF_GETTING_TRAITS____", stringsAsFactors=FALSE) #death to factors.
#continuous.data <- read.csv(file="____PATH_TO_DATA_OR_SOME_OTHER_WAY_OF_GETTING_TRAITS____", stringsAsFactors=FALSE) #death to factors.

cleaned.continuous <- CleanData(tree, continuous.data)
cleaned.discrete <- CleanData(tree, discrete.data)
VisualizeData(tree, cleaned.continuous)
VisualizeData(tree, cleaned.discrete)

#First, start basic. What is the rate of evolution of your trait on the tree? 

BM1 <- fitContinuous(cleaned.continuous$phy, cleaned.continuous$data, model="BM")
print(paste("The rate of evolution is", BM1$con.1$opt$sigsq, "in units of", "???"))
#Important: What are the rates of evolution? In what units?
OU1 <- fitContinuous(cleaned.continuous$phy, cleaned.continuous$data, model="OU")

#what does this do?  I can't use mfcol for some reason
par(mfcol(c(1,2)))


plot(cleaned.continuous$phy, show.tip.label=FALSE)
ou.tree <- rescale(tree, model="OU", OU1$con.1$opt$alpha)
plot(ou.tree)
#How are the trees different?
# they have different branch lengths


#Compare trees
AIC.BM1 <- BM1$con.1$opt$aic
AIC.OU1 <- OU1$con.1$opt$aic
delta.AIC.BM1 <- AIC.BM1-min(c(AIC.BM1,AIC.OU1))
delta.AIC.OU1 <- AIC.OU1-min(c(AIC.BM1,AIC.OU1))

print(delta.AIC.BM1)
print(delta.AIC.OU1)
#whichever = 0 is the better model (lower AIC)


#OUwie runs:
#This takes longer than you may be used to. 
#We're a bit obsessive about doing multiple starts and in general
#performing a thorough numerical search. It took you 3+ years
#to get the data, may as well take an extra five minutes to 
#get an accurate answer

#First, we need to assign regimes. The way we do this is with ancestral state estimation of a discrete trait.
#We can do this using ace() in ape, or similar functions in corHMM or diversitree. Use only one discrete char
one.discrete.char <- cleaned.discrete$data[,1]
#this only works as TRUE
reconstruction.info <- ace(one.discrete.char, cleaned.continuous$phy, type="discrete", method="ML", CI=TRUE)
best.states <- apply(reconstruction.info$lik.anc, 1, which.max)




# This part was very confusing to get it to run.  I used a lot of help to make this part work.




#NOW ADD THESE AS NODE LABELS TO YOUR TREE

labeled.tree <- cleaned.continuous$phy
labeled.tree$node.label<-best.states
tips<-rownames(cleaned.continuous$data)
c.continuous<-data.frame(tips,cleaned.discrete$data[,1],cleaned.continuous$data[,1])
colnames(c.continuous) <- c("tips","regime","data")

nodeBased.OUMV <- OUwie(labeled.tree, c.continuous,model="OUMV", simmap.tree=FALSE, diagn=FALSE)
print(nodeBased.OUMV)
#What do the numbers mean?
# I'm not 100% sure on this.
# it seems they are giving an AIC score for this model along with rates of some sort?

#Now run all OUwie models:
models <- c("BM1","BMS","OU1","OUMV","OUMA","OUMVA")
# I cannot get this to work at all
# I was told to take out "OUM" since it won't run and ends up crashing the rest of the models
results <- lapply(models, RunSingleOUwieModel, phy=labeled.tree, data=c.continuous)

AICc.values<-sapply(results, "[[", "AICc")
names(AICc.values)<-models
AICc.values<-AICc.values-min(AICc.values)


print(AICc.values) #The best model is the one with smallest AICc score

best<-results[[which.min(AICc.values)]] #store for later

print(best) #prints info on best model







#We get SE for the optima (see nodeBased.OUMV$theta) but not for the other parameters. Let's see how hard they are to estimate. 
#First, look at ?OUwie.fixed to see how to calculate likelihood at a single point.
?OUwie.fixed

#Next, keep all parameters but alpha at their maximum likelihood estimates (better would be to fix just alpha and let the others optimize given this constraint, but this is harder to program for this class). Try a range of alpha values and plot the likelihood against this.
alpha.values<-seq(from= 1 , to= 10 , length.out=50)

#keep it simple (and slow) and do a for loop:
likelihood.values <- rep(NA, length(alpha.values))
for (iteration in sequence(length(alpha.values))) {
  likelihood.values[iteration] <- OUwie.fixed(labeled.tree, c.continuous, model="OUMV", alpha=rep(alpha.values[iteration],2), sigma.sq=best$solution[2,], theta=best$theta[,1])$loglik
}

#Error message Error in Rate.mat[1, 1:k] <- alpha : 
#number of items to replace is not a multiple of replacement length


plot(x=best$solution[1,1] , y= best$loglik, xlab="Alpha Values", ylab="Loglik", type="l", bty="n")
#filled in x and y with same thing as below

points(x=best$solution[1,1], y=best$loglik, pch=16, col="red")
text(x=best$solution[1,1], y=best$loglik, "unconstrained best", pos=4, col="red")

#a rule of thumb for confidence for likelihood is all points two log likelihood units worse than the best value. Draw a dotted line on the plot to show this
abline(h=best$loglik-2, lty="dotted") #Two log-likelihood 



#everything below this works
#Now, let's try looking at both theta parameters at once, keeping the other parameters at their MLEs
require("akima")

nreps<-400
theta1.points<-c(best$theta[1,1], rnorm(nreps-1, best$theta[1,1], 5*best$theta[1,2])) #center on optimal value, have extra variance
theta2.points<-c(best$theta[2,1], rnorm(nreps-1, best$theta[2,1], 5*best$theta[2,2])) #center on optimal value, have extra variance
likelihood.values<-rep(NA,nreps)

for (iteration in sequence(nreps)) {
  likelihood.values[iteration] <- OUwie.fixed(labeled.tree, c.continuous, model="OUMV", alpha=best$solution[1,], sigma.sq=best$solution[2,], theta=c(theta1.points[iteration], theta2.points[iteration]))$loglik
}
#think of how long that took to do 400 iterations. Now remember how long the search took (longer).

likelihood.differences<-(-(likelihood.values-max(likelihood.values)))

#We are interpolating here: contour wants a nice grid. But by centering our simulations on the MLE values, we made sure to sample most thoroughly there
interpolated.points<-interp(x=theta1.points, y=theta2.points, z= likelihood.differences, linear=FALSE, extrap=TRUE, xo=seq(min(theta1.points), max(theta1.points), length = 400), yo=seq(min(theta2.points), max(theta2.points), length = 400))

contour(interpolated.points, xlim=range(c(theta1.points, theta2.points)),ylim=range(c(theta1.points, theta2.points)), xlab="Theta 1", ylab="Theta 2", levels=c(2,5,10),add=FALSE,lwd=1, bty="n", asp=1)

points(x=best$theta[1,1], y=best$theta[2,1], col="red", pch=16)

points(x=trait$X[which(trait$Reg==1)],y=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==1))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 1, plotted along x axis
points(y=trait$X[which(trait$Reg==2)],x=rep(min(c(theta1.points, theta2.points)), length(which(trait$Reg==2))), pch=18, col=rgb(0,0,0,.3)) #the tip values in regime 2, plotted along y axis


#The below only works if the discrete trait rate is low, so you have a good chance of estimating where the state is.
#If it evolves quickly, hard to estimate where the regimes are, so some in regime 1 are incorrectly mapped in
#regime 2 vice versa. This makes the models more similar than they should be.
#See Revell 2013, DOI:10.1093/sysbio/sys084 for an exploration of this effect.
library(phytools)
trait.ordered<-data.frame(c.continuous[,2], c.continuous[,2],row.names=c.continuous[,1])
trait.ordered<- trait.ordered[tree$tip.label,]
z<-trait.ordered[,1]
names(z)<-rownames(trait.ordered)
tree.mapped<-make.simmap(tree,z,model="ER",nsim=1)
leg<-c("black","red")
names(leg)<-c(1,2)
plotSimmap(tree.mapped,leg,pts=FALSE,ftype="off", lwd=1)

simmapBased<-OUwie(tree.mapped,trait,model="OUMV", simmap.tree=TRUE, diagn=FALSE)
print(simmapBased)
#How does this compare to our best model from above? Should they be directly comparable?
print(best)