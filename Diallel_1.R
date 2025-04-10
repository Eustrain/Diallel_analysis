############################################
##  this  run the diallel analysis 1     ##
###########################################

library(tidyverse)
library(ggplot2)



df0 <- read.csv("diallel_1.csv")


n<- length(levels(as.factor(df0$p1)))
r <- length(levels(as.factor(df0$rep)))
c <- length(levels(as.factor(df0$cross))) 
s <- c-n

###################################################
### anova #######################################
m0 <- lm(y~cross+rep,data = df0)
anv1<- anova(m0)

rSS <- anv1$`Sum Sq`[2]
css <- anv1$`Sum Sq`[1]
ESS <- anv1$`Sum Sq`[3]

rMS <- anv1$`Mean Sq`[2]
cMS <- anv1$`Mean Sq`[1]
EMS <- anv1$`Mean Sq`[3]

###################################################
## make matrix to diallelic anova
data1 <- data.frame(aggregate(y ~ p1:p2, data = df0, FUN= sum))
data1 <- merge(data1,expand.grid(p1= 6,p2=6),
               by=c("p1","p2"),all=TRUE)

data1 <- data1[order(data1$p2,data1$p1),]

mydf <- aggregate(y ~ p2, data1, "c",na.action=na.pass)
myMatrix <- as.matrix(mydf[, -1])
colnames(myMatrix) <- paste("MALE",mydf$p2)
rownames(myMatrix) <- paste("FEMALE",mydf$p2)

###############################################
fc <- sum(diag(myMatrix))^2/(n*r)
nSS <- (sum(diag(myMatrix^2))/3)-fc
nMS <-  nSS/(n-1)
##############################################
rf1t<- sum(myMatrix[upper.tri(myMatrix)])
f1t<- sum(myMatrix[lower.tri(myMatrix)])
xbar <- sum(myMatrix)
##############################################
CF <- ((f1t+rf1t)^2)/(s*r)

sMM<- ((sum(myMatrix[upper.tri(myMatrix)]^2)+sum(myMatrix[lower.tri(myMatrix)]^2))/r)-CF
sMS <- sMM/(s-1)
###############################################
GCF <- (xbar^2)/(c*r)
ncSS <-fc+CF-GCF
#################################################
f<- length(myMatrix[lower.tri(myMatrix)])
fCF<- f1t^2/(f*r)
fss<- (sum(myMatrix[lower.tri(myMatrix)]^2))/r-fCF
fms <- fss/(f-1)
####################################################
l<- length(myMatrix[upper.tri(myMatrix)])
icf <- rf1t^2/(l*r)
iss <- (sum(myMatrix[upper.tri(myMatrix)]^2))/r-icf
ims <- iss/(l-1)
fiss <- fCF+icf-CF
###################################################
####### anova 

fv <- c("Rep","progenies","Parents(n)","Crosses","F1´s","rF1s","f vs l","Parents vs crosses","Error")
df <- c((r-1),(c-1),(n-1),(s-1),(f-1),(l-1),1,1,((c-1)*(r-1)))
SS <- c(rSS,css,nSS,sMM,fss,iss,fiss,ncSS,ESS)
anova_df <- data.frame("Source of variation"= fv,"df" = df,"SS"=SS,"MS"=SS/df )

####################################################
### estimation components genetics
###################################################

## make matrix to diallelic anova
data2 <- data.frame(aggregate(y ~ p1:p2, data = df0, FUN= mean))
data2 <- merge(data2,expand.grid(p1= n,p2=n),
               by=c("p1","p2"),all=TRUE)

data2 <- data2[order(data2$p2,data2$p1),]

mydf2 <- aggregate(y ~ p2, data2, "c",na.action=na.pass)
myMatrix2 <- t(round(as.matrix(mydf2[, -1]),1))
colnames(myMatrix2) <- paste("MALE",mydf2$p2)
rownames(myMatrix2) <- paste("FEMALE",mydf2$p2)

Xi. <-rowSums(myMatrix2) 
Xj. <- colSums(myMatrix2)

Xbar <- sum(myMatrix2)
G<- (1/(2*n))*(sum((Xi.+Xj.)^2))

U<- 1/(n^2)*(sum(myMatrix2)^2)
##### ss of gca
SSg <- G-(2*U)
 
sum(myMatrix2 * (myMatrix2 + t(myMatrix2)))/2  -G+U

