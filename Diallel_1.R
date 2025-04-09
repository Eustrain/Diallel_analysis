##########################
## 
#########################


df<- read.csv("metodo1.csv",header = TRUE)
head(df)
df_1


df0<- df %>% pivot_longer(cols = starts_with("r"),names_to = "rep",values_to = "y") %>% 
  mutate(cross=paste("cr",X,sep = ""), ENV= 1) %>% 
  select(-X,-total,-Mean) %>% 
  mutate(cross=as.factor(cross))


n<- length(levels(as.factor(df0$p1)))
r <- length(levels(as.factor(df0$rep)))

m0 <- lm(y~cross+rep,data = df0)
anv1<- anova(m0)

rSS <- anv1$`Sum Sq`[2]
css <- anv1$`Sum Sq`[1]
ESS <- anv1$`Sum Sq`[3]

rMS <- anv1$`Mean Sq`[2]
cMS <- anv1$`Mean Sq`[1]
EMS <- anv1$`Mean Sq`[3]

data1 <- data.frame(aggregate(y ~ p1:p2, data = df0, FUN= sum))
data1 <- merge(data1,expand.grid(p1= 6,p2=6),
               by=c("p1","p2"),all=TRUE)

data1 <- data1[order(data1$p2,data1$p1),]

mydf <- aggregate(y ~ p2, data1, "c",na.action=na.pass)

myMatrix <- as.matrix(mydf[, -1])
colnames(myMatrix) <- paste("MALE",mydf$p2)
rownames(myMatrix) <- paste("FEMALE",mydf$p2)


pt<- sum(diag(myMatrix))

fc <- pt^2/(n*r)

nSS <- (sum(diag(myMatrix^2))/3)-fc
nMS <-  nSS/(n-1)

rf1t<- sum(myMatrix[upper.tri(myMatrix)])

f1t<- sum(myMatrix[lower.tri(myMatrix)])
xbar <- sum(myMatrix)

CF <- ((f1t+rf1t)^2)/90

sMM<- ((sum(myMatrix[upper.tri(myMatrix)]^2)+sum(myMatrix[lower.tri(myMatrix)]^2))/3)-CF
sMS <- sMM/29

GCF <- (xbar^2)/(36*3)
ncSS <-fc+CF-GCF




fCF<- fit^2/(15*3)
fss<- (sum(myMatrix[lower.tri(myMatrix)]^2))/3-fCF

fms <- fss/(14)

icf <- rf1t^2/(15*3)

iss <- (sum(myMatrix[upper.tri(myMatrix)]^2))/3-icf
ims <- iss/14

fiss <- fCF+icf-CF









