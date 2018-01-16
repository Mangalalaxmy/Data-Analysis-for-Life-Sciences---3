install.packages("devtools")
library(devtools)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)##this loads the three tables

head(geneExpression)
head(geneAnnotation)
head(sampleInfo)
sum(sampleInfo$date=="2005-06-27")

table(geneAnnotation$CHR=="chrY")

sampleInfo[sampleInfo$date=="2005-06-10",]
geneAnnotation[geneAnnotation$SYMBOL=="ARPC1A",]
sub1=geneExpression["200950_at","GSM136727.CEL.gz"]
sub1

i = which(geneAnnotation$SYMBOL=="ARPC1A")
j = which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]

meds = apply(geneExpression,2,median)
median(meds)

tester = function(e,group){
  x = e[group==1]
  y = e[group==0]
  return (t.test( x, y)$p.value)
}
g = factor(sampleInfo$group)
pvals = apply(geneExpression,1,tester, group=g)
min(pvals)

set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
  control = sample(population[,1],12)
  treatment = sample(population[,1],12)
  t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)
mean(pvals<0.05)
mean(pvals<0.01)

set.seed(100)
pvalsnull = replicate(20, {
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  t.test(cases,controls)$p.val
})
sum(pvalsnull<0.05)

set.seed(100)
pvals005 = replicate(1000,{
  pvalsnull = replicate(20, {
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  t.test(cases,controls)$p.val
  })
  sum(pvalsnull<=0.05)
})
table(pvals005)
mean(pvals005)
(1000-354)/1000

alphas <- seq(0,0.25,0.01)
par(mfrow=c(2,2))
for (m in c(2,10,100, 1000)){
  plot(alphas, alphas/m - (1-(1-alphas)^(1/m)))
  abline(h=0,col=2,lty=2)
}

fwer <- 0.05/8793
set.seed(1)
pvalsfwer = replicate(10000,{
  pvals = runif(8793,0,1)
  sum(pvals < fwer)
})
mean(pvalsfwer>0)

k = (1-(1-0.05)^(1/8793))
set.seed(1)
pvalsfwer = replicate(10000,{
  pvals = runif(8793,0,1)
  sum(pvals < k)
})
mean(pvalsfwer>0)


set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- alpha/m
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)  

library(devtools)
library(rafalib)
install_github("genomicsclass/GSE5859Subset")
install_bioc("genefilter")
install_bioc("qvalue")

library(GSE5859Subset)
data(GSE5859Subset)

head(sampleInfo)
sampleInfo$group
library(genefilter)
?rowttests
dim(geneExpression)

g = factor(sampleInfo$group)
alpha = 0.05
pvals = rowttests(geneExpression, g)$p.value
sum(pvals < alpha)

fweralpha = 0.05
m = 8793
BCalpha = fweralpha/m
sum(pvals < BCalpha)

?p.adjust
fdr = p.adjust(pvals, method = "fdr")
sum(fdr < 0.05)

h = hist(pvals)
lambda = 0.3
pi0 = sum(pvals > lambda)/ ((1-lambda)*m)
abline(h = pi0)

library(qvalue)
res = qvalue(pvals)
qvals = res$qvalues
sum(qvals < 0.05)

res$pi0
plot(pvals, qvals)

set.seed(1)
n <- 24
m <- 8793
B = 1000
delta <- 2
positives <- 500
g = factor(rep(c(0,1), each=12))
result = replicate(B, {
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.value
  FP1 = sum(pvals[-(1:positives)]<0.05/m) 
  FP1
})
mean(result/(m-positives))

set.seed(1)
n <- 24
m <- 8793
B = 1000
delta <- 2
positives <- 500
g = factor(rep(c(0,1), each=12))
result = replicate(B, {
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.value
  FP1 = sum(pvals[-(1:positives)]<0.05/m) 
  FN1 = sum(pvals[1:positives]>0.05/m)
  FN1
})
mean(result/(positives))

set.seed(1)
n <- 24
m <- 8793
B = 1000
delta <- 2
positives <- 500
g = factor(rep(c(0,1), each=12))
result = replicate(B, {
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.value
  BHvals = p.adjust(pvals, method="fdr")
  FP1 = sum(BHvals[-(1:positives)]<0.05) 
  FN1 = sum(BHvals[1:positives]>0.05)
  c(FP1,FN1)
})
mean(result[1,]/(m-positives))
mean(result[2,]/(positives))

set.seed(1)
n <- 24
m <- 8793
B = 1000
delta <- 2
positives <- 500
g = factor(rep(c(0,1), each=12))
result = replicate(B, {
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.value
  BHvals = p.adjust(pvals, method="fdr")
  res = qvalue(pvals)
  qvals = res$qvalues
  FP1 = sum(qvals[-(1:positives)]<0.05) 
  FN1 = sum(qvals[1:positives]>0.05)
  c(FP1,FN1)
})
mean(result[1,]/(m-positives))
mean(result[2,]/(positives))

dbinom(2, 4, 0.49)

dbinom(4, 10, 0.49)

1-pbinom(10,20,0.4) 
pbinom(10,20,0.4, lower.tail = FALSE)

1-dbinom(0, 189000000, (1/175223510))

1-pbinom(1, 189000000, (1/175223510))

pbinom(9,20,0.4) - pbinom(7,20,0.4)

b <- (9 - 20*.4)/sqrt(20*.4*.6)
a <- (7 - 20*.4)/sqrt(20*.4*.6)
pnorm(b)-pnorm(a)

pbinom(450,1000,0.4) - pbinom(350,1000,0.4)
b <- (450 - 1000*.4)/sqrt(1000*.4*.6)
a <- (350 - 1000*.4)/sqrt(1000*.4*.6)
pnorm(b)-pnorm(a)
0.9987609-0.9987512

exact = dbinom(k,N,p)
a <- (k+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (k-0.5 - N*p)/sqrt(N*p*(1-p))
approx = pnorm(a) - pnorm(b)
Ns <- c(5,10,30,100)
ps <- c(0.01,0.10,0.5,0.9,0.99)
mypar(4,5)
for (n in Ns){
  for (p in ps) {
    exact = dbinom((1:(n-1)),n,p)
    a <- ((1:(n-1))+0.5 - n*p)/sqrt(n*p*(1-p))
    b <- ((1:(n-1))-0.5 - n*p)/sqrt(n*p*(1-p))
    approx = pnorm(a) - pnorm(b)
    plot(exact,approx,main=paste("N =",n," p = ",p))
    abline(0,1)
  }
}

N <- 189000000
p <- 1/175223510
dbinom(2,N,p)
a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)
dpois(2,N*p)
ppois(1,N*p, lower.tail = FALSE)

library(devtools)
install_github("genomicsclass/dagdata")
library(dagdata)
data(hcmv)
library(rafalib)
mypar()
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")
breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
hist(counts)
probs <- dpois(counts,4)
likelihood <- prod(probs)
likelihood
logprobs <- dpois(counts,4,log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood
lambdas = seq(0,15,len=300)
logl = function(lambda, x) {
  logprobs = dpois(x,lambda,log=TRUE)
  sum(logprobs)
}
est = sapply(lambdas, function(lambda) logl(lambda, counts))
plot(lambdas, est)
mle = lambdas[which.max(est)]
mle
abline(v=mle, col="red")

breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab=)
bl = binLocation[which.max(counts)]
abline(v=bl)

lambda = mean(counts[- which.max(counts)])
1-ppois(13,lambda)

0.05/57

ps <- (seq(along=counts) - 0.5)/length(counts)
lambda <- mean( counts[ -which.max(counts)])
poisq <- qpois(ps,lambda)
qqplot(poisq,counts)
abline(0,1)

library(devtools)
install_github("genomicsclass/tissuesGeneExpression",force=TRUE)
library(tissuesGeneExpression)
data("tissuesGeneExpression")
library(genefilter)
y = e[,which(tissue=="endometrium")]

variance = rowVars(y)
qqnorm(variance)
qqline(variance)
qqnorm(sqrt(variance))
qqline(sqrt(variance))

library(limma)
estimates=fitFDist(variance,14)
print( estimates$scale )

tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
##this shows us files
filenames <- unzip(tmpfile,list=TRUE)
players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

library(dplyr)
filter(players,yearID==2012) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)
data = filter(players,yearID==2010 | yearID==2011 | yearID==2012) %>% mutate(AVG=H/AB) %>% filter(AB>500) %>% select(AVG)
mean(data$AVG)
sd(data$AVG)
hist(data$AVG)
qqnorm(data$AVG)
qqline(data$AVG)
sqrt(.45*(1-0.45)/20)
B = (0.11)^2/((0.11^2)+(0.027^2))
estavg = B*(0.275)+(1-B)*0.450
estavg

## to install:
library(rafalib)
install_bioc("SpikeInSubset")
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)
pData(rma95)
g <- factor(rep(0:1,each=3))
spike <- rownames(y) %in% colnames(pData(rma95))

pvals = rowttests(y,g)$p.value
sigpval = pvals < 0.01
print(mean(!spike[sigpval]))

library(genefilter)
sds <- rowSds(y[,g==0])
index <- paste0( as.numeric(spike), as.numeric(sigpval))
index <- factor(index,levels=c("11","01","00","10"),labels=c("TP","FP","TN","FN"))
boxplot(split(sds,index))

library(limma)
fit <- lmFit(y, design=model.matrix(~ g))
colnames(coef(fit))
fit <- eBayes(fit)
sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)
LIM = range( c(posteriorSD,sampleSD))
plot(sampleSD, posteriorSD,ylim=LIM,xlim=LIM)
abline(0,1)
abline(v=sqrt(fit$s2.prior))

library(limma)
fit = lmFit(y, design=model.matrix(~ g))
fit = eBayes(fit)
##second coefficient relates to diffences between group
pvals = fit$p.value[,2] 
index = pvals < 0.01
mean(!spike[index])

source("http://www.bioconductor.org/biocLite.R")
biocLite("SpikeInSubset")
library(SpikeInSubset)
data(mas133)
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")
mean(e[,1]< k & e[,2]< k)

plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)
sd(e[,2]-e[,1])
sum(abs(e[,2]-e[,1])>1)


