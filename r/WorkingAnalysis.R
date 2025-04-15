library(readr)
library(vegan)
#import Data 
Fish <- read.csv("Data/eDNA_Data.csv")
View(Fish)
Chem18 <- read.csv("Data/Waterchem2018.csv", header=TRUE)
View(Chem18)
head(Chem18)
Nutrients18 <- read.csv("Data/Nutrients2018.csv")
Pesticides18 <- read_csv("Data/Pesticides2018.csv")
Sediments18 <- read_csv("Data/Sediment2018.csv")
Chem20 <- read_csv("Data/Waterchem2020.csv")
Nutrients20 <- read_csv("Data/Nutrients2020.csv")
Pesticides20 <- read_csv("Data/Pesticides2020.csv")
Sediments20 <- read_csv("Data/Sediments2020.csv")

#SHOULD I TRANSPOSE THE EDNA MATRIX? 
#FishT=t(Fish) # transpose matrix
#CAN I STILL USE THE ANALYSES WE DID THE FIRST TIME- SMTH WITH COLS

head(Fish)
Fish=Fish[,-c(1:11)] #drop first 11 columns
tail(Fish) #shows last few rows
View(Fish)
rowSums(Fish) #weird distribution- binomial distribution 

#BELOW PROBS DOESN'T WORK if TRANSPOSED??
source("r/functions/gsindex.R") #knows about function now
#species index for entire dataset
gsindex(rowSums(Fish)) #234 species. suggests missing data as gsindex provides estimated sp pool size. assumes variation in abundace unlike shannons. but not published so don't use
fisher.alpha(rowSums(Fish)) #fishers- off by factor of 6 to gsindex- much more accurate than shannons na dsimpson- low sampling variance. issue is that doesn;t give richness estimate- only abstract estimate
dim(Fish) #shows how many species- diff to above indiactes missing data
colnames(Fish) # not even- some have 8 some have 10
cols=matrix(NA,8,2) # 8 rows by 2 columns
cols[1,]=c(1,10) #populate each batch
cols[2,]=c(11,20)
cols[3,]=c(21,30)
cols[4,]=c(31,34)
cols[5,]=c(35,44)
cols[6,]=c(45,54)
cols[7,]=c(55,64)
cols[8,]=c(65,68)
#keep filling in- see if g value repeat later. Line 14 gsindex sould be larger than below as its for whole world rather than certain places
#breaks it up into sites




#species estimate for eac estuary
g=array() # blank array
bootstrapgsi=array()
r=array()
h=array() #shannons
fa=array()

for(i in 1:8) { #looking through first 2 rows, opening for loop
  n=rowSums(Fish[,cols[i,1]:cols[i,2]]) #make count vector- first to last collumn of row i
  print(n) #goes from 0-10. diversity function won't work on seros so git rid of them
  n=n[n>0] #modify object to only keep numbers greater than 0
  print(n) # now no zeros
  g[i]=gsindex(n) # put first value in first cell, second in 2nd etc
  fa[i]=fisher.alpha(n)
  r[i]= length(n) #tracking the raw number of detected species- what you foudn and not what was there
  f<-n/sum(n) #do before bootstrap to avoid messing with data
  h[i]=-sum(f*log(f)) #runs shannons
  n=sample(n,replace = T) #resampling with replacememt 
  bootstrapgsi[i]=gsindex(n) #runs boottsrao
}
fa
g #suggests varuation bw sites and diversity 
# in wet season- fish dissapear (last 4 number)
#answering q of whether sp richnessestimates are replicable 
plot(g[1:4],g[5:8])
plot(fa,g) #not what is expected. fishers is getting confused by several points- should be linear trend. means either GS or Fishers is wrong.
#gs is geomtric index- good way of estimating species richness
#gs is assuming that sampling is about 1 individual per species (simplification- not true), eac species has own abundance (random set of abundances out in the world), assuming that distribution of abundances is evenley spaced, assumes exponential distribution (chau assumes unifrom dist which makes no sense)
# if above is true then gs is very good richness estimate
#Fishers assumes data follows log series (can be a disiater)
# report both gsindex and fishers. 
#refer to kerr and alroy

#wilcoxon rank sum test to compare diversity index
wilcox.test(fa[1:4], fa[5:8], paired=TRUE) #not signficant
wilcox.test(g[1:4], g[5:8], paired=TRUE) 
#perfect divide insifcates that there is a suggestive difference between the wet and dry season due to low sample size. want at least 20 site vs year combinations (400 samples)




r 
plot(r,g) #they somehwat agree- relatively linear
#could do bootstrap on g by resampling n
#not enough data to do sig testing

bootstrapgsi #suggests not alot of error- only chnages last number- similar to g above otherwise

h #very small number- make exponential below
h=exp(h)

sd(log(g))
sd(log(h)) # inla diff between standard dev
sd(log(r)) #well behaved data

plot(g,h) #linear but flattens out bc limited by h- loose the variation




sum(is.na(Chem18))/length(unlist(Chem18))
sum(is.na(Nutrients18))/length(unlist(Nutrients18))
sum(is.na(Pesticides18))/length(unlist(Pesticides18))
sum(is.na(Sediments18))/length(unlist(Sediments18))

#removed first two collumns + new matrix of numbers
Chem18=as.data.frame(Chem18)
Chem18a=matrix(as.numeric(unlist(Chem18[,-(1:2)])),35,21) #matrix of numners
is.numeric(Chem18a) #make sure its numeric
rownames(Chem18a)=Chem18[,2]
Chem18a=matrix(as.numeric(unlist(Chem18a)),nrow(Chem18a),ncol(Chem18a))

Nutrients18a=Nutrients18[,-(1:2)]
Nutrients18a=matrix(as.numeric(unlist(Nutrients18a)),nrow(Nutrients18a),ncol(Nutrients18a))
Pesticides18a=Pesticides18[,-(1:2)]
Pesticides18a=matrix(as.numeric(unlist(Pesticides18a)),nrow(Pesticides18a),ncol(Pesticides18a))
Sediments18a=Sediments18[,-(1:2)]
Sediments18a=matrix(as.numeric(unlist(Sediments18a)),nrow(Sediments18a),ncol(Sediments18a))
length(Chem18[,2])


#replace NA's with median
clean<- function(x){
  for(i in 1:ncol(x)) {
    x[is.na(x[,i]),i]=median(x[,i],na.rm=TRUE) #which rows in the collumn i have NAs and then remove NAs from caluctaion of medians
  }
  x #returns x
}

#should return 0 NAs
Chem18a=clean(Chem18a)
sum(is.na(Chem18a))/length(unlist(Chem18a))

Nutrients18a=clean(Nutrients18a)
sum(is.na(Nutrients18a))/length(unlist(Nutrients18a))

Pesticides18a=clean(Pesticides18a)
sum(is.na(Pesticides18a))/length(unlist(Pesticides18a))

Sediments18a=clean(Sediments18a)
sum(is.na(Sediments18a))/length(unlist(Sediments18a))




#waterchem Factor analysis
factanal(na.omit(Chem18a),factors=5)
cor(na.omit(Chem18a))
#unqueness tells us if it worked- is there unexplaied variation (doesnt correspond to variables)- bad if 1, 0 is explained
# loading correlation- needs to add up to 1
#GREAT 5 FACTOR analysis
#can do regression with 5 variables
Chem18FA<-factanal(na.omit(Chem18a),factors=5,scores="regression")$scores
plot(Chem18FA,
     main = "Factor Analysis Scores: 2018 Water Chemistry Data",
     xlab = "Factor 1",
     ylab = "Factor 2",
     col = "red") #Do not want all data to be in one corner (clustered) with some outliers
Chem18FA=rbind(Chem18FA[1:10,],NA,Chem18FA[11:33,])
Chem18FA=rbind(Chem18FA,Chem18FA)
summary(lm(log(Fish)~Chem18FA))
Chem18FA
Chem18


#nutrients Factor analysis
factanal(na.omit(Nutrients18a),factors=5)
cor(na.omit(Nutrients18a))
Nutr18FA<-factanal(na.omit(Nutrients18a),factors=5,scores="regression")$scores
plot(Nutr18FA, 
     main = "Factor Analysis Scores: 2018 Nutrient Data",
     xlab = "Factor 1",
     ylab = "Factor 2",
     col = "red")
summary(lm(log(Fish)~Nutr18FA))


#pesticides Factor analysis
factanal(na.omit(Pesticides18a),factors=6)
View(Pesticides18)
cor(na.omit(Pesticides18a))
Pest18FA<-factanal(na.omit(Pesticides18a),factors=5,scores="regression")$scores
plot(Pest18FA)
summary(lm(log(Fish)~Pest18FA))

#sediments 
factanal(na.omit(Sediments18a),factors=3)
cor(na.omit(Sediments18a))
Sed18FA<-factanal(na.omit(Sediments18a),factors=3,scores="regression")$scores
plot(Sed18FA)
summary(lm(log(Fish)~Sed18FA))





SitesSeason=array(dim=68)
SitesSeason[1:10]=1
SitesSeason[11:20]=2
SitesSeason[21:30]=3
SitesSeason[31:34]=4
SitesSeason[35:44]=5
SitesSeason[45:54]=6
SitesSeason[55:64]=7
SitesSeason[65:68]=8

#PCoa with bray curtis
B<-cmdscale(vegdist(t(Fish)))
Season=c(rep(1,34),rep(2,34)) #
Season # usee to make characters different
plot(B,cex=Season) #big circles=2020 data, small=2018
plot(B,col=hsv(h=SitesSeason/9),pch=19) #sites are different colours. whopping huge difference
colnames(Fish)
summary(lm(B[,1]~as.factor(SitesSeason)))
dim(B)        

summary(lm(B[,1]~Chem18FA[,1]+ Nutr18FA[,1]+Sed18FA[,1]+as.factor(SitesSeason)))

row.names(Chem18FA)

#regression for each principle coordinate axis (Pcoa)
summary(lm(B[,1]~Chem18FA+as.factor(SitesSeason)))
summary(lm(B[,2]~Chem18FA+as.factor(SitesSeason)))


