install.packages("readr")
install.packages("vegan")
install.packages("paran")
library(paran)
library(readr)
library(vegan)
library(lme4)
library(lmerTest)
library(MuMIn)


#import Data 
Fish <- read.csv("Data/eDNA_Data_2018+2020.csv")
Chem18 <- read.csv("Data/Waterchem2018.csv", header=TRUE)
Nutrients18 <- read.csv("Data/Nutrients2018.csv",header=TRUE)
Pesticides18 <- read.csv("Data/Pesticides2018.csv",header=TRUE)
Sediments18 <- read.csv("Data/Sediment2018.csv",header=TRUE)
Chem20 <- read.csv("Data/Waterchem2020.csv",header=TRUE)
Nutrients20 <- read.csv("Data/Nutrients2020.csv",header=TRUE)
Pesticides20 <- read.csv("Data/Pesticides2020.csv",header=TRUE)
Sediments20 <- read.csv("Data/Sediments2020.csv",header=TRUE)
View(Chem18)



head(Fish)
tail(Fish) #shows last few rows
View(Fish)
Fish1=Fish[,-c(1:3)] #drop first 11 column- DO I STILL NEED THIS??
rowSums(Fish1) #weird distribution- binomial distribution 

#BELOW PROBS DOESN'T WORK if TRANSPOSED??
source("r/functions/gsindex.R") #knows about function now
#species index for entire dataset
gsindex(colSums(Fish1)) #234 species. suggests missing data as gsindex provides estimated sp pool size. assumes variation in abundace unlike shannons. but not published so don't use
fisher.alpha(colSums(Fish1)) #fishers- off by factor of 6 to gsindex- much more accurate than shannons na dsimpson- low sampling variance. issue is that doesn;t give richness estimate- only abstract estimate
#alpha and richness have different sclaed
#report thatv sampling of the species pool is pretty good (divide 194/250= 80% of specis pool so bpretty good pick up) but eDNA is not perfect
dim(Fish) #shows how many species- diff to above indiactes missing data

#NEED THIS?
rownames(Fish1)=paste(Fish[,1],Fish[,3]) # adds row names

cols=matrix(NA,8,2) # 8 rows by 2 columns
cols[1,]=c(1,10) #populate each batch
cols[2,]=c(11,20)
cols[3,]=c(21,30)
cols[4,]=c(31,34)
cols[5,]=c(35,44)
cols[6,]=c(45,54)
cols[7,]=c(55,64)
cols[8,]=c(65,68)

cols

#keep filling in- see if g value repeat later. Line 14 gsindex sould be larger than below as its for whole world rather than certain places
#breaks it up into sites




#species estimate for eac estuary
g=array() # blank array
bootstrapgsi=array()
r=array()
h=array() #shannons
fa=array()

for(i in 1:8) { #looking through first 2 rows, opening for loop
  n=colSums(Fish1[cols[i,1]:cols[i,2],]) #make count vector- first to last collumn of row i
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
fa # sensible numbers
g #suggests varuation bw sites and diversity 
# in wet season- fish dissapear (last 4 number)
#answering q of whether sp richnessestimates are replicable 
#report FA and g
plot(g[1:4],g[5:8])
plot(fa,g) #not what is expected. fishers is getting confused by several points- should be linear trend. means either GS or Fishers is wrong.
#gs is geomtric index- good way of estimating species richness
#gs is assuming that sampling is about 1 individual per species (simplification- not true), eac species has own abundance (random set of abundances out in the world), assuming that distribution of abundances is evenley spaced, assumes exponential distribution (chau assumes unifrom dist which makes no sense)
# if above is true then gs is very good richness estimate
#Fishers assumes data follows log series (can be a disiater)
# report both gsindex and fishers. 
#refer to kerr and alroy

r #raw richness
plot(r,g,xlim=c(50,250), ylim=c(50,250))  #they somehwat agree- relatively linear
abline(0,1) #x=creates line of unity
#hada real effect 
#could do bootstrap on g by resampling n
#not enough data to do sig testing


#wilcoxon rank sum test to compare diversity index
wilcox.test(fa[1:4], fa[5:8], paired=TRUE) #not signficant
wilcox.test(g[1:4], g[5:8], paired=TRUE) 
wilcox.test(r[1:4], r[5:8], paired=TRUE) 
#0.125 is mnathetically lowest number so very clean seperation as this is the lower limit for a 4-p[air comaprison- worked perfectly asd far as it can
t.test(log(fa[1:4]), log(fa[5:8]), paired=TRUE)
t.test(log(g[1:4]), log(g[5:8]), paired=TRUE) 
t.test(log(r[1:4]), log(r[5:8]), paired=TRUE) #r= raw richness
#data is 0-bounded and would be skewed if we had lots of points so logged the data, allows for a bell curve to be created
#very signficant diffeence between wet and dry season, 2020 signficantly lower
#report wilcozx and t.test
#GSI and fishers modifications did change the outcome by greatly increasing the number of species and decreasing the p vlue changed in ttest. shows asampling is good but not perfect


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
Chem18a=matrix(as.numeric(unlist(Chem18[,-(1:2)])),34,21) #matrix of numners
is.numeric(Chem18a) #make sure its numeric
rownames(Chem18a)=Chem18[,2]
Chem18a=matrix(as.numeric(unlist(Chem18a)),nrow(Chem18a),ncol(Chem18a))

Nutrients18a=Nutrients18[,-(1:2)]
Nutrients18a=matrix(as.numeric(unlist(Nutrients18a)),nrow(Nutrients18a),ncol(Nutrients18a))
Pesticides18a=Pesticides18[,-(1:2)]
Pesticides18a=matrix(as.numeric(unlist(Pesticides18a)),nrow(Pesticides18a),ncol(Pesticides18a))
Sediments18a=Sediments18[,-(1:2)]
Sediments18a=matrix(as.numeric(unlist(Sediments18a)),nrow(Sediments18a),ncol(Sediments18a))

Chem20a=Chem20[,-(1:2)]
Chem20a=matrix(as.numeric(unlist(Chem20a)),nrow(Chem20a),ncol(Chem20a))
Nutrients20a=Nutrients20[,-(1:2)]
Nutrients20a=matrix(as.numeric(unlist(Nutrients20a)),nrow(Nutrients20a),ncol(Nutrients20a))
Pesticides20a=Pesticides20[,-(1:2)]
Pesticides20a=matrix(as.numeric(unlist(Pesticides20a)),nrow(Pesticides20a),ncol(Pesticides20a))
Sediments20a=Sediments20[,-(1:2)]
Sediments20a=matrix(as.numeric(unlist(Sediments20a)),nrow(Sediments20a),ncol(Sediments20a))




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

Chem20a=clean(Chem20a)
sum(is.na(Chem20a))/length(unlist(Chem20a))

Nutrients20a=clean(Nutrients20a)
sum(is.na(Nutrients20a))/length(unlist(Nutrients20a))

Pesticides20a=clean(Pesticides20a)
sum(is.na(Pesticides20a))/length(unlist(Pesticides20a))

Sediments20a=clean(Sediments20a)
sum(is.na(Sediments20a))/length(unlist(Sediments20a))


###2018##
#waterchem Factor analysis
paran(Chem18a) #horn's parallel analysis and tells you how many compenents you have-use in factor ana
factanal(na.omit(Chem18a),factors=5) #won;t run unbless you have at least 5
cor(na.omit(Chem18a)) # see if raw data has problem-i.e. correlation of 1 will crash it
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
Chem18FA=rbind(Chem18FA[1:10,],NA,Chem18FA[11:33,]) #some bad data so we inserted a row so it would lineup
Chem18FA=rbind(Chem18FA,Chem18FA)
summary(lm(log(Fish1)~Chem18FA))
Chem18FA
Chem18
dim(Chem18a)
dim(na.omit(Chem18a))

dim(Nutrients18a)
dim(na.omit(Nutrients18a))

#nutrients Factor analysis
paran(Nutrients18a)
factanal(na.omit(Nutrients18a),factors=4)
cor(na.omit(Nutrients18a))
Nutr18FA<-factanal(na.omit(Nutrients18a),factors=5,scores="regression")$scores
plot(Nutr18FA, 
     main = "Factor Analysis Scores: 2018 Nutrient Data",
     xlab = "Factor 1",
     ylab = "Factor 2",
     col = "red")
summary(lm(log(g)~Nutr18FA))

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


    ###2020##
# waterchem Factor analysis
#undo binding changes
paran(ChemB)
factanal(na.omit(ChemB),factors=5)
cor(na.omit(ChemB))
ChemBFA<-factanal(na.omit(ChemB),factors=5,scores="regression")$scores
plot(Chem20FA,
     main = "Factor Analysis Scores: 2020 Water Chemistry Data",
     xlab = "Factor 1",
     ylab = "Factor 2",
     col = "red") #Do not want all data to be in one corner (clustered) with some outliers
summary(lm(log(Fish)~Chem20FA))


#nutrients Factor analysis
factanal(na.omit(Nutrients20a),factors=5)
cor(na.omit(Nutrients20a))
Nutr20FA<-factanal(na.omit(Nutrients20a),factors=5,scores="regression")$scores
plot(Nutr20FA, 
     main = "Factor Analysis Scores: 2020 Nutrient Data",
     xlab = "Factor 1",
     ylab = "Factor 2",
     col = "red")
summary(lm(log(Fish)~Nutr20FA))

#pesticides Factor analysis
factanal(na.omit(Pesticides20a),factors=6)
View(Pesticides20)
cor(na.omit(Pesticides20a))
Pest20FA<-factanal(na.omit(Pesticides20a),factors=5,scores="regression")$scores
plot(Pest20FA)
summary(lm(log(Fish)~Pest20FA))

#sediments 
factanal(na.omit(Sediments20a),factors=3)
cor(na.omit(Sediments20a))
Sed20FA<-factanal(na.omit(Sediments20a),factors=3,scores="regression")$scores
plot(Sed20FA)
summary(lm(log(Fish)~Sed20FA))



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
B<-cmdscale(vegdist(Fish1))
Season=c(rep(1,34),rep(2,34)) #
Site=SitesSeason
Site[35:68]=Site[35:68]-4
Site
Season # usee to make characters different
plot(B,cex=Season) #big circles=2020 data, small=2018
plot(B,col=hsv(h=SitesSeason/9),pch=19) #sites are different colours. whopping huge difference
colnames(Fish)
summary(lm(B[,1]~as.factor(SitesSeason))) # see if fish comp is affected by environemtn
dim(B)        

summary(lm(B[,1]~Chem18FA[,1]+ Nutr18FA[,1]+Sed18FA[,1]+as.factor(SitesSeason)))

row.names(Chem18FA)

#regression for each principle coordinate axis (Pcoa)
summary(lm(B[,1]~Chem18FA+as.factor(SitesSeason)))
summary(lm(B[,2]~Chem18FA+as.factor(SitesSeason)))
ChemFA=rbind(Chem18FA,Chem18FA)

#
m=lmer(B[,1]~ChemFA+(1|Site)+(1|Season))
summary(m)
rand(m) #season importnat and site isn't based on pvalues
ranef(m) #Site doesn;t do much and season does alot
r.squaredGLMM(m) #no variance from checm but alot of variance for whole model and chem is hudgely important 


#add each matrices- niutrients etc 
#undo changes 



