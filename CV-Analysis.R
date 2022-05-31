
#################
################# This file performs the cross validation analysis, outputs results,
################# and creates figures.

cor.set1 = read.csv(file =   "data/Correlationsset1.csv", check.names = FALSE, header = TRUE)[,-1]

cor.set2 = read.csv(file =   "data/Correlationsset2.csv", check.names = FALSE, header = TRUE)[,-1]

cor.set3 = read.csv(file =   "data/Correlationsset3.csv", check.names = FALSE, header = TRUE)[,-1]
cor.set4 = read.csv(file =   "data/Correlations S-10-9-02-truncated.csv", check.names = FALSE, header = TRUE)[,-1]

cor.set5 = read.csv(file =   "data/Correlations K-9-9-truncated.csv", check.names = FALSE, header = TRUE)[,-1]
cor.set6 = read.csv(file =   "data/Correlations S-10-9-10-truncated.csv", check.names = FALSE, header = TRUE)[,-1]

cor.set1$set = "Rot 1"
cor.set2$set = "Rot 2"
cor.set3$set = "Rot 3"

library('scales')
library('ggplot2')
library('lme4')
library('psych')
library('reshape2')
library('MASS')
devtools::load_all('MixMatrix_0.1.0')

# library('MixMatrix')
library('CholWishart')
library('dplyr')
library('cowplot')
invlogit = function(x) exp(x)/(1+exp(x))
logit = function(x) log(x/(1-x))

#totalcorset = rbind(cor.set1,cor.set2,cor.set3,cor.set5,cor.set4, cor.set6)
totalcorset = rbind(cor.set1,cor.set5, cor.set4, cor.set6)
totalcorset$zfive = fisherz(totalcorset$`5-10`)
totalcorset$zten = fisherz(totalcorset$`10-20`)
#totalcorset$zfive = (totalcorset$`5-10`)
#totalcorset$zten = (totalcorset$`10-20`)
totalcorset$repetition = factor(totalcorset$set)
totalcorset$img = factor(totalcorset$img)
matchset = totalcorset[totalcorset$match == "match",]
matchset = droplevels(matchset)
fisher_trans = function() trans_new("fisher", function(x) psych::fisherz(x), 
                                    function(x) psych::fisherz2r(x))


keep.vars = c("knife","img","match","zfive","zten","repetition")
corset.shaped = reshape(totalcorset[,keep.vars], idvar = c("knife","match","repetition"), timevar = "img", direction = "wide")
corset.subset = corset.shaped[corset.shaped$repetition == "Rot 1",]
matchlabels = corset.shaped$match[corset.shaped$repetition == "Rot 1"]
#matchlabels = corset.shaped$match
dims = dim(corset.shaped)
df = 3
knifematrix = matrix(unlist(corset.shaped[,4:21]),ncol = 18)
knifearray = array(t(knifematrix),dim = c(2,9,dims[1]))[,c(1,3,5,7,9),]
knifematrix.part = knifematrix[corset.shaped$repetition == "Rot 1",]

knifematrix.partarray = array(t(knifematrix.part),dim = c(2,9,81))[,c(1,3,5,7,9),]
matchingknives = knifematrix.partarray[,,matchlabels == "match"]
nonmatchingknives = knifematrix.partarray[,,matchlabels != "match"]





### sets 4 and 6 are the ones for investigation with 9
### 
### 
### 

nameone = c("E10", "T10", "E48", "E39", "E36", "E34", "E28", "E27", "E24")
nametwo = c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10")
namethree = c("S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20")
namefour = c("T01", "T02", "T03", "T04", "T05", "T06", "T07", "T08", "T09")

sets = c("K-9-9-01"     , "Rot 1",         "S-10-09-02-01", "S-10-9-10-01" )

tensets = c("S-10-09-02-01", "S-10-9-10-01" )

train_test = function(train, test, matchlabels, df){

  matchingknives = train[,,matchlabels == "match"]
  nonmatchingknives = train[,,matchlabels != "match"]
  
  matching.fit = MLmatrixt(df = df,  matchingknives,
                                              row.mean = TRUE, col.variance = "AR(1)")
  unmatched.fit = MLmatrixt(df = df,  nonmatchingknives,
                                                row.mean = TRUE, col.variance = "AR(1)")
  
  classes =  dmatrixt(test,mean = matching.fit$mean,
                      U = matching.fit$U, V = matching.fit$var*matching.fit$V, log = TRUE, df = df) - 
    dmatrixt(test,mean = unmatched.fit$mean,
             U = unmatched.fit$U, V = unmatched.fit$var*unmatched.fit$V, log = TRUE, df = df)
  
  classes
}

traintest_name = function(matrix, array, name, matchlabels, df){
  knifeout = !grepl(name, matrix[, "knife"])
  train = array[,,knifeout]
  test = array[,,!knifeout]
  train_test(train, test, matchlabels[knifeout], df)
  
}

setoneres = numeric(0)
classone = numeric(0)
for(name in nameone) {
 boole = grepl(name, corset.subset[, "knife"])
 classtmp = corset.subset$match[boole]
  tmpres= traintest_name(corset.subset, knifematrix.partarray, name, matchlabels, df = 3)
 setoneres = c(setoneres, tmpres)
 classone = c(classone, classtmp)
}

fivecvdf = tibble(matching = classone, logodds = setoneres, df = 3, set = "Rot 1")

# dfchoose is 5, 10, 15, 20, after initially doing for 3

for (dfchoose in c(5, 10, 15, 20, 30)){
for(set in sets){
  if(set == "Rot 1") namelist = nameone
  if(set == "K-9-9-01") namelist = namefour
  if(set == "S-10-09-02-01") namelist = nametwo
  if(set == "S10-9-10-01") namelist = namethree
  knifematrix = matrix(unlist(corset.shaped[,4:21]),ncol = 18)
  knifearray = array(t(knifematrix),dim = c(2,9,dims[1]))[,c(1,3,5,7,9),]
  knifematrix.part = knifematrix[corset.shaped$repetition == set,]
  
  knifematrix.partarray = array(t(knifematrix.part),dim = c(2,9,dim(knifematrix.part)[1]))[,c(1,3,5,7,9),]
  corset.subset = corset.shaped[corset.shaped$repetition == set,]
  matchlabels = corset.shaped$match[corset.shaped$repetition == set]
  
  setoneres = numeric(0)
  classone = numeric(0)
  for(name in namelist) {
    boole = grepl(name, corset.subset[, "knife"])
    classtmp = corset.subset$match[boole]
    tmpres= traintest_name(corset.subset, knifematrix.partarray, name, matchlabels, df = dfchoose)
    setoneres = c(setoneres, tmpres)
    classone = c(classone, classtmp)
  }
  
  tempdf = tibble(matching = classone, logodds = setoneres, df = dfchoose, set = set)
#  if(set != "Rot 1") 
    fivecvdf = rbind(fivecvdf, tempdf)
}
}

qplot(matching, logodds, fill = factor(df), data = fivecvdf, geom = "boxplot")+
  theme_bw()+
  xlab(NULL)+
  ylab("Log odds (base 10)")+
  labs(fill = "Degrees of \nfreedom")+
  NULL

ggsave("five-image-CV.pdf", width = 5, height = 4)

### stp here


write.csv(fivecvdf, "five-image-CV-results.csv", row.names = FALSE)




############ all 9 images


tencvdf = tibble(matching = character(0), logodds = numeric(0), df = numeric(0), set = character(0))

### manually go 3, 5, 10, 15, 20
for (dfchoose in c(3, 5, 10, 15, 20, 30)){
for(set in tensets){
  if(set == "Rot 1") namelist = nameone
  if(set == "K-9-9-01") namelist = namefour
  if(set == "S-10-09-02-01") namelist = nametwo
  if(set == "S10-9-10-01") namelist = namethree
  knifematrix = matrix(unlist(corset.shaped[,4:21]),ncol = 18)
  knifearray = array(t(knifematrix),dim = c(2,9,dims[1]))
  knifematrix.part = knifematrix[corset.shaped$repetition == set,]
  
  knifematrix.partarray = array(t(knifematrix.part),dim = c(2,9,dim(knifematrix.part)[1]))
  corset.subset = corset.shaped[corset.shaped$repetition == set,]
  matchlabels = corset.shaped$match[corset.shaped$repetition == set]
  
  setoneres = numeric(0)
  classone = numeric(0)
  for(name in namelist) {
    boole = grepl(name, corset.subset[, "knife"])
    classtmp = corset.subset$match[boole]
    tmpres= traintest_name(corset.subset, knifematrix.partarray, name, matchlabels, df = dfchoose)
    setoneres = c(setoneres, tmpres)
    classone = c(classone, classtmp)
  }
  
  tempdf = tibble(matching = classone, logodds = setoneres, df = dfchoose, set = set)
  #  if(set != "Rot 1") 
  tencvdf = rbind(tencvdf, tempdf)
}
}
  
qplot(matching, logodds/log(10), fill = factor(df), data = tencvdf, geom = "boxplot")+
  theme_bw()+
  xlab(NULL)+
  ylab("Log odds (base 10)")+
  labs(fill = "Degrees of \nfreedom")+
  NULL

ggsave("nine-image-CV.pdf", width = 5, height = 4)

write.csv(tencvdf, "nine-image-CV-results.csv", row.names = FALSE)






threecvdf = tibble(matching = character(0), logodds = numeric(0), df = numeric(0), set = character(0))

### manually go 3, 5, 10, 15, 20
for (dfchoose in c(3, 5, 10, 15, 20, 30)){
for(set in sets){
  if(set == "Rot 1") namelist = nameone
  if(set == "K-9-9-01") namelist = namefour
  if(set == "S-10-09-02-01") namelist = nametwo
  if(set == "S10-9-10-01") namelist = namethree
  knifematrix = matrix(unlist(corset.shaped[,4:21]),ncol = 18)
  knifearray = array(t(knifematrix),dim = c(2,9,dims[1]))[,c(1,5,9),]
  knifematrix.part = knifematrix[corset.shaped$repetition == set,]
  
  knifematrix.partarray = array(t(knifematrix.part),dim = c(2,9,dim(knifematrix.part)[1]))[,c(1,5,9),]
  corset.subset = corset.shaped[corset.shaped$repetition == set,]
  matchlabels = corset.shaped$match[corset.shaped$repetition == set]
  
  setoneres = numeric(0)
  classone = numeric(0)
  for(name in namelist) {
    boole = grepl(name, corset.subset[, "knife"])
    classtmp = corset.subset$match[boole]
    tmpres= traintest_name(corset.subset, knifematrix.partarray, name, matchlabels, df = dfchoose)
    setoneres = c(setoneres, tmpres)
    classone = c(classone, classtmp)
  }
  
  tempdf = tibble(matching = classone, logodds = setoneres, df = dfchoose, set = set)
  #  if(set != "Rot 1") 
  threecvdf = rbind(threecvdf, tempdf)
}
}
  
qplot(matching, logodds, fill = factor(df), data = threecvdf, geom = "boxplot")+
  theme_bw()+
  xlab(NULL)+
  ylab("Log odds (base 10)")+
  labs(fill = "Degrees of \nfreedom")+
  NULL

ggsave("three-image-CV.pdf", width = 5, height = 4)



write.csv(threecvdf, "three-image-CV-results.csv", row.names = FALSE)







#ninecvdf = read.csv('nine-image-CV-results.csv')
ninecvdf$images = "9 images"

#threecvdf = read.csv('three-image-CV-results.csv')
threecvdf$images = "3 images"

#fivecvdf = read.csv('five-image-CV-results.csv')
fivecvdf$images = "5 images"

combcvdf = rbind(ninecvdf, fivecvdf, threecvdf)

library(ggplot2)


qplot(matching, logodds, fill = factor(df), data = combcvdf, geom = "boxplot", facets = ~images)+
  theme_bw()+
  xlab(NULL)+
  ylab(expression(paste(Log[10]," odds of match"))) +
  labs(fill = "Degrees of \nfreedom")+
  geom_hline(yintercept=0)+
  scale_fill_brewer(palette = "Dark2")+
  NULL

ggsave("combined-image-CV.pdf", width = 6, height = 4)

