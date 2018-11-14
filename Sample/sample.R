library(foreach)
library(doParallel)
registerDoParallel(40)

source("GEECC_GEECO_20171220_8.r")

load("AMASS.test.V2.rdata")

outfile.margin = "Sample.Margin.txt"
errorfile.margin = "Sample.Margin.error"

outfile.inter = "Sample.inter.txt"
errorfile.inter = "Sample.inter.error"

fam =  as.integer(factor(sarcoid.ped$FID,labels = 1:length(unique(sarcoid.ped$FID))))

Kin.invert = solve(sarcoid.Kin)


initials.margin=foreach(j = 1:ncol(sarcoid.g),.combine=cbind) %dopar% {
	initial<-glm(as.factor(sarcoid.ped$PHENOTYPE) ~ as.factor(sarcoid.g[,j]) + sarcoid.covar,family=binomial) #initial value
        b0 = initial$coeff
	return(b0)
}

initials.margin = as.matrix(initials.margin)


GEE_Cpp_CO(outfile.margin,errorfile.margin, initials.margin,sarcoid.ped$PHENOTYPE,sarcoid.g,sarcoid.covar,fam,sarcoid.Kin,Kin.invert,1)


initials.inter =foreach(j = 1:ncol(sarcoid.g),.combine=cbind) %dopar% {
	initial<-glm(as.factor(sarcoid.ped$PHENOTYPE) ~ as.factor(sarcoid.g[,j]) * as.factor(sarcoid.ped$Exposure) + sarcoid.covar ,family=binomial) #initial value
        b0 = initial$coeff
	b.temp = b0[7]
	b0[5:7] = b0[4:6]
	b0[4] = b.temp
	#if(j==1) initials = b0 else initials = cbind(initials,b0)
	return(b0)
}
initials.inter = as.matrix(initials.inter)

GEE_Cpp_CC(outfile.inter,errorfile.inter, initials.inter,sarcoid.ped$PHENOTYPE,sarcoid.ped$Exposure,sarcoid.g,sarcoid.covar,fam,sarcoid.Kin,Kin.invert,1)

#save(list=c("sarcoid.ped","sarcoid.d","sarcoid.g","sarcoid.covar","sarcoid.Kin"),file = "AMASS.test.V2.rdata")

