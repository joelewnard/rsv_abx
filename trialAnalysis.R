
load('dat.Rdata')
dx = read.csv(file='dx.csv',header=T)

resp = dat$dx%in%dx$dx[which(dx$resp==1)]
lrti = dat$dx%in%dx$dx[which(dx$lrti==1)]
aom = dat$dx%in%dx$dx[which(dx$om==1)]

####################################
#### all new abx prescriptions #####
####################################
ids = unique(dat$id)

dat.fn = function(criterion=rep(1,dim(dat)[1]),endTime=1e6){
  
  dat$lastTime =  dat$lastTime+1 ### consistent  with +1 for drug dates
  dat$lastTime[dat$lastTime>endTime] = endTime
  dat$startTime[dat$startTime>endTime] = NA
  dat$startTime[dat$startTime<0] = NA

  out = c()
  for (i in 1:length(ids)){
  
    country = dat$country[which(dat$id==ids[i])][1]
    child = dat$child[which(dat$id==ids[i])][1]
    mom = dat$mom[which(dat$id==ids[i])][1]
    sex = dat$sex[which(dat$id==ids[i])][1]
    momAge = dat$momAge[which(dat$id==ids[i])][1]
    ppChild = dat$ppChild[which(dat$id==ids[i])][1]
    ppMom = dat$ppMom[which(dat$id==ids[i])][1]
    gestRand = dat$gestRand[which(dat$id==ids[i])][1]
    gestBirth = dat$gestBirth[which(dat$id==ids[i])][1]
    actArm = dat$actArm[which(dat$id==ids[i])][1]
    intendArm = dat$intendArm[which(dat$id==ids[i])][1]
    smoke = dat$smoke[which(dat$id==ids[i])][1]
    alc2y = dat$alc2y[which(dat$id==ids[i])][1]
    alcNow = dat$alcNow[which(dat$id==ids[i])][1]
    recDrug = dat$recDrug[which(dat$id==ids[i])][1]
    
    dates = unique(dat$startTime[which(dat$id==ids[i]&criterion==1&dat$startTime<dat$lastTime)])
    dates = dates[is.na(dates)==F&dates>0]
    
    firstTime = 0
    lastTime = dat$lastTime[which(dat$id==ids[i])][1]
  
    obs = 0
    if (length(dates)>0){
      for (j in 1:length(dates)){
        obs = c(obs,dates[j])
      }
    }
    obs = c(obs,lastTime)
    outcome = c(0,rep(1,length(obs)-2),0)
    
    if (obs[length(obs)]==obs[length(obs)-1]){
      obs = obs[1:(length(obs)-1)]
      outcome = outcome[1:(length(obs)-1)]
    }
    
    newrows = cbind(ids[i],
                    country,child,mom,sex,momAge,ppChild,ppMom,gestRand,gestBirth,actArm,intendArm,
                    smoke,alc2y,alcNow,recDrug,
                    obs[1:(length(obs)-1)],obs[2:length(obs)],outcome[2:length(obs)])
    
    out = rbind(out,newrows)
    print(i)
  }
  
  out = as.data.frame(out)
  names(out) = c('id','country','child','mom','sex','momAge','ppChild','ppMom','gestRand','gestBirth',
                 'actArm','intendArm',
                 'smoke','alc2y','alcNow','recDrug',
                 'start','end','outcome')
  
  out$child = as.logical(out$child)
  out$mom = as.logical(out$mom)
  out$momAge = as.numeric(out$momAge)
  out$ppChild = as.logical(out$ppChild)
  out$ppMom = as.logical(out$ppMom)
  out$gestRand = as.numeric(out$gestRand)
  out$gestBirth = as.numeric(out$gestBirth)
  out$start = as.numeric(out$start)
  out$end = as.numeric(out$end)
  out$outcome = as.numeric(out$outcome)
  out$lmic = out$country%in%c('AR','BD','MX','PH','ZA')
  out$hic = out$country%in%c('AU','CL','ES','NZ','UK','US')
  return(out)
}

datAll90 = dat.fn(endTime=90)
datAllLrti90 = dat.fn(lrti==1,endTime=90)
datAllAom90 = dat.fn(aom==1,endTime=90)
datAllResp90 = dat.fn(resp==1,endTime=90)


datAll180 = dat.fn(endTime=180)
datAllLrti180 = dat.fn(lrti==1,endTime=180)
datAllAom180 = dat.fn(aom==1,endTime=180)
datAllResp180 = dat.fn(resp==1,endTime=180)


datAllTot = dat.fn()
datAllLrtiTot = dat.fn(lrti==1)
datAllAomTot = dat.fn(aom==1)
datAllRespTot = dat.fn(resp==1)




datAntib90 = dat.fn(dat$drugClass1=='antibiotic',endTime=90)

datPen90 = dat.fn(dat$drugClass3=='penicillin',endTime=90)
datCeph90 = dat.fn(dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'),endTime=90)
datBeta90 = dat.fn(dat$drugClass2=='betaLactam',endTime=90)
datMacro90 = dat.fn(dat$drugClass3=='macrolide',endTime=90)
datAmino90 = dat.fn(dat$drugClass3=='aminoglycoside',endTime=90)
datTb90 = dat.fn(dat$drugClass3=='antiTB',endTime=90)
datOthAntib90 = dat.fn(dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F),endTime=90)

datOth90 = dat.fn(dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'),endTime=90)
datAntiV90 = dat.fn(dat$drugClass1=='antiviral',endTime=90)
datAntiF90 = dat.fn(dat$drugClass1=='antifungal',endTime=90)


datAntibLrti90 = dat.fn(lrti==1&dat$drugClass1=='antibiotic',endTime=90)

datPenLrti90 = dat.fn(lrti==1&dat$drugClass3=='penicillin',endTime=90)
datCephLrti90 = dat.fn(lrti==1&dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'),endTime=90)
datBetaLrti90 = dat.fn(lrti==1&dat$drugClass2=='betaLactam',endTime=90)
datMacroLrti90 = dat.fn(lrti==1&dat$drugClass3=='macrolide',endTime=90)
datAminoLrti90 = dat.fn(lrti==1&dat$drugClass3=='aminoglycoside',endTime=90)
datOthAntibLrti90 = dat.fn(lrti==1&dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F),endTime=90)

datOthLrti90 = dat.fn(lrti==1&dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'),endTime=90)
datAntiVLrti90 = dat.fn(dat$drugClass1=='antiviral'&lrti==1,endTime=90)
datAntiFLrti90 = dat.fn(dat$drugClass1=='antifungal'&lrti==1,endTime=90)


datAntibAom90 = dat.fn(aom==1&dat$drugClass1=='antibiotic',endTime=90)

datPenAom90 = dat.fn(aom==1&dat$drugClass3=='penicillin',endTime=90)
datCephAom90 = dat.fn(aom==1&dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'),endTime=90)
datBetaAom90 = dat.fn(aom==1&dat$drugClass2=='betaLactam',endTime=90)
datMacroAom90 = dat.fn(aom==1&dat$drugClass3=='macrolide',endTime=90)
datAminoAom90 = dat.fn(aom==1&dat$drugClass3=='aminoglycoside',endTime=90)
datOthAntibAom90 = dat.fn(aom==1&dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F),endTime=90)

datOthAom90 = dat.fn(aom==1&dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'),endTime=90)


datAntibResp90 = dat.fn(resp==1&dat$drugClass1=='antibiotic',endTime=90)

datPenResp90 = dat.fn(resp==1&dat$drugClass3=='penicillin',endTime=90)
datCephResp90 = dat.fn(resp==1&dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'),endTime=90)
datBetaResp90 = dat.fn(resp==1&dat$drugClass2=='betaLactam',endTime=90)
datMacroResp90 = dat.fn(resp==1&dat$drugClass3=='macrolide',endTime=90)
datAminoResp90 = dat.fn(resp==1&dat$drugClass3=='aminoglycoside',endTime=90)
datOthAntibResp90 = dat.fn(resp==1&dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F),endTime=90)

datOthResp90 = dat.fn(resp==1&dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'),endTime=90)




datAntibTot = dat.fn(dat$drugClass1=='antibiotic')

datPenTot = dat.fn(dat$drugClass3=='penicillin')
datCephTot = dat.fn(dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'))
datBetaTot = dat.fn(dat$drugClass2=='betaLactam')
datMacroTot = dat.fn(dat$drugClass3=='macrolide')
datAminoTot = dat.fn(dat$drugClass3=='aminoglycoside')
datOthAntibTot = dat.fn(dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F))

datOthTot = dat.fn(dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'))
datAntiVTot = dat.fn(dat$drugClass1=='antiviral')
datAntiFTot = dat.fn(dat$drugClass1=='antifungal')



datAntibLrtiTot = dat.fn(lrti==1&dat$drugClass1=='antibiotic')

datPenLrtiTot = dat.fn(lrti==1&dat$drugClass3=='penicillin')
datCephLrtiTot = dat.fn(lrti==1&dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'))
datBetaLrtiTot = dat.fn(lrti==1&dat$drugClass2=='betaLactam')
datMacroLrtiTot = dat.fn(lrti==1&dat$drugClass3=='macrolide')
datAminoLrtiTot = dat.fn(lrti==1&dat$drugClass3=='aminoglycoside')
datOthAntibLrtiTot = dat.fn(lrti==1&dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F))

datOthLrtiTot = dat.fn(lrti==1&dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'))
datAntiVLrtiTot = dat.fn(dat$drugClass1=='antiviral'&lrti==1)
datAntiFLrtiTot = dat.fn(dat$drugClass1=='antifungal'&lrti==1)



datAntibAomTot = dat.fn(aom==1&dat$drugClass1=='antibiotic')

datPenAomTot = dat.fn(aom==1&dat$drugClass3=='penicillin')
datCephAomTot = dat.fn(aom==1&dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'))
datBetaAomTot = dat.fn(aom==1&dat$drugClass2=='betaLactam')
datMacroAomTot = dat.fn(aom==1&dat$drugClass3=='macrolide')
datAminoAomTot = dat.fn(aom==1&dat$drugClass3=='aminoglycoside')
datOthAntibAomTot = dat.fn(aom==1&dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F))

datOthAomTot = dat.fn(aom==1&dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'))


datAntibRespTot = dat.fn(resp==1&dat$drugClass1=='antibiotic')

datPenRespTot = dat.fn(resp==1&dat$drugClass3=='penicillin')
datCephRespTot = dat.fn(resp==1&dat$drugClass3%in%c('cephalosporin1','cephalosporin2','cephalosporin3','cephalosporin4'))
datBetaRespTot = dat.fn(resp==1&dat$drugClass2=='betaLactam')
datMacroRespTot = dat.fn(resp==1&dat$drugClass3=='macrolide')
datAminoRespTot = dat.fn(resp==1&dat$drugClass3=='aminoglycoside')
datOthAntibRespTot = dat.fn(resp==1&dat$drugClass2!='betaLactam'&(dat$drugClass3%in%c('macrolide','aminoglycoside')==F))

datOthRespTot = dat.fn(resp==1&dat$drugClass1%in%c('antifungal','antiprotozoan','antiviral'))



save(datAll90,file='datAll90.Rdata'); save(datAllTot,file='datAllTot.Rdata')
save(datAllLrti90,file='datAllLrti90.Rdata'); save(datAllLrtiTot,file='datAllLrtiTot.Rdata')
save(datAntib90,file='datAntib90.Rdata'); save(datAntibTot,file='datAntibTot.Rdata')
save(datAntibLrti90,file='datAntibLrti90.Rdata'); save(datAntibLrtiTot,file='datAntibLrtiTot.Rdata')
save(datPen90,file='datPen90.Rdata'); save(datPenTot,file='datPenTot.Rdata')
save(datPenLrti90,file='datPenLrti90.Rdata'); save(datPenLrtiTot,file='datPenLrtiTot.Rdata')
save(datCeph90,file='datCeph90.Rdata'); save(datCephTot,file='datCephTot.Rdata')
save(datCephLrti90,file='datCephLrti90.Rdata'); save(datCephLrtiTot,file='datCephLrtiTot.Rdata')
save(datMacro90,file='datMacro90.Rdata'); save(datMacroTot,file='datMacroTot.Rdata')
save(datMacroLrti90,file='datMacroLrti90.Rdata'); save(datMacroLrtiTot,file='datMacroLrtiTot.Rdata')
save(datAmino90,file='datAmino90.Rdata'); save(datAminoTot,file='datAminoTot.Rdata')
save(datAminoLrti90,file='datAminoLrti90.Rdata'); save(datAminoLrtiTot,file='datAminoLrtiTot.Rdata')
save(datOthAntib90,file='datOthAntib90.Rdata'); save(datOthAntibTot,file='datOthAntibTot.Rdata')
save(datOthAntibLrti90,file='datOthAntibLrti90.Rdata'); save(datOthAntibLrtiTot,file='datOthAntibLrtiTot.Rdata')
save(datOth90,file='datOth90.Rdata'); save(datOthTot,file='datOthTot.Rdata')
save(datOthLrti90,file='datOthLrti90.Rdata'); save(datOthLrtiTot,file='datOthLrtiTot.Rdata')

save(datAllAom90,file='datAllAom90.Rdata'); save(datAllAomTot,file='datAllAomTot.Rdata')
save(datAntibAom90,file='datAntibAom90.Rdata'); save(datAntibAomTot,file='datAntibAomTot.Rdata')
save(datPenAom90,file='datPenAom90.Rdata'); save(datPenAomTot,file='datPenAomTot.Rdata')
save(datCephAom90,file='datCephAom90.Rdata'); save(datCephAomTot,file='datCephAomTot.Rdata')
save(datMacroAom90,file='datMacroAom90.Rdata'); save(datMacroAomTot,file='datMacroAomTot.Rdata')
save(datAminoAom90,file='datAminoAom90.Rdata'); save(datAminoAomTot,file='datAminoAomTot.Rdata')
save(datOthAntibAom90,file='datOthAntibAom90.Rdata'); save(datOthAntibAomTot,file='datOthAntibAomTot.Rdata')
save(datOthAom90,file='datOthAom90.Rdata'); save(datOthAomTot,file='datOthAomTot.Rdata')

save(datAllResp90,file='datAllResp90.Rdata'); save(datAllRespTot,file='datAllRespTot.Rdata')
save(datAntibResp90,file='datAntibResp90.Rdata'); save(datAntibRespTot,file='datAntibRespTot.Rdata')
save(datPenResp90,file='datPenResp90.Rdata'); save(datPenRespTot,file='datPenRespTot.Rdata')
save(datCephResp90,file='datCephResp90.Rdata'); save(datCephRespTot,file='datCephRespTot.Rdata')
save(datMacroResp90,file='datMacroResp90.Rdata'); save(datMacroRespTot,file='datMacroRespTot.Rdata')
save(datAminoResp90,file='datAminoResp90.Rdata'); save(datAminoRespTot,file='datAminoRespTot.Rdata')
save(datOthAntibResp90,file='datOthAntibResp90.Rdata'); save(datOthAntibRespTot,file='datOthAntibRespTot.Rdata')
save(datOthResp90,file='datOthResp90.Rdata'); save(datOthRespTot,file='datOthRespTot.Rdata')

save(datAntiVLrtiTot,file='datAntiVLrtiTot.Rdata')
save(datAntiFLrtiTot,file='datAntiFLrtiTot.Rdata')
save(datAntiVTot,file='datAntiVTot.Rdata')
save(datAntiFTot,file='datAntiFTot.Rdata')
save(datAntiVLrti90,file='datAntiVLrti90.Rdata')
save(datAntiFLrti90,file='datAntiFLrti90.Rdata')
save(datAntiV90,file='datAntiV90.Rdata')
save(datAntiF90,file='datAntiF90.Rdata')

################################################################################
###### Analysis here ###########################################################
################################################################################



load('datAll90.Rdata'); load('datAllLrti90.Rdata'); load('datAllTot.Rdata'); load('datAllLrtiTot.Rdata')
load('datAntib90.Rdata'); load('datAntibLrti90.Rdata'); load('datAntibTot.Rdata'); load('datAntibLrtiTot.Rdata')
load('datPen90.Rdata'); load('datPenLrti90.Rdata'); load('datPenTot.Rdata'); load('datPenLrtiTot.Rdata')
load('datCeph90.Rdata'); load('datCephLrti90.Rdata'); load('datCephTot.Rdata'); load('datCephLrtiTot.Rdata')
load('datMacro90.Rdata'); load('datMacroLrti90.Rdata'); load('datMacroTot.Rdata'); load('datMacroLrtiTot.Rdata')
load('datAmino90.Rdata'); load('datAminoLrti90.Rdata'); load('datAminoTot.Rdata'); load('datAminoLrtiTot.Rdata')
load('datOthAntib90.Rdata'); load('datOthAntibLrti90.Rdata'); load('datOthAntibTot.Rdata'); load('datOthAntibLrtiTot.Rdata')
load('datOth90.Rdata'); load('datOthLrti90.Rdata'); load('datOthTot.Rdata'); load('datOthLrtiTot.Rdata')

load('datAllAom90.Rdata');   load('datAllAomTot.Rdata')
load('datAntibAom90.Rdata');   load('datAntibAomTot.Rdata')
load('datPenAom90.Rdata');   load('datPenAomTot.Rdata')
load('datCephAom90.Rdata');   load('datCephAomTot.Rdata')
load('datMacroAom90.Rdata');   load('datMacroAomTot.Rdata')
load('datAminoAom90.Rdata');   load('datAminoAomTot.Rdata')
load('datOthAntibAom90.Rdata');  load('datOthAntibAomTot.Rdata')
load('datOthAom90.Rdata');  load('datOthAomTot.Rdata')

load(file='datAntiVLrtiTot.Rdata')
load(file='datAntiFLrtiTot.Rdata')
load(file='datAntiVTot.Rdata')
load(file='datAntiFTot.Rdata')
load(file='datAntiVLrti90.Rdata')
load(file='datAntiFLrti90.Rdata')
load(file='datAntiV90.Rdata')
load(file='datAntiF90.Rdata')

library(survival)

childMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAll90); 1-summary(childMod90)[[9]][1,c(1,4,3)]
childModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllLrti90); 1-summary(childModLrti90)[[9]][1,c(1,4,3)]

childMod90hic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAll90); 1-summary(childMod90hic)[[9]][1,c(1,4,3)]
childModLrti90hic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllLrti90); 1-summary(childModLrti90hic)[[9]][1,c(1,4,3)]

childMod90lmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAll90); 1-summary(childMod90lmic)[[9]][1,c(1,4,3)]
childModLrti90lmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllLrti90); 1-summary(childModLrti90lmic)[[9]][1,c(1,4,3)]



childModResp90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllResp90); 1-summary(childModResp90)[[9]][1,c(1,4,3)]
childModResp90hic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllResp90); 1-summary(childModResp90hic)[[9]][1,c(1,4,3)]
childModResp90lmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllResp90); 1-summary(childModResp90lmic)[[9]][1,c(1,4,3)]

childModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllAom90); 1-summary(childModAom90)[[9]][1,c(1,4,3)]
childModAom90hic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllAom90); 1-summary(childModAom90hic)[[9]][1,c(1,4,3)]
childModAom90lmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllAom90); 1-summary(childModAom90lmic)[[9]][1,c(1,4,3)]


childModRespTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllRespTot); 1-summary(childModRespTot)[[9]][1,c(1,4,3)]
childModRespTothic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllRespTot); 1-summary(childModRespTothic)[[9]][1,c(1,4,3)]
childModRespTotlmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllRespTot); 1-summary(childModRespTotlmic)[[9]][1,c(1,4,3)]

childModAomTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllAomTot); 1-summary(childModAomTot)[[9]][1,c(1,4,3)]
childModAomTothic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllAomTot); 1-summary(childModAomTothic)[[9]][1,c(1,4,3)]
childModAomTotlmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllAomTot); 1-summary(childModAomTotlmic)[[9]][1,c(1,4,3)]









childModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllTot); 1-summary(childModTot)[[9]][1,c(1,4,3)]
childModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllLrtiTot); 1-summary(childModLrtiTot)[[9]][1,c(1,4,3)]
childModAomTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAllAomTot); 1-summary(childModAomTot)[[9]][1,c(1,4,3)]

childModTothic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllTot); 1-summary(childModTothic)[[9]][1,c(1,4,3)]
childModLrtiTothic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllLrtiTot); 1-summary(childModLrtiTothic)[[9]][1,c(1,4,3)]
childModAomTothic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&hic==1,data=datAllAomTot); 1-summary(childModAomTothic)[[9]][1,c(1,4,3)]

childModTotlmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllTot); 1-summary(childModTotlmic)[[9]][1,c(1,4,3)]
childModLrtiTotlmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllLrtiTot); 1-summary(childModLrtiTotlmic)[[9]][1,c(1,4,3)]
childModAomTotlmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child&lmic==1,data=datAllAomTot); 1-summary(childModAomTotlmic)[[9]][1,c(1,4,3)]







motherModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAllTot); 1-summary(motherModTot)[[9]][1,c(1,4,3)]
motherModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAllLrtiTot); 1-summary(motherModLrtiTot)[[9]][1,c(1,4,3)]

motherModTotHic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom&hic==1,data=datAllTot); 1-summary(motherModTotHic)[[9]][1,c(1,4,3)]
motherModLrtiTotHic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom&hic==1,data=datAllLrtiTot); 1-summary(motherModLrtiTotHic)[[9]][1,c(1,4,3)]

motherModTotLmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom&lmic==1,data=datAllTot); 1-summary(motherModTotLmic)[[9]][1,c(1,4,3)]
motherModLrtiTotLmic = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom&lmic==1,data=datAllLrtiTot); 1-summary(motherModLrtiTotLmic)[[9]][1,c(1,4,3)]





childAntibMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntib90); 1-summary(childAntibMod90)[[9]][1,c(1,4,3)]

childPenMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datPen90); 1-summary(childPenMod90)[[9]][1,c(1,4,3)]
childCephMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datCeph90); 1-summary(childCephMod90)[[9]][1,c(1,4,3)]
childMacroMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datMacro90); 1-summary(childMacroMod90)[[9]][1,c(1,4,3)]
childAminoMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAmino90); 1-summary(childAminoMod90)[[9]][1,c(1,4,3)]
childOthAntibMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthAntib90); 1-summary(childOthAntibMod90)[[9]][1,c(1,4,3)]

childOthMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOth90); 1-summary(childOthMod90)[[9]][1,c(1,4,3)]




childAntibModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntibLrti90); 1-summary(childAntibModLrti90)[[9]][1,c(1,4,3)]

childPenModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datPenLrti90); 1-summary(childPenModLrti90)[[9]][1,c(1,4,3)]
childCephModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datCephLrti90); 1-summary(childCephModLrti90)[[9]][1,c(1,4,3)]
childMacroModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datMacroLrti90); 1-summary(childMacroModLrti90)[[9]][1,c(1,4,3)]
childAminoModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAminoLrti90); 1-summary(childAminoModLrti90)[[9]][1,c(1,4,3)]
childOthAntibModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthAntibLrti90); 1-summary(childOthAntibModLrti90)[[9]][1,c(1,4,3)]

childOthModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthLrti90); 1-summary(childOthModLrti90)[[9]][1,c(1,4,3)]


childAntibModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntibAom90); 1-summary(childAntibModAom90)[[9]][1,c(1,4,3)]

childPenModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datPenAom90); 1-summary(childPenModAom90)[[9]][1,c(1,4,3)]
childCephModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datCephAom90); 1-summary(childCephModAom90)[[9]][1,c(1,4,3)]
childMacroModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datMacroAom90); 1-summary(childMacroModAom90)[[9]][1,c(1,4,3)]
childAminoModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAminoAom90); 1-summary(childAminoModAom90)[[9]][1,c(1,4,3)]
childOthAntibModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthAntibAom90); 1-summary(childOthAntibModAom90)[[9]][1,c(1,4,3)]

childOthModAom90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthAom90); 1-summary(childOthModAom90)[[9]][1,c(1,4,3)]




childAntibModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntibTot); 1-summary(childAntibModTot)[[9]][1,c(1,4,3)]

childPenModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datPenTot); 1-summary(childPenModTot)[[9]][1,c(1,4,3)]
childCephModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datCephTot); 1-summary(childCephModTot)[[9]][1,c(1,4,3)]
childMacroModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datMacroTot); 1-summary(childMacroModTot)[[9]][1,c(1,4,3)]
childAminoModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAminoTot); 1-summary(childAminoModTot)[[9]][1,c(1,4,3)]
childOthAntibModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthAntibTot); 1-summary(childOthAntibModTot)[[9]][1,c(1,4,3)]

childOthModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthTot); 1-summary(childOthModTot)[[9]][1,c(1,4,3)]


childAntibModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntibLrtiTot); 1-summary(childAntibModLrtiTot)[[9]][1,c(1,4,3)]

childPenModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datPenLrtiTot); 1-summary(childPenModLrtiTot)[[9]][1,c(1,4,3)]
childCephModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datCephLrtiTot); 1-summary(childCephModLrtiTot)[[9]][1,c(1,4,3)]
childMacroModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datMacroLrtiTot); 1-summary(childMacroModLrtiTot)[[9]][1,c(1,4,3)]
childAminoModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAminoLrtiTot); 1-summary(childAminoModLrtiTot)[[9]][1,c(1,4,3)]
childOthAntibModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthAntibLrtiTot); 1-summary(childOthAntibModLrtiTot)[[9]][1,c(1,4,3)]

childOthModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datOthLrtiTot); 1-summary(childOthModLrtiTot)[[9]][1,c(1,4,3)]







motherAntibModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAntibTot); 1-summary(motherAntibModTot)[[9]][1,c(1,4,3)]

motherPenModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datPenTot); 1-summary(motherPenModTot)[[9]][1,c(1,4,3)]
motherCephModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datCephTot); 1-summary(motherCephModTot)[[9]][1,c(1,4,3)]
motherMacroModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datMacroTot); 1-summary(motherMacroModTot)[[9]][1,c(1,4,3)]
motherAminoModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAminoTot); 1-summary(motherAminoModTot)[[9]][1,c(1,4,3)]
motherOthAntibModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datOthAntibTot); 1-summary(motherOthAntibModTot)[[9]][1,c(1,4,3)]

motherOthModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datOthTot); 1-summary(motherOthModTot)[[9]][1,c(1,4,3)]




motherAntibModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAntibLrtiTot); 1-summary(motherAntibModLrtiTot)[[9]][1,c(1,4,3)]

motherPenModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datPenLrtiTot); 1-summary(motherPenModLrtiTot)[[9]][1,c(1,4,3)]
motherCephModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datCephLrtiTot); 1-summary(motherCephModLrtiTot)[[9]][1,c(1,4,3)]
motherMacroModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datMacroLrtiTot); 1-summary(motherMacroModLrtiTot)[[9]][1,c(1,4,3)]
motherAminoModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAminoLrtiTot); 1-summary(motherAminoModLrtiTot)[[9]][1,c(1,4,3)]
motherOthAntibModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datOthAntibLrtiTot); 1-summary(motherOthAntibModLrtiTot)[[9]][1,c(1,4,3)]

motherOthModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datOthLrtiTot); 1-summary(motherOthModLrtiTot)[[9]][1,c(1,4,3)]


childAntiVMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiV90); 1-summary(childAntiVMod90)[[9]][1,c(1,4,3)]
childAntiFMod90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiF90); 1-summary(childAntiFMod90)[[9]][1,c(1,4,3)]
childAntiVModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiVLrti90); 1-summary(childAntiVModLrti90)[[9]][1,c(1,4,3)]
childAntiFModLrti90 = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiFLrti90); 1-summary(childAntiFModLrti90)[[9]][1,c(1,4,3)]

childAntiVModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiVTot); 1-summary(childAntiVModTot)[[9]][1,c(1,4,3)]
childAntiFModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiFTot); 1-summary(childAntiFModTot)[[9]][1,c(1,4,3)]
childAntiVModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiVLrtiTot); 1-summary(childAntiVModLrtiTot)[[9]][1,c(1,4,3)]
childAntiFModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=child,data=datAntiFLrtiTot); 1-summary(childAntiFModLrtiTot)[[9]][1,c(1,4,3)]

motherAntiVModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAntiVTot); 1-summary(motherAntiVModTot)[[9]][1,c(1,4,3)]
motherAntiFModTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAntiFTot); 1-summary(motherAntiFModTot)[[9]][1,c(1,4,3)]
motherAntiVModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAntiVLrtiTot); 1-summary(motherAntiVModLrtiTot)[[9]][1,c(1,4,3)]
motherAntiFModLrtiTot = coxph(Surv(start,end,outcome)~intendArm+frailty(id)+strata(country),subset=mom,data=datAntiFLrtiTot); 1-summary(motherAntiFModLrtiTot)[[9]][1,c(1,4,3)]


################################################################################
##### incidence rates ##########################################################
################################################################################

datObj = datAntiFLrtiTot
age = datObj$mom==1&datObj$lmic==1
yrs1 = sum((datObj$end-datObj$start)[datObj$intendArm=='RSV F Vaccine'&age],na.rm=T)/365.25
yrs0 = sum((datObj$end-datObj$start)[datObj$intendArm=='Placebo'&age],na.rm=T)/365.25
events1 = sum(datObj$outcome[datObj$intendArm=='RSV F Vaccine'&age],na.rm=T)
events0 = sum(datObj$outcome[datObj$intendArm=='Placebo'&age],na.rm=T)
c(round(100*events1/yrs1,1),sum(events1),round(100*events0/yrs0,1),sum(events0))
sum(yrs1); sum(yrs0)




