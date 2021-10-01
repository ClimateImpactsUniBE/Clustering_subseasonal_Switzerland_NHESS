source('~/clustering_library.R')

##########################################################################################
# LOAD DATA
##########################################################################################

DF = read.table('/scratch3/atuel/RhiresD/Data/daily_pr_catchments_RhiresD_1961-2017.txt',header=T)
id_list = as.numeric(sapply(tail(colnames(DF),-1),function(x){return(substr(x,2,5))}))
DF$d = as.Date(DF$d)
numMths = aggregate(p~y+m,data.frame(p=rep(1,length(DF$d)),m=month(DF$d),y=year(DF$d)),mean)
numMths = aggregate(p~m,numMths,sum)$p

##########################################################################################
# 21-DAY CLUSTER ANALYSIS
##########################################################################################

L = 21
QQ = 0.99
mth_seq = list(c(12,1,2),3:5,6:8,9:11)
fracInClus = fracDays = clusLenMean = clusTotPrMean = clusContPr = array(0,c(ncol(DF)-1,4))
sindex = 1
for (season in c('DJF','MAM','JJA','SON')){
  
  mth_sea = mth_seq[[sindex]]
  data = DF
  data[!month(DF$d)%in%mth_sea,2:ncol(data)] = 0

  for (j in 2:ncol(DF)){
    prVec = data[,j]
    ww = 0*prVec
    for (mm in mth_sea){
      ww[month(DF$d)==mm] = 1*(prVec>=quantile(prVec[month(DF$d)==mm],QQ,na.rm=T))[month(DF$d)==mm]
    }
    res = myDecluster(prVec,ww,r=2)
    X = 0*prVec
    X[res[[1]]] = 1 # position of events
    eec = movingSum(X,n=L)
    pra = movingSum(prVec,n=L) # precipitation accumulation
    clus_loc = clus_length = clus_tot_pr = c()
    if (max(eec,na.rm=T)>1){
      while(max(eec,na.rm=T)>1){
        inds = which(eec==max(eec,na.rm=T))
        clus_loc = c(clus_loc,inds[which.max(pra[inds])])
        clus_tot_pr = c(clus_tot_pr,max(pra[inds]))
        clus_length = c(clus_length,max(eec,na.rm=T))
        X[inds[which.max(pra[inds])]:(inds[which.max(pra[inds])]+L-1)] = 0
        prVec[inds[which.max(pra[inds])]:(inds[which.max(pra[inds])]+L-1)] = 0
        eec = movingSum(X,n=L) # extreme event count
        pra = movingSum(prVec,n=L) # precipitation accumulation
      }
    }
    # Fraction of extreme events included in clusters
    fracInClus[j-1,sindex] = (length(res[[1]])-sum(X))/length(res[[1]])
    if (length(clus_loc)>0){
      # Cluster length
      clusLenMean[j-1,sindex] = mean(clus_length)
      # Cluster fraction
      fracDays[j-1,sindex] = length(clus_loc)*L/sum(month(DF$d)%in%mth_sea)
      # Cluster total precipitation
      clusTotPrMean[j-1,sindex] = mean(clus_tot_pr)
      # Cluster contribution to precipitation
      clusContPr[j-1,sindex] = sum(clus_tot_pr)/sum(data[,j])
    }
  }
  sindex = sindex+1
}

write.table(fracInClus,row.names=F,file='~/seasonal_frac_in_clus_L21_monthly.txt')
write.table(fracDays,row.names=F,file='~/seasonal_frac_days_L21_monthly.txt')
write.table(clusLenMean,row.names=F,file='~/seasonal_clus_len_mean_L21_monthl.txt')
write.table(clusTotPrMean,row.names=F,file='~/seasonal_clus_tot_pr_mean_L21_monthly.txt')
write.table(clusContPr,row.names=F,file='~/seasonal_clus_cont_pr_L21_monthly.txt')
write.table(clusContPr,row.names=F,file='~/seasonal_clus_cont_pr_L21_monthly.txt')
