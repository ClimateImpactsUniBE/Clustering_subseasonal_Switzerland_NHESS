
source('~/clustering_library.R')

########################################################################
## LOAD DATA ##
########################################################################

DF = read.table('~/Data/daily_pr_catchments_RhiresD_1961-2017.txt',header=T)
DF$d = as.Date(DF$d)
id_list = as.numeric(sapply(tail(colnames(DF),-1),function(x){return(substr(x,2,5))}))

#######################################################################
## SEASONAL RIPLEY'S K ##
#######################################################################

# Percentile threshold for extreme precipitation
QQ = 0.99
mth_seq = list(c(12,1,2),3:5,6:8,9:11)
matK = array(0,c(ncol(DF)-1,4,length(seq(2,50,2))))
matSumX = array(0,c(ncol(DF)-1,4))
for (season in c('DJF','MAM','JJA','SON')){
  
  sindex = which(c('DJF','MAM','JJA','SON')==season)
  mths = mth_seq[[sindex]]
  print(season)
  print(mths)
  indSeason = which(month(DF$d)%in%mths)
  
  for (j in 2:ncol(DF)){
    
    prVec = DF[,j]
    prVec[-indSeason] = 0
    ww = 0*prVec
    for (mm in mths){
      ww[month(DF$d)==mm] = 1*(prVec>=quantile(prVec[month(DF$d)==mm],QQ,na.rm=T))[month(DF$d)==mm]
    }
    res = myDecluster(prVec,ww,r=2)
    X = 0*ww
    X[res[[1]]] = 1
    matK[j-1,sindex,] = sapply(X=seq(2,50,2),FUN=ripleyK,x=X)
    matSumX[j-1,sindex] = sum(X)
  }
}

#######################################################################
## FDR SIGNIFICANCE ##
#######################################################################

kVec = 5:60
qVec = c(seq(0.005,0.95,0.005),seq(0.951,0.999,0.001))

for (season in c('DJF','MAM','JJA','SON')){
  
  load(paste('~/Data/mat_all_simulated_percentiles_RhiresD_kVec_1_60_',season,'.RData',sep=''))
  s = which(c('DJF','MAM','JJA','SON')==season)
  
  matQ = 0*matK[,s,]
  for (i in 1:nrow(matQ)){
    for (tt in 2:dim(matQ)[2]){
      if (!is.na(matK[i,s,tt]) & matK[i,s,tt]>0 & matSumX[i,s]%in%kVec){
        qq = matSimu[kVec==matSumX[i,s],,tt]
        if (max(qq)<matK[i,s,tt]){
          matQ[i,tt] = 1
        }
        if (min(qq)<matK[i,s,tt] & max(qq)>matK[i,s,tt]){
          matQ[i,tt] = qVec[tail(which(qq<matK[i,s,tt]),1)]
        }
      }
    }
  }
  
  alpha = 0.05
  p_star = 1:dim(matQ)[2]
  for (tt in 2:dim(matQ)[2]){
    # Sort p-values
    sorted_p = 1-sort(c(matQ[,tt]))
    # Updated p-value
    res = mapply(function(x,y){return(x<=y)},sorted_p,alpha*(1:length(sorted_p))/length(sorted_p))
    p_star[tt] = sorted_p[which(res==TRUE)[1]]
  }
  
  A = 0*matQ
  for (tt in 2:ncol(matQ)){
    A[,tt] = 1*((1-matQ)[,tt]<=p_star[tt])
  }
  
  # Extract information about significance
  A1 = 1*apply(A[,seq(2,50,2)%in%6:14],1,function(x){return(sum(x>0)>2)})
  A2 = 1*apply(A[,seq(2,50,2)%in%16:24],1,function(x){return(sum(x>0)>2)})
  A3 = 1*apply(A[,seq(2,50,2)%in%26:34],1,function(x){return(sum(x>0)>2)})
  A4 = 1*apply(A[,seq(2,50,2)%in%36:44],1,function(x){return(sum(x>0)>2)})
  
  # write.table(data.frame(x1=A1,x2=A2,x3=A3,x4=A4),row.names=F,file=paste('~/matClus_FDR_0.05_',season,'_Q99.txt',sep=''))
}
