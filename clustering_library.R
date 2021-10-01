require(ncdf4)
require(lubridate)
require(extRemes)

# Runs declustering

myDecluster = function(v,x,r){
  
  # v: data series
  # x: binary time series of event occurrences (raw)
  # r: run length for declustering
  y = x
  w = rle(x)
  # Set runs of zeros with length<r to 1
  final_pos = cumsum(w$lengths)
  for (uu in 1:(length(w$lengths)-1)){
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] = 1
    }
  }
  if (x[length(x)]==1){
    uu = w$lengths[length(w$lengths)]
    if ((w$lengths[uu]<r) & (w$values[uu]==0)){
      y[(final_pos[uu]-w$lengths[uu]+1):final_pos[uu]] = 1
    }
  }
  w = rle(y)
  # maximum element of each run of '1's
  final_pos = cumsum(w$lengths)
  out_max_idx = out_max = out_clust_len = out_clust_evts = out_clust_tot = 0*1:sum(w$values==1)
  index = 1
  for (i in 1:length(final_pos)){
    if (w$values[i]==1){
      vec = v[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      vec2 = x[(final_pos[i]-w$lengths[i]+1):final_pos[i]]
      out_max[index] = max(vec)
      out_max_idx[index] = which.max(vec)+final_pos[i]-w$lengths[i]
      out_clust_len[index] = length(vec)
      out_clust_evts[index] = sum(vec2)
      out_clust_tot[index] = sum(vec)
      index = index+1
    }
  }
  return(list(out_max_idx,out_max,out_clust_len,out_clust_tot))
}

# Ripley's K function

ripleyK = function(x,t){
  # x: binary time series of event occurrences
  # t: window width
  S = 0
  inds = which(x==1)
  inds = inds[inds>t & inds<(length(x)-t+1)]
  for (i in inds){
    S = S+sum(x[max((i-t),0):(i-1)])+sum(x[(i+1):min((i+t),length(x))])
  }
  return(S/sum(x[inds]))
}
