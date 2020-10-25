R_dichotomization=function(lith)
{
  
  lith_binary=foreach(k=1:nrow(lith),.combine=rbind) %dopar%
  {
    dichotomise=function(x)
    {
      uu=rep(0,length(names(table(lith[,4]))))
      uu[x]=1
      return(uu)
    }
    dichotomise(lith[k,4])
  }
  colnames(lith_binary)=1:length(names(table(lith[,4])))
  return(lith_binary)
}