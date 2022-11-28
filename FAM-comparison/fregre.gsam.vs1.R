
fregre.gsam.vs1=function(resp,ldata,basis.x,alpha=0.05,limdcor=0.01){
  n=length(ldata$df[[resp]])
  nvar=length(names(ldata$df))-1+length(ldata)-1 #Global number of variates
  nesc=length(names(ldata$df))-1 #Number of scalar variables
  namesc=names(ldata$df)[which(names(ldata$df)!=resp)]
  namfun=names(ldata)[which(names(ldata)!="df")]
  dcor=matrix(0,nrow=nvar,ncol=nvar)
  ind=numeric(nvar)
  colnames(dcor)=c(namesc,namfun)
  names(ind)=colnames(dcor)
  rownames(dcor)=1:nvar
  form=as.formula(paste0(resp,"~1"))
  
  j=1
  Mset=c()
  Sset=c(namesc,namfun)
  modelo=fregre.gsam(formula=form,data=ldata,basis.x=basis.x)
  res=ldata$df[,resp]-mean(ldata$df[,resp])
  modfinish=FALSE
  
  while(!modfinish){
    if (nesc>0) {
      for (i in 1:length(namesc)){
        if (namesc[i] %in% Sset) {
          if (is.factor(ldata$df[[namesc[i]]])) {
            y=model.matrix(as.formula(paste0("~-1+",namesc[i])),ldata$df)
          } else {
            y=ldata$df[[namesc[i]]]
          }
          tt=dcor.xy(res,y,n=length(res))
          dcor[j,i]=tt$estimate*(tt$p.value<alpha)
          #print(paste0("Variable:",names(ldata$df)[namesc[i]]))
        } else {dcor[j,i]=NA}
      }
    }
    #for (i in (nesc+1):nvar){
    if (length(namfun)>0){
      for (i in 1:length(namfun)){
        nf=i-nesc
        if (namfun[i] %in% Sset){
          #print(paste0("Variable:",names(ldata)[nf+1]))
          dd=dcor.xy(ldata[[namfun[i]]],res,n=n)
          dcor[j,nesc+i]=dd$estimate*(dd$p.value<0.05)
        } else {dcor[j,nesc+i]=NA}
      }
    }
    #    dcor[j,colnames(dcor) %in% Mset]=-dcor[j,colnames(dcor) %in% Mset]
    jj=which(dcor[j,]==max(dcor[j,],na.rm=TRUE))
    cat(paste0("Iter:",j," dcor(",colnames(dcor)[jj],")= ",round(dcor[j,jj],4),"\n"))
    cand=colnames(dcor)[jj]
    if (max(dcor[j,],na.rm=TRUE)>limdcor) {
      modant=modelo
      formant=form
      if (jj<=nesc){
        Mset=c(Mset,cand)
        Sset=setdiff(Sset,cand)
        if (is.factor(ldata$df[,namesc[jj]])){
          form=update.formula(form, paste0(".~.+",namesc[jj]))
        }else{
          form=update.formula(form, paste0(".~.+s(",namesc[jj],")"))
        }
      } else {
        jj2=jj-nesc
        Mset=c(Mset,cand)
        Sset=setdiff(Sset,cand)
        form=update.formula(form, paste0(".~.+s(",namfun[jj2],")"))
      }
      j=j+1
    } else {modfinish=TRUE}
    modelo=fregre.gsam(formula=form,data=ldata,basis.x=basis.x)
    aa=summary(modant)$sp.criterion
    ab=summary(modelo)$sp.criterion
    #    if (ab>aa){
    if ((ab/aa)>.99){
      Mset=setdiff(Mset,cand)
      form=formant
      modelo=modant	
      #	modfinish=TRUE
    } else {res=residuals(modelo)}
    print(Mset)
  }
  return(list(form=form, ldata=ldata, basis.x=if (length(namfun)>0) {basis.x}else{NULL}, model=modelo,dcor=dcor))
}
