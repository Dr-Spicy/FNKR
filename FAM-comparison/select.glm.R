# funciones del Script
# dcor.y
# dist.list
# select.gls
# select.glm

################################################################################
dcor.y<-function(ldist,response,n,bcdcor=TRUE){
  lenldist<-length(ldist)
  namldist<-names(ldist)
  if (missing(n)) {
    if (is.fdata(ldist[[1]]) | is.matrix(ldist[[1]]|is.data.frame(ldist[[1]])) )
      n<-nrow(ldist[[1]])
    if (is.vector(ldist[[1]]))    n <- length(response)
  }
  
  if (missing(response)) {print("se calculan todas las distancia")
    dst<-diag(lenldist)
    for (i in 1:(lenldist-1)) {
      for (j in i:(lenldist)) {
        if (bcdcor)     dst[i,j]<-dst[j,i]<- bcdcor.dist(ldist[[i]],ldist[[j]],n=n)
        else            dst[i,j]<-dst[j,i]<-cor.fdist(ldist[[i]],ldist[[j]])
      }}
    
    rownames(dst)<-colnames(dst)<-namldist
  }
  else{                     #se calculan todas las distancias respecto la respuesta
    if (length(response)>1) stop("Incorrect response label")
    ii<-response==namldist
    if (sum(ii)==0) stop("Incorrect response label")
    iresponse<-which(ii)
    dst<-numeric(lenldist-1)
    predictors<-setdiff(namldist,response)
    
    for (i in 1:(lenldist-1)) {
      #           dst[i]<-cor.fdist(ldist[[response]],ldist[[predictors[i]]])
      if (bcdcor)  dst[i]<-bcdcor.dist(ldist[[response]],ldist[[predictors[i]]],n=n)
      else dst[i]<-dcor.dist(ldist[[response]],ldist[[predictors[i]]])
    }
    names(dst)<-predictors
  }
  dst
}
################################################################################
################################################################################
#se calculan todas las distancias respecto as? mismas
dist.list<-function(ldata,...){
  lenldata<-length(ldata)
  ldist<-list()
  for (i in 1:(lenldata)) {
    if (is.factor(ldata[[i]])) ldata[[i]]<-model.matrix(~ldata[[i]]) #transforma el factor en dummyies
    
    if (is.fdata(ldata[[i]]))      ldist[[i]]<-metric.lp(ldata[[i]],...)
    else   ldist[[i]]<-as.matrix(dist(ldata[[i]]), diag =TRUE, upper = TRUE,     p = 2)
  }
  names(ldist)<-names(ldata)
  ldist
}
################################################################################
select.gls<- 
  function(data,y,x,pvalue,dcor,type.basis="pc",
           par.basis=list(),tol=.05,par.model,xydist){
    #0-
    if (missing(y)) {stop("Hay que especificar el nombre de la variable respuesta")  }
    if (missing(pvalue)) pvalue<-.05
    if (missing(par.model)) par.model<-list()   
    n.edf.correction <- FALSE
    namdata<-names(data)
    idf<-which("df"==namdata )
    ydat<-data$df[,y]
    namfunc<-names(data[-idf])
    namnfunc<-names(data$df)
    #namnfunc<-setdiff(namnfunc,y)
    if (missing(x)) {  
      #    ldata<-as.list(data[idf][,y,drop=F])
      xdatos<-data[-idf]    
      x<-names(data[-idf])
      warning("Solo las variables funcionales son consideradas como covariables")
    }
    else {
      ifunc<-intersect(x,namfunc)
      infunc<-intersect(c(y,x),namnfunc)
      xdatos<-as.list(data$df[,infunc,drop=F])
      xdatos<-c(xdatos,data[ifunc])
    } 
    ldata0<-xdatos
    
    #  ldata0<-c(xdatos,data$df[,y,drop=F])
    #xdatos<-c(xdatos,data[namfunc])
    #    print(names(xdatos))
    #    print(names(ldata0))
    resp<-y
    
    xentran<-NULL
    tipoentran<-NULL
    #1- c?clulo distancias de cada objeto consigo mismo
    #  print("1 calculando distancias")
    xynam<-names(ldata0)
    if (missing(xydist))     xydist<-ldist0<-dist.list(ldata0)
    else ldist0<-xydist[xynam]
    # dist_ij<-dcor.y(ldist0) 
    parar<-FALSE
    it<-1
    if (is.null(par.basis$kmax) & type.basis=="pc") par.basis$kmax<-kmax<-8
    basisx<-NULL
    npredictors<-length(ldata0)-1
    ipredictors<-numeric(npredictors)
    names(ipredictors)<- setdiff(names(ldata0),resp)
    #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
    # print("2 calculando correlaciones")
    dist_resp<-dcor.y(ldist0,resp)
    #dist_resp
    n.edf<-length(ydat)
    fpredictors.nl<-""
    form.nl<-paste(resp,"~",sep="")
    basis2<-list()  
    ycen<-data$df[,y]-mean(data$df[,y])
    gof<-NULL
    anyfdata<-FALSE
    while (!parar){
      # print("3 Seleccion variable-Regresion")
      #      print("bucle")    
      esfactor<-FALSE #SE UTILIZA PARA NO PONER s(factor)
      
      #3- selecci?n dcor m?s elevada
      nam<-  names(dist_resp)
      ind.xentra<-which.max(dist_resp)
      xentra<-nam[ind.xentra]
      #  dd<-dcor.ttest(ldist0[[resp]],ldist0[[xentra]],distance=TRUE)
      
      dd<-dcor.test(ldist0[[resp]],ldist0[[xentra]],R=R)
      if (par.fda.usc$verbose) {
        #       print(it)
        #        print(dist_resp)
        #	 print(n.edf)
        #      print(nam)
        cat("Entra la variable ",xentra)
        #        print(dd)
        #print(names(dd))
      }
      
      #      if (dd$p.value>tol) {parar=TRUE
      #print(dd$estimate)
      
      if (dd$p.value>tol) {parar=TRUE
      #print("para pq ninguna variable es significativa")
      }
      else{
        #print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        #print(xentra)
        if (is.fdata(ldata0[[xentra]])) {
          anyfdata<-TRUE
          par.basis$fdataobj<-ldata0[[xentra]]
          #entran las Componentes significativas en la regressi?n lineal #si se hace con Gam son muchas combinaciones    
          nam<-paste("fregre.",type.basis,".cv",sep="")
          par.basis$y<-ldata0[[resp]]
          # if par.basis$kmax ya kmax o basis.x ya vendran datdos
          #print(nam)
          #print(names(par.basis))
          #print(names(ldata0))
          res<- do.call(nam,par.basis)
          #print(names(res))
          
          #      res2<-fregre.basis.cv(ldata0[[xentra]],ldata0[[resp]],kmax)
          #          res<-fregre.pc.cv(ldata0[[xentra]],ldata0[[resp]],kmax)
          switch(type.basis,"pc"={
            best.pc<-res$pc.opt
            basis1<-create.pc.basis(ldata0[[xentra]],best.pc)
          },
          "pls"={
            best.pc<-res$pls.opt
            basis1<-create.pls.basis(ldata0[[xentra]],ldata0[[resp]],best.pc)
          },"basis"={
            best.pc<-res$basis.x.opt$nbasis     
            basis1<-res$basis.x.opt
            basis2<-res$basis.b.opt          
          })
          if (par.fda.usc$verbose) {
            print(paste("Selecci?n de las bases por CV: regresion(e~",xentra,")",sep=""))
            cat("Indice Bases ?ptimas: ",best.pc,"\n")          
          }         
          
          
          if (is.null(basisx)) { 
            basisb<-basisx<-list()
          }
          basisx[[xentra]]<-basis1
          basisb[[xentra]]<-basis2
          #    print(names(basisx))
        }
        #     print(xentra)
        #     print(class(ldata0[[xentra]]))          
        if (!parar) {
          xentran<-c(xentran,xentra)
          ind.xentra2<-which(xentra==names(ldist0))
          ldist0<-ldist0[-ind.xentra2]
          ipredictors[xentra]<-ipredictors[xentra]+1
          #4- contrucci?n del modelo para esta variable #consido un cat?logo de 4 posibilidades (lineal/nolineal, funcional/scalar)
          #
          #  form.l<-as.formula(paste(resp,"~",fpredictors.l,sep="",collapse="+"))
          # el modelo funcional fregre tiene como caso particular el glm o gam (ver para factores)/iteracciones/etc
          #print(linear)
          fpredictors.nl<-paste(xentran,sep="",collapse="+")
          form.nl<-as.formula(paste(resp,"~",fpredictors.nl,sep="",collapse="+"))
          #          if (type.basis=="basis")         res.nl<-fregre.glm(form.nl,data=data,basis.x=basisx,basis.b=basisb)
          #          else        res.nl<-fregre.glm(form.nl,data=data,basis.x=basisx)
          # sino hay correlation es un lm
          if (!anyfdata){
            #	    par.model2<-par.model	
            nam.model<-"gls"           #estaba glse
            par.model$model<-form.nl
            par.model$data<-data$df
            if (is.null(par.model$correlation)) par.model[["correlation"]]<-NULL	
            #print(names(par.model))
            res.nl<-do.call(nam.model,par.model)
            #print(names(par.model))
            
          }
          else{
            #print("1")            
            names(par.model)[which(names(par.model)=="model")]<-"formula"
            nam.model<-"fregre.gls"          
            par.model$formula<-form.nl
            par.model$data<-data
            par.model[["basis.x"]]<-basisx
            if (is.null(par.model$correlation)) par.model[["correlation"]]<-NULL
            if (type.basis=="basis")  par.model$basis.b<-basisb 
            #print(nam.model)
            #print(basisx)
            #print(names(par.model))
            #print(form.nl)
            res.nl<-do.call(nam.model,par.model) 
            #print("2")          
            
          }
          # ajuste del modelo GLS
          #print(nam.model)
          #print(names(par.model))
          
          suma<-summary(res.nl)          
          dd <- res.nl$dims
          df <- dd[["p"]]
          edf<-dd[["N"]] -         dd[["p"]]	
          #         sr2 <- sum(res.nl$residuals^2)/edf
          r2 <- 1 - sum(res.nl$residuals^2)/sum(ycen^2)
          if (n.edf.correction) n.edf<- edf
          ldata0[[resp]]<-res.nl$residuals
          ldist0[[resp]]<-as.matrix(dist(ldata0[[resp]]), diag =TRUE, upper = TRUE,     p = 2)
          #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
          if (it==npredictors) {
            parar=TRUE
            gof<-rbind(gof,c(suma$logLik,suma$BIC,suma$AIC,edf,r2))
            #else gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre))
          }
          else{
            it<-it+1
            dist_resp<-dcor.y(ldist0,resp)
          }
          #     print(summary(res.nl))
          if (!parar){ 
            gof<-rbind(gof,drop(c(suma$logLik,suma$BIC,suma$AIC,edf,r2)))
          }
        }  }
      
    }  
    # print(dim(gof));print(gof)
    gof<-data.frame(xentran,(gof))
    #  print(gof)
    #  print("aaa ver ohhh")
    #gof<-as.data.frame(gof)
    if (is.null(xentran)) {
      warning("ninguna variable seleccionada, se estima un LM a pelo")    
      res.nl<-fregre.glm(as.formula(paste(resp,"~1",sep="")),data=ldata)
      suma<-summary.lm(res.nl)
      gof<-data.frame(rbind(c(1,AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.squared)))
    }  
    # print(gof)
    # print(dim(gof))
    names(gof)<-c("xentra","logLik","BIC","AIC","edf","r.sq")
    #else names(gof)<-c("xentra","AIC","deviance","df.residual","r.sq","dev.expl","GCV.ubre")
    out<-list("model"=res.nl,"gof"=gof,"i.predictor"=ipredictors,"xydist"=xydist)
    return(out)
  }
################################################################
#vsgls<-select.gls(lftrain,y="ly0t1",x=xvar1,par.model=par.cor,tol=dcor.tol,xydist=vsglm$xydist ,par.basis=par.base)      

# hacer un par.control q vayan las bases
select.glm<-
  function(data,y,x,pvalue,dcor,type.basis="pc",par.basis=list(),tol=.05,xydist,R=100){
    #0-
    if (missing(y)) {stop("Hay que especificar el nombre de la variable respuesta")  }
    if (missing(pvalue)) pvalue<-.05
    namdata<-names(data)
    idf<-which("df"==namdata )
    ydat<-data$df[,y]
    namfunc<-names(data[-idf])
    namnfunc<-names(data$df)
    #namnfunc<-setdiff(namnfunc,y)
    if (missing(x)) {  
      #    ldata<-as.list(data[idf][,y,drop=F])
      xdatos<-data[-idf]    
      x<-names(data[-idf])
      warning("Solo las variables funcionales son consideradas como covariables")
    }
    else {
      ifunc<-intersect(x,namfunc)
      infunc<-intersect(c(y,x),namnfunc)
      xdatos<-as.list(data$df[,infunc,drop=F])
      xdatos<-c(xdatos,data[ifunc])
    } 
    ldata0<-xdatos
    #  ldata0<-c(xdatos,data$df[,y,drop=F])
    #xdatos<-c(xdatos,data[namfunc])
    resp<-y
    
    xentran<-NULL
    tipoentran<-NULL
    #1- c?clulo distancias de cada objeto consigo mismo
    #  print("1 calculando distancias")
    #  ldist0<-dist.list(ldata0)
    xynam<-names(ldata0)
    if (missing(xydist))     xydist<-ldist0<-dist.list(ldata0)
    else ldist0<-xydist[xynam]
    
    # dist_ij<-dcor.y(ldist0) 
    parar<-FALSE
    it<-1
    if (is.null(par.basis$kmax) & type.basis=="pc") par.basis$kmax<-kmax<-8
    basisx<-NULL
    npredictors<-length(ldata0)-1
    ipredictors<-numeric(npredictors)
    names(ipredictors)<- setdiff(names(ldata0),resp)
    #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
    #  print("2 calculando correlaciones")
    dist_resp<-dcor.y(ldist0,resp)
    #dist_resp
    n.edf<-length(ydat)
    fpredictors.nl<-""
    form.nl<-paste(resp,"~",sep="")
    basis2<-list()  
    gof<-NULL
    while (!parar){
      #print("VS-Regresion GLM")    
      esfactor<-FALSE #SE UTILIZA PARA NO PONER s(factor)
      
      #3- selecci?n dcor m?s elevada
      nam<-  names(dist_resp)
      ind.xentra<-which.max(dist_resp)
      xentra<-nam[ind.xentra]
      #dd<-dcor.ttest(ldist0[[resp]],ldist0[[xentra]],distance=TRUE)
      #dd<-dcor.test(ldist0[[resp]],ldist0[[xentra]],R=R)
      dd<-fda.usc:::dcor.test(ldist0[[resp]],ldist0[[xentra]],n=n.edf)
      
      if (par.fda.usc$verbose) {
        print(it)
        print(dist_resp)
        #      print(nam)
        cat("Entra la variable ",xentra)
        print(dd)
      }
      if (dd$p.value>tol) {parar=TRUE
      #	if (dd$estimate<tol) {parar=TRUE
      #print("para pq ninguna variable es significativa")
      }
      else{
        if (is.fdata(ldata0[[xentra]])) {
          par.basis$fdataobj<-ldata0[[xentra]]
          #entran las Componentes significativas en la regressi?n lineal #si se hace con Gam son muchas combinaciones    
          nam<-paste("fregre.",type.basis,".cv",sep="")
          par.basis$y<-ldata0[[resp]]
          # if par.basis$kmax ya kmax o basis.x ya vendran datdos
          res<- do.call(nam,par.basis)
          #      res2<-fregre.basis.cv(ldata0[[xentra]],ldata0[[resp]],kmax)
          #          res<-fregre.pc.cv(ldata0[[xentra]],ldata0[[resp]],kmax)
          switch(type.basis,"pc"={
            best.pc<-res$pc.opt
            basis1<-create.pc.basis(ldata0[[xentra]],best.pc)
          },
          "pls"={
            best.pc<-res$pls.opt
            basis1<-create.pls.basis(ldata0[[xentra]],ldata0[[resp]],best.pc)
          },"basis"={
            best.pc<-res$basis.x.opt$nbasis     
            basis1<-res$basis.x.opt
            basis2<-res$basis.b.opt          
          })
          if (par.fda.usc$verbose) {
            print(paste("Selecci?n de las bases por CV: regresion(e~",xentra,")",sep=""))
            cat("Indice Bases ?ptimas: ",best.pc,"\n")          
          }       
          if (is.null(basisx)) { 
            basisb<-basisx<-list()
          }
          basisx[[xentra]]<-basis1
          basisb[[xentra]]<-basis2
          #    print(names(basisx))
        }
        #     print(xentra)
        #     print(class(ldata0[[xentra]]))          
        if (!parar) {
          xentran<-c(xentran,xentra)
          ind.xentra2<-which(xentra==names(ldist0))
          ldist0<-ldist0[-ind.xentra2]
          ipredictors[xentra]<-ipredictors[xentra]+1
          #4- contrucci?n del modelo para esta variable #consido un cat?logo de 4 posibilidades (lineal/nolineal, funcional/scalar)
          #
          #  form.l<-as.formula(paste(resp,"~",fpredictors.l,sep="",collapse="+"))
          # el modelo funcional fregre tiene como caso particular el glm o gam (ver para factores)/iteracciones/etc
          #print(linear)
          fpredictors.nl<-paste(xentran,sep="",collapse="+")
          form.nl<-as.formula(paste(resp,"~",fpredictors.nl,sep="",collapse="+"))
          if (type.basis=="basis")         res.nl<-fregre.glm(form.nl,data=data,basis.x=basisx,basis.b=basisb)
          else        res.nl<-fregre.glm(form.nl,data=data,basis.x=basisx)
          #    print(summary.lm(res.nl))
          suma<-summary.lm(res.nl)
          #  if (n.edf.correction) n.edf<-res.nl$df.residual
          ldata0[[resp]]<-res.nl$residuals
          ldist0[[resp]]<-as.matrix(dist(ldata0[[resp]]), diag =TRUE, upper = TRUE,     p = 2)
          #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
          #   res.final<-fregre.gsam(form.nl,data=ldata)
          if (it==npredictors) {parar=TRUE
          gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.squared))
          #                              else gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre))
          }
          else{
            it<-it+1
            dist_resp<-dcor.y(ldist0,resp)
          }
          #     print(summary(res.nl))
          if (!parar){ 
            gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.squared))
          }
        }  }
      
    }  
    # print(dim(gof));print(gof)
    gof<-data.frame(xentran,(gof))
    #  print(gof)
    #  print("aaa ver ohhh")
    #gof<-as.data.frame(gof)
    if (is.null(xentran)) {
      warning("ninguna variable seleccionada")    
      res.nl<-fregre.glm(as.formula(paste(resp,"~1",sep="")),data=ldata0)
      suma<-summary.lm(res.nl)
      gof<-data.frame(rbind(c(1,AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.squared)))
    }  
    # print(gof)
    # print(dim(gof))
    names(gof)<-c("xentra","AIC","deviance","edf","r.sq")
    #else names(gof)<-c("xentra","AIC","deviance","df.residual","r.sq","dev.expl","GCV.ubre")
    out<-list("model"=res.nl,"gof"=gof,"i.predictor"=ipredictors,"xydist"=xydist)
    return(out)
  }
################################################################

print.select.glm <- function(x,...){
  cat("\n Model selected:\n")
  print(x[[1]])
  cat("\n Stepwise GoF:\n")
  print(x[[2]])
  cat("\n Predictor selected:\n")
  print(x[[3]])  
}


print.select.gsam <- function(x,...){
  cat("\n Model selected:\n")
  print(x[[1]])
  cat("\n Stepwise GoF:\n")
  print(x[[2]])
  cat("\n Predictor selected:\n")
  print(x[[3]])  
}

print.select.bgsam <- function(x,...){
  #cat("\n Model selected:\n")
  #print(x[[1]])
  #cat("\n Stepwise GoF:\n")
  #print(x[[2]])
  cat("\n Predictor selected (by order):\n")
  print(x[[4]])
  cat("\n Model type selected (by order):\n")
  summary(m3[[1]])
}


################################################################################
select.gsam.cv<-function(dat,y,x,pvalue,dcor,type.basis="pc",criterio="AICc",kmax=8,
                         par.basis=list(),tol=.1,par.model,kbs,xydist,trace=TRUE){
  #0-
  if (missing(y)) {stop("Hay que especificar el nombre de la variable respuesta")  }
  if (missing(pvalue)) pvalue<-.05
  if (missing(par.model)) par.model<-list() 
  if (missing(kbs)) kbs<--1
  n.edf.correction<-FALSE
  namdata<-names(dat)
  idf<-which(namdata=="df")
  ydat<-dat$df[,y]
  namfunc<-names(dat[-idf])
  namnfunc<-names(dat$df)
  # print("***************     nombre variables *************")  
  #print(namfunc)
  # print(namnfunc)
  #namnfunc<-setdiff(namnfunc,y)
  if (missing(x)) {  
    #    ldata<-as.list(data[idf][,y,drop=F])
    xdatos<-dat[-idf]    
    x<-names(dat[-idf])
    warning("Solo las variables funcionales son consideradas como covariables")
  }
  else {
    ifunc<-intersect(x,namfunc)
    infunc<-intersect(c(y,x),namnfunc)
    xdatos<-as.list(dat$df[,infunc,drop=F])
    xdatos<-c(xdatos,dat[ifunc])
  } 
  
  ldata0<-xdatos
  #  ldata0<-c(xdatos,data$df[,y,drop=F])
  #xdatos<-c(xdatos,data[namfunc])
  #    print(names(xdatos))
  # print(names(ldata0))
  resp<-y
  
  xentran<-NULL
  tipoentran<-NULL
  #1- c?clulo distancias de cada objeto consigo mismo
  #print("1 calculando distancias")
  #    ldist0<-dist.list(ldata0)
  xynam<-names(ldata0)
  #print(xynam)
  
  if (missing(xydist))    {
    xydist<-ldist0<-dist.list(ldata0)
  }    else {
    ldist0<-xydist[xynam]
  }
  #  print(names(ldist0))
  # print("fin names(ldist0")
  # dist_ij<-dcor.y(ldist0) 
  parar<-FALSE
  it<-1
  basisx<-NULL
  npredictors<-length(ldata0)
  ipredictors<-numeric(npredictors-1)
  names(ipredictors)<- setdiff(names(ldata0),resp)
  # print("32222")
  # print(ipredictors)  
  
  #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
  # print("2 calculando correlaciones")
  dist_resp<-dcor.y(ldist0,resp)
  #dist_resp
  n.edf<-length(ydat)
  fpredictors.nl<-""
  form.nl<-paste(resp,"~",sep="")
  basis2<-list()  
  ycen<-dat$df[,y]-mean(dat$df[,y])
  gof<-NULL
  anyfdata<-FALSE
  
  while (!parar){
    
    #  print("3 Seleccion variable-Regresion")
    #      print("bucle")    
    esfactor<-FALSE #SE UTILIZA PARA NO PONER s(factor)
    
    #3- selecci?n dcor m?s elevada
    nam<-  names(dist_resp)
    ind.xentra<-which.max(dist_resp)
    xentra<-nam[ind.xentra]
    #dd<-dcor.ttest(ldist0[[resp]],ldist0[[xentra]],distance=TRUE)
    #print(n.edf)    
    #print(dd)
    dd<-fda.usc:::dcor.test(ldist0[[resp]],ldist0[[xentra]],n=n.edf)
    if (trace){
      print(dd)
      print(xentra)
    }   #print(dd2)
    #print("33333")
    #   dd=dcor.xy(ldata0[[resp]],ldata0[[xentra]],n=n.edf)
    #print(dd)
    #  dcor[j,i]=dd$estimate*(dd$p.value<0.05)
    #    print(names(ldist0))
    if (par.fda.usc$verbose) {
      print(it)
      print(dist_resp)
      print(n.edf)
      print(nam)
      cat("Entra la variable ",xentra)
      print(dd)
      print(names(dd))
    }
    #  print(dd)
    #      if (dd$p.value>tol) {parar=TRUE
    if (dd$estimate<tol) {
      parar=TRUE
      if (trace) print("para pq ninguna variable es significativa")
    }
    else{
      if (trace)      print("aaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
      if (trace)       print(xentra)
      
      if (is.fdata(ldata0[[xentra]])) {
        # print(xentra)
        if (xentra %in% names(par.basis)){
          itype<-which(xentra== names(par.basis))
          par.basis.ipred<-par.basis[[itype]]
        }    else {par.basis.ipred <- list()}
        anyfdata<-TRUE
        par.basis.ipred$fdataobj<-ldata0[[xentra]]
        
        #if (is.list(type.basis))
        if (xentra %in% names(type.basis)){
          itype<-which(xentra== names(type.basis))
          type.basis.ipred<-type.basis[itype]
        }    else {type.basis.ipred <- type.basis[1]}
        
        
        tbasis<-type.basis.ipred
        if (type.basis.ipred=="fourier") tbasis<-"basis"
        if (type.basis.ipred=="bspline") tbasis<-"basis"
        
        #################
        #nam<-paste("fregre.",tbasis,".cv",sep="")
        #par.basis.ipred$y<-ldata0[[resp]]
        # print("aaaaaaa nnnnnnnn tttttttt eeeeeeesssssssssss ggggggsssssssaaaaaaaammmmmmmm")        
        nam<-"fregre.gsam.cv"
        #  par.basis.ipred<-list()
        par.basis.ipred$y<-resp#entra la etiqueta y no toda la variable como en el basis.cv o pc.cv
        par.basis.ipred$x<-xentra
        par.basis.ipred$dat<-dat
        #      print(names(data))
        #      print(resp)
        #      print(xentra)
        # if par.basis$kmax ya kmax o basis.x ya vendran datdos
        #          print(nam)
        #print(names(par.basis))
        #print(names(ldata0))
        #   print("antes res")
        #res<- do.call(nam,par.basis.ipred)
        res<- fregre.gsam.cv(dat,resp,xentra,criterio=criterio,kmax=kmax)
        if (trace)   print("res fregre.gsam.cv")
        if (trace)   print(summary(res))
        
        #newww
        if (is.null(res$pc.opt)) {parar=TRUE}
        else{
          switch(tbasis,"pc"={
            best.pc<-res$pc.opt
            #############if (is.null(res$pc.opt)) {}
            basis1<-create.pc.basis(ldata0[[xentra]],best.pc) #es lo mismo ?
            ###          basis1<-z$basis.x[["xentra"]]
          },
          "pls"={
            best.pc<-res$pls.opt
            basis1<-create.pls.basis(ldata0[[xentra]],ldata0[[resp]],best.pc)
          },"basis"={
            best.pc<-res$basis.x.opt$nbasis     
            basis1<-res$basis.x.opt
            basis2<-res$basis.b.opt          
          })
          
          if (par.fda.usc$verbose) {
            print(paste("Selecci?n de las bases por CV: regresion(e~",xentra,")",sep=""))
            cat("Indice Bases optimas: ",best.pc,"\n")          
          }         
          
          if (is.null(basisx)) { 
            basisb<-basisx<-list()
          }
          basisx[[xentra]]<-basis1
          basisb[[xentra]]<-basis2
          #    print(names(basisx))
        }
      }
      #     print(xentra)
      #     print(class(ldata0[[xentra]]))          
      if (!parar) {
        xentran<-c(xentran,xentra)
        ind.xentra2<-which(xentra==names(ldist0))
        ldist0<-ldist0[-ind.xentra2]
        ipredictors[xentra]<-ipredictors[xentra]+1
        #4- contrucci?n del modelo para esta variable #consido un cat?logo de 4 posibilidades (lineal/nolineal, funcional/scalar)
        #
        #  form.l<-as.formula(paste(resp,"~",fpredictors.l,sep="",collapse="+"))
        # el modelo funcional fregre tiene como caso particular el glm o gam (ver para factores)/iteracciones/etc
        #print(linear)
        #          fpredictors.nl<-paste(xentran,sep="",collapse="+")
        #          form.nl<-as.formula(paste(resp,"~",fpredictors.nl,sep="",collapse="+"))
        if (is.factor(ldata0[[xentra]])) esfactor=TRUE                 
        # print("despues")
        if (!esfactor)      {
          #  print("no es factor")
          fpredictors.nl<-paste("s(",xentra,",k=",kbs,")",sep="",collapse="+")
          form.nl
          
        }
        else         {    #hay factor el algun momento
          #  print(" es factor")
          #  print(fpredictors.nl)
          fpredictors.nl<-xentra             
        }
        form.nl<-(paste(form.nl,"+",fpredictors.nl,sep=""))                               
        #          if (type.basis=="basis")         res.nl<-fregre.glm(form.nl,data=data,basis.x=basisx,basis.b=basisb)
        #          else        res.nl<-fregre.glm(form.nl,data=data,basis.x=basisx)
        # sino hay correlation es un lm
        #      if (!anyfdata){
        #            nam.model<-"gam"          
        #        par.model$formula<-form.nl
        #        par.model$data<-data$df
        #res.nl<-do.call(nam.model,par.model)            
        #          }
        #          else{
        #names(par.model)[which(names(par.model)=="model")]<-"formula"
        nam.model<-"fregre.gsam"          
        par.model$formula<-as.formula(form.nl)
        par.model$dat<-dat
        par.model[["basis.x"]]<-basisx
        ##            if (tb=="basis")  par.model$basis.b<-basisb 
        res.nl<-do.call(nam.model,par.model)             
        if (trace)         print("sale fregre.gsam") #         }
        if (trace) print(" ajuste del modelo GLS")
        #           print(nam.model)
        # print(names(par.model))          
        suma<-summary(res.nl)          
        dd <- res.nl$dims
        df <- dd[["p"]]
        edf<-dd[["N"]] -         dd[["p"]]    
        #         sr2 <- sum(res.nl$residuals^2)/edf
        r2 <- 1 - sum(res.nl$residuals^2)/sum(ycen^2)
        if (n.edf.correction) n.edf<- edf
        ldata0[[resp]]<-res.nl$residuals
        ldist0[[resp]]<-as.matrix(dist(ldata0[[resp]]), diag =TRUE, upper = TRUE,     p = 2)
        #2- c?clulo correlaci?n  de cada la distancia de la respuesta vs distancia del resto de objetos
        if (it==npredictors) {
          parar=TRUE
          #gof<-rbindgof,c(suma$logLik,suma$BIC,suma$AIC,edf,r2))
          gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre))
        }
        else{
          it<-it+1
          dist_resp<-dcor.y(ldist0,resp)
        }
        #     print(summary(res.nl))
        if (!parar){ 
          #            gof<-rbind(gof,drop(c(suma$logLik,suma$BIC,suma$AIC,edf,r2)))
          gof<-rbind(gof,c(AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre))
        }
      }  }
    
  }  
  # print(dim(gof));print(gof)
  gof<-data.frame(xentran,(gof))
  #  print(gof)
  #  print("aaa ver ohhh")
  #gof<-as.data.frame(gof)
  if (is.null(xentran)) {
    warning("ninguna variable seleccionada, se estima un LM a pelo")    
    res.nl<-fregre.gsam(as.formula(paste(resp,"~1",sep="")),data=ldata)
    suma<-summary(res.nl)
    gof<-data.frame(rbind(c(1,AIC(res.nl),deviance(res.nl),res.nl$df.residual,suma$r.sq,suma$dev.expl,res.nl$gcv.ubre)))      
  }  
  # print(gof)
  # print(dim(gof))"
  names(gof)<-c("xentra","AIC","deviance","df.residual","r.sq","dev.expl","GCV.ubre")
  out<-list("model"=res.nl,"gof"=gof,"i.predictor"=ipredictors,"xydist"=xydist)
  return(out)
}
################################################################################


