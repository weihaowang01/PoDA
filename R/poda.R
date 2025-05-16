#the imputation method proposed by Zhou et al.


.imputation=function(Y,Z,formulaa){
  n=ncol(Y)
  m=nrow(Y)
  adaptive=TRUE
  corr.cut = 0.1
  allvars <- all.vars(as.formula(formulaa))
  names(Z) <- allvars
  Y<-as.data.frame(Y)
  Z<-as.data.frame(Z)
  if (any(Y == 0)) {
    N <- colSums(Y)
    
    if (adaptive) {
      logN <- log(N)
      tmp <- lm(as.formula(paste0("logN", formulaa)),
                Z)
      
      corr.pval <- coef(summary(tmp))[-1, "Pr(>|t|)"]
      if (any(corr.pval <= corr.cut)) {
        imputation <- TRUE
      }else {
        imputation <- FALSE
      }
    }
    if (imputation) {
      N.mat <- matrix(rep(N, m), nrow = m, byrow = TRUE)
      N.mat[Y > 0] <- 0
      tmp <- N[max.col(N.mat)]
      Y <- Y + N.mat/tmp
    }else {
      Y <- Y + 0.5
    }
  }
  
  return(as.matrix(Y))
}





# Hard threshold function
.hardPODA = function(X, Y,standarad_sigma,lmd){
  r = Y
  N = length(Y)
  gamma = matrix(rep(0, N), N, 1)
  theta = lmd/standarad_sigma
  gamma[abs(r) > theta] = r[abs(r) > theta]
  return(list(gamma=gamma, ress = r))
  
}






# Helper function for tuning parameter selection

.PODAFUNnew = function(y_poda,omega.i,standarad_sigma1,lmd){
  m=length(y_poda)
  gamma.j<-matrix(0,m,2)
  
  ######## the first iteration #########
  test_x_1<-.hardPODA(NULL,y_poda,standarad_sigma1,lmd) #hard
  mean.alpha.est<-median( -(omega.i[which(test_x_1$gamma==0)]))
  if(is.na(mean.alpha.est)){
    return(c("fail"))
  }
  y_poda.new<-y_poda+mean.alpha.est/(standarad_sigma1)
  test_x_1.new<-.hardPODA(NULL,y_poda.new,standarad_sigma1,lmd)
  gamma.j[,1]<-test_x_1.new$gamma
  
  
  ######## the second iteration #########
  
  mean.alpha.est<-median( -(omega.i[which(test_x_1.new$gamma==0)]))
  if(is.na(mean.alpha.est)){
    return(c("fail"))
  }
  y_poda.new<-y_poda+mean.alpha.est/(standarad_sigma1)
  test_x_1.new<-.hardPODA(NULL,y_poda.new,standarad_sigma1,lmd)
  gamma.j[,2]<-test_x_1.new$gamma
  
  
  ############ iterative loop ############
  
  
  Nind.tmp.len = list()
  ttime=300
  diff=rep(0,ttime)
  for(i in 1:ttime){
    if(length(which(test_x_1.new$gamma==0))==0){
      return(c("fail"))
    }else{
      gamma.j[,1] <-gamma.j[,2]
      mean.alpha.est<-median( -(omega.i[which(test_x_1$gamma==0)]))
      if(is.na(mean.alpha.est)){
        return(c("fail"))
      }
    }
    y_poda.new<-y_poda+mean.alpha.est/(standarad_sigma1)
    test_x_1.new<-.hardPODA(NULL,y_poda.new,standarad_sigma1,lmd)
    gamma.j[,2]<-test_x_1.new$gamma
    if(norm(as.matrix(gamma.j[,2]-gamma.j[,1],type="I"))<1e-5){
      break
    }
  }
  if(norm(as.matrix(gamma.j[,2]-gamma.j[,1],type="I"))>1e-5){
    print("unconvergence")
  }
  ReList=list(gamma=test_x_1.new$gamma,ress=y_poda.new)
  return(ReList)
}




# Tuning parameter selection through BIC criterion

.PODAlmdSelection = function(y_poda,omega.i,standarad_sigma1){
  method1="hard"
  m=length(y_poda)
  ########PODA find outlier#########
  r = 1
  N = length(y_poda)
  ress = NULL
  gammas = NULL
  betaInit = NULL
  lambdas = seq(ceiling(norm(matrix(y_poda, N, 1), "I")/1),
                0.1, by = -0.01)
  lambdas_discard=NULL
  for (sigma in lambdas) {  
    result = .PODAFUNnew(y_poda,omega.i,(standarad_sigma1),sigma)
    if(class(result)!="list"){
      lambdas_discard=c(lambdas_discard,sigma)
    }else{
      gammas = as.matrix(cbind(gammas, result$gamma))
      ress = cbind(ress, result$ress) 
    }
  }
  if(length(lambdas_discard)>0){
    lambdas=lambdas[-which(lambdas%in%lambdas_discard)]
  }
  DF = colSums(abs(gammas) > 1e-05) 
  sigmaSqEsts = colSums(((y_poda) %*% matrix(rep(1, ncol(gammas)), 
                                             1, ncol(gammas)) - gammas)^2)/(length(y_poda) - DF)
  sigmaSqEsts[sigmaSqEsts < 0] = 0
  sigmaSqEsts[sigmaSqEsts == "NaN"] = 0
  riskEst = ((N - r) * log(sigmaSqEsts * (N - r - DF)/(N -r)) + (log(N - r)) * (DF + 1))/(N - r)   #BIC
  riskEst[which(riskEst=="-Inf")]=1e5
  optSet = which(riskEst == min(riskEst))
  gammaOptInd = optSet[which(DF[optSet] == min(DF[optSet]))[1]]  
  gammaOpt = gammas[, gammaOptInd]
  resOpt = ress[, gammaOptInd]
  tau = mad(ress[gammas[, gammaOptInd] == 0, gammaOptInd])
  lmd=lambdas[gammaOptInd]
  resOpt.scale = resOpt/tau
  p = 2 * pnorm(-abs(resOpt.scale))
  out=lmd
  return(out)
}










# Feature/taxa selection function

.PODAFeatureSelection = function(y_poda,omega.i,standarad_sigma1,lmd){
  m=length(y_poda)
  gamma.j<-matrix(0,m,2)
  
  ######## the first iteration #########
  test_x_1<-.hardPODA(NULL,y_poda,standarad_sigma1,lmd) #hard
  mean.alpha.est<-median( -(omega.i[which(test_x_1$gamma==0)]))
  
  
  y_poda.new<-y_poda+mean.alpha.est/(standarad_sigma1)
  test_x_1.new<-.hardPODA(NULL,y_poda.new,standarad_sigma1,lmd)
  gamma.j[,1]<-test_x_1.new$gamma
  
  
  ######## the second iteration #########
  
  if(length(which(test_x_1.new$gamma==0))==0){
    return(NULL)
    break
  }else{
    mean.alpha.est<-median( -(omega.i[which(test_x_1.new$gamma==0)]))
  }
  y_poda.new<-y_poda+mean.alpha.est/(standarad_sigma1)
  test_x_1.new<-.hardPODA(NULL,y_poda.new,standarad_sigma1,lmd)
  gamma.j[,2]<-test_x_1.new$gamma
  
  ############ iterative loop ############
  
  
  Nind.tmp.len = list()
  for(i in 1:50){
    if(length(which(test_x_1.new$gamma==0))==0){
      return(NULL)
      break
    }else{
      gamma.j[,1] <-gamma.j[,2]
      mean.alpha.est<-median( -(omega.i[which(test_x_1.new$gamma==0)]))
    }
    y_poda.new<-y_poda+mean.alpha.est/(standarad_sigma1)
    test_x_1.new<-.hardPODA(NULL,y_poda.new,standarad_sigma1,lmd)
    gamma.j[,2]<-test_x_1.new$gamma
    
    Nind.tmp.len[[i]]=which(test_x_1.new$gamma!=0)
    if(norm(as.matrix(gamma.j[,2]-gamma.j[,1],type="I"))<1e-5){
      out=which(test_x_1.new$gamma!=0)
      break
    }
  }
  Nind=out
  return(Nind)
}








# Main function for PoDA
# @param Y  matrix; observed abundance matrix. Row: taxa; column: samples. NAs are not expected in observed abundance matrix.
# @param u  vector; the phenotype of interest.
# @param confoun matrix; the other phenotypes to be adjusted.
# @param rep_time  a positive integer; the number of repeated information splitting (rep_time > 1)
# @param fdrnomial a real value between 0 and 1; the nominal level of FDR control
# @param fdrmethod character; p-value adjusting approach. 
# @param lib_cut a non-negative real value; samples with less than lib.cut read counts are excluded.
# @param pre_cut a real value between 0 and 1; taxa with prevalence (percentage of nonzeros) less than prev.cut are excluded. 
# @param nCore a positive integer; the number of cores used for computing

# @output rej  vector; the index of rejected taxa
# @output rejtaxon.name  vector; the name of rejected taxa
# @output rej_atleast vector; the name of taxa which are rejected at least 1 time in multiple information splitting
# @output all_freq  data frame; the rejection frequency of all taxa
# @output rej_freq  data frame; the rejection frequency of rejected taxa




PODA = function( Y = Y, u = u, confoun = NULL ,rep_time=50,tau = 1,fdrnomial=0.05,fdrmethod="BH",lib_cut=1000, 
                 pre_cut=0.2,nCore=9){
  ##################################################################
  ####################### preparation ##############################
  ##################################################################
  lmdpi=0.5
  total_ind<-NULL
  rej_ind<-NULL
  sample_cut=lib_cut
  original.name=rownames(Y)

  
  
  n=ncol(Y) 
  m=nrow(Y)
  
  ##################################################################
  ###################### data preprocessing ########################
  ##################################################################
  rownames(Y)=paste("taxon",1:m)
  ############# sample filtration ###############
  if(is.null(confoun)){
    sample_fil = which(colSums(Y)>=sample_cut)
    Y=Y[,sample_fil]
    u=u[sample_fil]
    Z=as.matrix(u) 
    colnames(Z)=c("u")
    formula="u"
    formulaa="~u"
  }else{
    confoun=as.matrix(confoun)
    sample_fil = which(colSums(Y)>=sample_cut)
    Y=Y[,sample_fil]
    confoun=as.matrix(confoun[sample_fil,])
    colnames(confoun)=paste0("z",1:ncol(confoun))
    u=u[sample_fil]
    Z=cbind(u,confoun)
    colnames(Z)=c("u",paste0("z",1:ncol(confoun)))
    formula=paste(colnames(Z), collapse = "+")
    formulaa=paste("~",formula)
  }

  ############# taxa filtration ################ 
  datasize <- ncol(Y)
  prevalence = apply(as(Y, "matrix"), 1, function(x) {
    return(sum(x > 0))
  })/(datasize)
  
  keepOTUs = which(prevalence> pre_cut)
  Y = Y[keepOTUs,]
  
  
  
  n=ncol(Y) 
  m=nrow(Y) 
  retain_index=as.numeric(gsub("taxon","",rownames(Y)))
  retain_name=original.name[retain_index]



  
  ########## imputation ###########
  x_implement<-.imputation(Y,Z,formulaa)
  ####### clr transformation #########
  wi<-t(t(log(x_implement))-rowMeans(t(log(x_implement))))
  
  
  ##################################################################
  ###################### variance estimation #######################
  ##################################################################  
  
  C<-cbind(matrix(1,nrow=n),confoun)
  X1<-cbind(u,C)
  beta_hat<-wi%*%X1%*%solve(t(X1)%*%X1)  
  epsilon<-wi-beta_hat%*%t(X1)
  p<-dim(Z)[2]
  Sigma_hat<-(1/(n-p-1))*diag(epsilon%*%t(epsilon))
  standarad_sigma1<-sqrt(Sigma_hat) 
  sigma1<-standarad_sigma1**2
  std.sigma = standarad_sigma1
  
  
  
  
  
  ##################################################################
  ################### project the confounder #######################
  ################################################################## 
  
  H0<-C%*%solve(t(C)%*%C)%*%t(C)
  wi.center<-(diag(n)-H0)%*%t(wi)
  ui.center<-(diag(n)-H0)%*%as.matrix(u)
  
  
  
  ##################################################################
  ########################## choose lambda #########################
  ################################################################## 
  numerator<-t(wi.center)%*%ui.center
  dominator<-sum(ui.center**2)
  omega.i<-numerator/dominator    
  
  y_poda<-(as.matrix(omega.i)+median( -omega.i))/(standarad_sigma1) 
  
  lmd = .PODAlmdSelection( y_poda,omega.i,standarad_sigma1 )

  
  numCores <- nCore
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  ##################################################################
  ################## multiple information splitting ################
  ##################################################################
  
  result = foreach(flag = 1:rep_time,.combine='cbind',.packages = c('mvtnorm','foreach','MASS','matrixStats'),.export= c(".imputation",".PODAFeatureSelection",".hardPODA",".PODAFUNnew")) %dopar% {
    n=ncol(Y)
    m=nrow(Y)
    ind_result=rep(3,m)
    
    x_implement<-.imputation(Y,Z,formulaa)
    ########clr transformation##########
    wi<-t(t(log(x_implement))-rowMeans(t(log(x_implement))))
    ########Correction for heteroskedasticity##########
    C<-cbind(matrix(1,nrow=n),confoun)   
    X1<-cbind(u,C)
    beta_hat<-wi%*%X1%*%solve(t(X1)%*%X1)
    epsilon<-wi-beta_hat%*%t(X1)
    p<-dim(Z)[2]
    Sigma_hatt<-(1/(n-p-1))*(epsilon%*%t(epsilon))
    Sigma_hat<-(1/(n-p-1))*diag(epsilon%*%t(epsilon))
    standarad_sigma1<-sqrt(Sigma_hat)
    sigma1<-standarad_sigma1**2
    std.sigma = standarad_sigma1
    
    ############ information splitting ############
    wi_f = wi
    wi_g = wi
    
    set.seed(flag*100)
    all_fission_z <- t(rmvnorm(n, mean = rep(0, m), sigma = Sigma_hatt))
    wi_f <- wi_f + tau * all_fission_z
    wi_g <- wi_g - (1 / tau) * all_fission_z

    ######## project the confounder #########
    H0<-C%*%solve(t(C)%*%C)%*%t(C)
    wi.center_f<-(diag(n)-H0)%*%t(wi_f)
    wi.center_g<-(diag(n)-H0)%*%t(wi_g)
    ui.center<-(diag(n)-H0)%*%as.matrix(u)
    
    
    ######## variance estimation ##########
    beta_hat_g<-wi_g%*%X1%*%solve(t(X1)%*%X1)
    epsilon_g<-wi_g-beta_hat_g%*%t(X1)
    p<-dim(Z)[2]
    Sigma_hat_g<-(1/(n-p-1))*diag(epsilon_g%*%t(epsilon_g))
    standarad_sigma1_g<-sqrt(Sigma_hat_g)
    sigma1_g<-standarad_sigma1_g**2
    std.sigma_g = standarad_sigma1_g
    
    

    ######## taxon selection ##########
    wi.center = wi.center_f
    standarad_sigma1_f = sqrt(1+tau**2) * standarad_sigma1 
    Ys_tidle<-wi.center
    kappa<-c(t(ui.center)%*%(diag(n)-H0)%*%ui.center/(sum(ui.center**2))**2)
    numerator<-t(wi.center)%*%ui.center
    dominator<-sum(ui.center**2)
    omega.i<-numerator/dominator   
    y_poda<-(as.matrix(omega.i)+median( -omega.i))/(standarad_sigma1)
    tryselection <- try({
      Nind <- .PODAFeatureSelection(y_poda, omega.i, standarad_sigma1, lmd)
    }, silent = TRUE)

    if (inherits(tryselection, "try-error")) {
      return(NULL)
    }
  

    wi.center = wi.center_g
    standarad_sigma1 = standarad_sigma1_g 
    sigma1=standarad_sigma1**2
    
    
    if(length(Nind)==0){
      Nind=NULL 
      pp_kp=rep(0,m)
    }else{
      alpha_bar = -(sum(c(t(ui.center)%*%wi.center[,-Nind])/sigma1[-Nind]))/(sum(c(t(ui.center)%*%ui.center)/sigma1[-Nind])) 
      TEST = c(t(ui.center)%*%wi.center[,Nind])/c(t(ui.center)%*%ui.center)+alpha_bar
      
      TEST_var =  ((standarad_sigma1[Nind])**2/c(t(ui.center)%*%ui.center))
      pp_kp<-2*pnorm(-abs(TEST/sqrt(TEST_var))) 
      p.adj_kp<-p.adjust(pp_kp,fdrmethod) 
      
      ind_result[Nind]=pp_kp 
      ind_result[Nind[which(p.adj_kp<=fdrnomial)]]=1+ind_result[Nind[which(p.adj_kp<=fdrnomial)]] 
    }
    
    ind_result=c(length(Nind),ind_result)
    
    
    return(ind_result)
  }
  
  
  
  
  
  ##################################################################
  ########################## aggregation ###########################
  ##################################################################
  
  
  stopCluster(cl)
  result=as.matrix(result[,which(result[1,]!=0)]) 
  
  if(length(result[1,])==0){
    rej=NULL
    return(list(rej=NULL,rej_freq=NULL,rej_index=NULL,all_freq=NULL))
  }else{
    pp.value=as.matrix(result[-1,]) 
    Nind.mat = matrix(0,ncol=ncol(pp.value),nrow=nrow(pp.value)) 
    Nind.mat[(which(result[-1,]>=1&result[-1,]<2.5))]=1 
    pp.value[which( Nind.mat==1)]=pp.value[which( Nind.mat==1)]-1 
    len.pval = apply(pp.value,2,function(x){length(which(x<=(1-lmdpi)))}) 
    testfun=function(x){
      x[which(x<2.5)] = p.adjust(x[which(x<2.5)],fdrmethod)
      return(x)
    }
    qq.value = apply(pp.value,2,testfun) 
    qq.value[which(qq.value>=fdrnomial)]=NA 
    Nind =  which(rowSums(!is.na(qq.value))>0) 
    choose.matrix = matrix(0,nrow=nrow(result)-1,ncol=ncol(result)) 
    choose.matrix[ which(qq.value<fdrnomial) ]=1 
    choose.matrix=as.matrix(choose.matrix)
    inclu.rate = t(t(choose.matrix)/result[1,]) 
    inclu.rate = rowMeans(inclu.rate)
    tmp.ind = which(sort(inclu.rate)>0)[1]
    if(is.na(tmp.ind)){
      rej=NULL
    }else{
      tmp.sum = rep(0,m)
      tmp.sum2 = 0
      for(flag in tmp.ind:m){ 
        tmp.sum2 = tmp.sum2+sort(inclu.rate)[flag]
        tmp.sum[flag] = tmp.sum2 
      }
      tmp=min((1/(lmdpi)-mean( (len.pval)/(result[1,]*lmdpi))),1) 
      tmp.ind = sort(which(tmp.sum<=(tmp)*fdrnomial),decreasing=T)[1]
      tmp.sum2 = sort(inclu.rate)[tmp.ind]
      rej = which(inclu.rate>tmp.sum2)
      rej = as.numeric(gsub('[taxon]',' ',rownames(Y)[rej]))
    }
    rej_atleast=original.name[as.numeric(gsub('[taxon]',' ',rownames(Y)[Nind]))]
    rej_freq=data.frame(taxon=rej_atleast,frequency=rowSums(choose.matrix)[Nind])
    all_freq=data.frame(taxon=retain_name,frequency=rowSums(choose.matrix))
    rejtaxon.name=original.name[rej]
    
    return(list(rej=rej,rejtaxon.name=rejtaxon.name,rej_freq=rej_freq,rej_atleast=rej_atleast,all_freq=all_freq)) 
  }

 
}

