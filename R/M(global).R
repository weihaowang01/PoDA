####################################
############# hard IPOD ############
####################################
.hardIPOD = function(X, Y,standarad_sigma,lmd){
  r = Y
  N = length(Y)
  gamma = matrix(rep(0, N), N, 1)
  theta = lmd/standarad_sigma
  gamma[abs(r) > theta] = r[abs(r) > theta]
  return(list(gamma=gamma, ress = r))
  
}





# ###################################################
# ##### 应该是在这个函数里进行Algorithm 1的循环 #####
# ###################################################


IPODFUNnew = function(y_ipod,omega.i,standarad_sigma1,lmd){
  m=length(y_ipod)
  yinz = 1
  yinzi=1
  reg.inv<-1
  gamma.i<-matrix(0,m,2)
  
  
  ########IPOD find outlier#########
  test_x_1<-.hardIPOD(NULL,y_ipod,standarad_sigma1,lmd) #hard
  mean.alpha.est<-median( -(omega.i[which(test_x_1$gamma==0)]))
  if(is.na(mean.alpha.est)){
    return(c("nope"))
  }
  y_ipod.new<-y_ipod+mean.alpha.est/(standarad_sigma1)
  test_x_1.new<-.hardIPOD(NULL,y_ipod.new,standarad_sigma1,lmd)
  reg.inv=1
  gamma.i[,1]<-test_x_1.new$gamma
  
  
  #######################################迭代第二次
  
  mean.alpha.est<-median( -(omega.i[which(test_x_1.new$gamma==0)]))
  if(is.na(mean.alpha.est)){
    return(c("nope"))
  }
  y_ipod.new<-y_ipod+mean.alpha.est/(standarad_sigma1*yinzi)
  test_x_1.new<-.hardIPOD(NULL,y_ipod.new,standarad_sigma1,lmd)
  gamma.i[,2]<-test_x_1.new$gamma
  
  
  ####################################开始循环
  
  
  Nind.tmp.len = list()
  ttime=300
  diff=rep(0,ttime)
  for(i in 1:ttime){
    
    if(length(which(test_x_1.new$gamma==0))==0){
      return(c("nope"))
    }else{
      gamma.i[,1] <-gamma.i[,2]
      mean.alpha.est<-median( -(omega.i[which(test_x_1$gamma==0)]))
      if(is.na(mean.alpha.est)){
        return(c("nope"))
      }
    }
    y_ipod.new<-y_ipod+mean.alpha.est/(standarad_sigma1*yinzi)
    test_x_1.new<-.hardIPOD(NULL,y_ipod.new,standarad_sigma1,lmd)
    reg.inv=1
    gamma.i[,2]<-test_x_1.new$gamma
    
    if(norm(as.matrix(gamma.i[,2]-gamma.i[,1],type="I"))<1e-5){
      break
    }
  }
  
  if(norm(as.matrix(gamma.i[,2]-gamma.i[,1],type="I"))>1e-5){
    print("unconvergence")
  }
  
  
  
  ReList=list(gamma=test_x_1.new$gamma,ress=y_ipod.new)
  
  return(ReList)
}




######################################################
################## IPODlmdSelection ##################
######################################################

.IPODlmdSelection = function(y_ipod,omega.i,standarad_sigma1){
  method1="hard"
  yinz = 1
  yinzi=1
  
  m=length(y_ipod)
  gamma.i<-matrix(0,m,2)
  ########IPOD find outlier#########
  
  #我的X就是NULL
  r = 1
  N = length(y_ipod)
  ress = NULL
  gammas = NULL
  
  betaInit = NULL
  
  #给lambda的取值打网格
  lambdas = seq(ceiling(norm(matrix(y_ipod, N, 1), "I")/1),
                0.1, by = -0.01)
  
  lambdas_discard=NULL
  
  for (sigma in lambdas) {    #把lambda的值全部取一遍
    result = IPODFUNnew(y_ipod,omega.i,(standarad_sigma1),sigma)
    if(class(result)!="list"){
      lambdas_discard=c(lambdas_discard,sigma)
    }else{
      gammas = as.matrix(cbind(gammas, result$gamma))
      ress = cbind(ress, result$ress) #把每一次的结果放在一起
    }
  }
  if(length(lambdas_discard)>0){
    lambdas=lambdas[-which(lambdas%in%lambdas_discard)]
  }
  
  
  DF = colSums(abs(gammas) > 1e-05) #计数每一列有多少个gamma不为0 计true的个数
  
  
  sigmaSqEsts = colSums(((y_ipod) %*% matrix(rep(1, ncol(gammas)), 
                                             1, ncol(gammas)) - gammas)^2)/(length(y_ipod) - DF)
  
  sigmaSqEsts[sigmaSqEsts < 0] = 0
  sigmaSqEsts[sigmaSqEsts == "NaN"] = 0
  riskEst = ((N - r) * log(sigmaSqEsts * (N - r - DF)/(N -r)) + (log(N - r)) * (DF + 1))/(N - r)   #BIC
  riskEst[which(riskEst=="-Inf")]=1e5
  
  
  
  
  optSet = which(riskEst == min(riskEst))
  gammaOptInd = optSet[which(DF[optSet] == min(DF[optSet]))[1]]   #取risk最小的
  gammaOpt = gammas[, gammaOptInd]
  resOpt = ress[, gammaOptInd]
  tau = mad(ress[gammas[, gammaOptInd] == 0, gammaOptInd])
  lmd=lambdas[gammaOptInd]
  resOpt.scale = resOpt/tau
  p = 2 * pnorm(-abs(resOpt.scale))
  
  out=lmd
  
  return(out)
}










######################################################
################# IPODFeatureSelection ###############
######################################################


.IPODFeatureSelection = function(y_ipod,omega.i,standarad_sigma1,lmd){
  m=length(y_ipod)
  yinz = 1
  yinzi=1
  reg.inv<-1
  gamma.i<-matrix(0,m,2)
  
  ########IPOD find outlier#########
  test_x_1<-.hardIPOD(NULL,y_ipod,standarad_sigma1,lmd) #hard
  mean.alpha.est<-median( -(omega.i[which(test_x_1$gamma==0)]))
  
  
  y_ipod.new<-y_ipod+mean.alpha.est/(standarad_sigma1)
  test_x_1.new<-.hardIPOD(NULL,y_ipod.new,standarad_sigma1,lmd)
  reg.inv=1
  gamma.i[,1]<-test_x_1.new$gamma
  
  
  #######################################迭代第二次
  
  if(length(which(test_x_1.new$gamma==0))==0){
    return(NULL)
    break
  }else{
  mean.alpha.est<-median( -(omega.i[which(test_x_1.new$gamma==0)]))
  }
  y_ipod.new<-y_ipod+mean.alpha.est/(standarad_sigma1*yinzi)
  test_x_1.new<-.hardIPOD(NULL,y_ipod.new,standarad_sigma1,lmd)
  gamma.i[,2]<-test_x_1.new$gamma
  
  ####################################开始循环
  
  
  Nind.tmp.len = list()
  diff=rep(0,50)
  for(i in 1:50){
    
    if(length(which(test_x_1.new$gamma==0))==0){
      return(NULL)
      break
    }else{
      gamma.i[,1] <-gamma.i[,2]
      mean.alpha.est<-median( -(omega.i[which(test_x_1.new$gamma==0)]))
    }
    y_ipod.new<-y_ipod+mean.alpha.est/(standarad_sigma1*yinzi)
    test_x_1.new<-.hardIPOD(NULL,y_ipod.new,standarad_sigma1,lmd)
    reg.inv=1
    gamma.i[,2]<-test_x_1.new$gamma
    
    Nind.tmp.len[[i]]=which(test_x_1.new$gamma!=0)
    diff[i]=norm(gamma.i[,2]-gamma.i[,1],type="2")
    if(norm(as.matrix(gamma.i[,2]-gamma.i[,1],type="I"))<1e-5){
      out=which(test_x_1.new$gamma!=0)
      break
    }
  }
  Nind=out
  return(Nind)
}

