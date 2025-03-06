# Y = Y
# u = u
# confun = confun
# lmdpi=0.5
# rep_time=500
# tau = 1
# fdrnomial=0.05
# lib_cut=0
# pre_cut=0
# nCore=9

## Y: relative abundance matrix row: col
## u: covariates of interest.
## confun: other confounders
## lmdpi
## rep_time
## tau
## fdrnomial
## lib_cut
## pre_cut
# ## nCore
# Y = Y
# u = u
# confun = NULL 
# lmdpi=0.5
# rep_time=50
# tau = 1
# fdrnomial=0.05 
# lib_cut=1000
# pre_cut=0.2
# nCore=9



PODA = function( Y = Y, u = u, confun = NULL ,lmdpi=0.5,rep_time=50,tau = 1,fdrnomial=0.05 ,lib_cut=1000, 
                 pre_cut=0.2,nCore=9){
  ##################################################################
  ######################### 准备工作 ###############################
  ##################################################################
  total_ind<-NULL
  rej_ind<-NULL
  sample_cut=lib_cut
  original.name=rownames(Y)
  
  
  
  n=ncol(Y) #样本数
  m=nrow(Y) #taxa数
  
  ##################################################################
  ########################## 处理数据 ##############################
  ##################################################################
  
  ############# 过滤sample ###############
  if(is.null(confun)){
    sample_fil = which(colSums(Y)>=sample_cut)
    Y=Y[,sample_fil]
    u=u[sample_fil]
    Z=as.matrix(u) 
    colnames(Z)=c("u")
    formula="u"
    formulaa="~u"
  }else{
    confun=as.matrix(confun)
    sample_fil = which(colSums(Y)>=sample_cut)
    Y=Y[,sample_fil]
    confun=as.matrix(confun[sample_fil,])
    colnames(confun)=paste0("z",1:ncol(confun))
    u=u[sample_fil]
    Z=cbind(u,confun)
    colnames(Z)=c("u",paste0("z",1:ncol(confun)))
    formula=paste(colnames(Z), collapse = "+")
    formulaa=paste("~",formula)
  }
  
  ############# 过滤taxa ################ 
  datasize <- ncol(Y)
  prevalence = apply(as(Y, "matrix"), 1, function(x) {
    return(sum(x > 0))
  })/(datasize)
  
  keepOTUs = which(prevalence> pre_cut)
  Y = Y[keepOTUs,]
  
  
  
  n=ncol(Y) #新的sample数
  m=nrow(Y) #新的taxa数
  
  rownames(Y)=paste("taxon",1:m)
  
  
  ########## 0的imputation ###########
  x_implement<-.imputation(Y,Z,formulaa)
  
  
  ####### clr transformation #########
  ws<-t(t(log(x_implement))-rowMeans(t(log(x_implement))))
  
  
  ##################################################################
  ##################### 估计方差bar{sigma} #########################
  ##################################################################  
  
  C<-cbind(matrix(1,nrow=n),confun)      #confun
  X1<-cbind(u,C)
  beta_hat<-ws%*%X1%*%solve(t(X1)%*%X1)  
  epsilon<-ws-beta_hat%*%t(X1)
  p<-dim(Z)[2]
  Sigma_hat<-(1/(n-p-1))*diag(epsilon%*%t(epsilon))
  standarad_sigma1<-sqrt(Sigma_hat) ######## 原始标准差
  sigma1<-standarad_sigma1**2 ######## 原始方差
  std.sigma = standarad_sigma1 #########原始标准差存档
  
  
  
  
  
  ##################################################################
  ################### project the confounder #######################
  ################################################################## 
  
  H0<-C%*%solve(t(C)%*%C)%*%t(C)
  ws.center<-(diag(n)-H0)%*%t(ws)
  us.center<-(diag(n)-H0)%*%as.matrix(u)
  
  
  
  ##################################################################
  ############## 开始初次outlier detection 选择lambda ##############
  ################################################################## 
  numerator<-t(ws.center)%*%us.center
  dominator<-sum(us.center**2)
  omega.i<-numerator/dominator    #我理论里的y_ipod
  
  y_ipod<-(as.matrix(omega.i)+median( -omega.i))/(standarad_sigma1) #初值选所有估计值的中值
  
  lmd = .IPODlmdSelection( y_ipod,omega.i,standarad_sigma1 )
  print(lmd)
  
  
  
  ##################################################################
  ###############开始初次outlier detection 选择lambda###############
  ################################################################## 
  
  
  
  # variables_to_export <- ls()
  
  numCores <- nCore
  cl <- makeCluster(numCores)
  registerDoParallel(cl)
  ##################################################################
  ########################### 开始MDF ##############################
  ##################################################################
  
  result = foreach(flag = 1:rep_time,.combine='cbind',.packages = c('mvtnorm','foreach','MASS','matrixStats'),.export= c(".imputation",".IPODFeatureSelection",".hardIPOD","IPODFUNnew")) %dopar% {
    n=ncol(Y)
    m=nrow(Y)
    ind_result=rep(3,m)
    
    
    ########################################## 就多了这些 ################################################
    x_implement<-.imputation(Y,Z,formulaa)
    ########clr transformation##########
    ws<-t(t(log(x_implement))-rowMeans(t(log(x_implement))))
    ########Correction for heteroskedasticity##########
    C<-cbind(matrix(1,nrow=n),confun)      #confun
    X1<-cbind(u,C)
    beta_hat<-ws%*%X1%*%solve(t(X1)%*%X1)
    epsilon<-ws-beta_hat%*%t(X1)
    p<-dim(Z)[2]
    Sigma_hatt<-(1/(n-p-1))*(epsilon%*%t(epsilon))
    Sigma_hat<-(1/(n-p-1))*diag(epsilon%*%t(epsilon))
    standarad_sigma1<-sqrt(Sigma_hat)
    sigma1<-standarad_sigma1**2
    std.sigma = standarad_sigma1
    
    ############ 用standarad sigma 做fission ############
    ws_f = ws
    ws_g = ws
    
    set.seed(flag*100)
    # all_fission_z <- t(rmvnorm(n, mean = rep(0, m), sigma = diag(sigma1)))
    all_fission_z <- t(rmvnorm(n, mean = rep(0, m), sigma = Sigma_hatt))# Generate all random samples at once
    ws_f <- ws_f + tau * all_fission_z
    ws_g <- ws_g - (1 / tau) * all_fission_z
    
    # for(i in 1:m){
    #   set.seed(i+100*flag)
    #   # set.seed(i+(100*flag+1000*term))
    #   fission.z = rmvnorm(1,mean=rep(0,n),sigma1[i]*(diag(n)))
    #   ws_f[i,]= ws_f[i,] + tau *fission.z
    #   ws_g[i,]= ws_g[i,] - (1/tau) *fission.z
    # }

    ########project the confounder#########
    H0<-C%*%solve(t(C)%*%C)%*%t(C)
    ws.center_f<-(diag(n)-H0)%*%t(ws_f)
    ws.center_g<-(diag(n)-H0)%*%t(ws_g)
    us.center<-(diag(n)-H0)%*%as.matrix(u)
    
    
    ######## 用ws_g估计方差 ##########
    beta_hat_g<-ws_g%*%X1%*%solve(t(X1)%*%X1)
    epsilon_g<-ws_g-beta_hat_g%*%t(X1)
    p<-dim(Z)[2]
    Sigma_hat_g<-(1/(n-p-1))*diag(epsilon_g%*%t(epsilon_g))
    standarad_sigma1_g<-sqrt(Sigma_hat_g)
    sigma1_g<-standarad_sigma1_g**2
    std.sigma_g = standarad_sigma1_g
    
    
    # ######## 用ws_f估计方差 ##########
    # beta_hat_f<-ws_f%*%X1%*%solve(t(X1)%*%X1)
    # epsilon_f<-ws_f-beta_hat_f%*%t(X1)
    # p<-dim(Z)[2]
    # Sigma_hat_f<-(1/(n-p-1))*diag(epsilon_f%*%t(epsilon_f))
    # standarad_sigma1_f<-sqrt(Sigma_hat_f)
    # sigma1_f<-standarad_sigma1_f**2
    # std.sigma_f = standarad_sigma1_f   
    
    
    
    ######## 用ws.center_f选择变量 ##########
    ws.center = ws.center_f
    
    
    
    
    standarad_sigma1_f = sqrt(1+tau**2) * standarad_sigma1 ############### 这里！！！
    # standarad_sigma1 = standarad_sigma1_f
    
    ########Change loss function#########
    Ys_tidle<-ws.center
    kappa<-c(t(us.center)%*%(diag(n)-H0)%*%us.center/(sum(us.center**2))**2)
    
    #########mean(alpha)初值取0
    yinz<-sqrt(2 * log(m))
    yinz = 1
    
    yinzi<-1
    numerator<-t(ws.center)%*%us.center
    dominator<-sum(us.center**2)
    omega.i<-numerator/dominator    #我理论里的y_ipod
    y_ipod<-(as.matrix(omega.i)+median( -omega.i))/(standarad_sigma1)
    tryselection <- try({
      Nind <- .IPODFeatureSelection(y_ipod, omega.i, standarad_sigma1, lmd)
    }, silent = TRUE)

    if (inherits(tryselection, "try-error")) {
      return(NULL)
    }
    # Nind = .IPODFeatureSelection(y_ipod,omega.i,standarad_sigma1,lmd)
    
    # ##########data fission###############
    ws.center = ws.center_g
    standarad_sigma1 = standarad_sigma1_g #检验时用从g中估计出的方差保证独立性
    sigma1=standarad_sigma1**2
    
    
    if(length(Nind)==0){
      Nind=NULL #如果集合为空，都是non-DA
      pp_kp=rep(0,m)
    }else{
      alpha_bar = -(sum(c(t(us.center)%*%ws.center[,-Nind])/sigma1[-Nind]))/(sum(c(t(us.center)%*%us.center)/sigma1[-Nind])) #联合估计
      TEST = c(t(us.center)%*%ws.center[,Nind])/c(t(us.center)%*%us.center)+alpha_bar
      
      TEST_var =  ((standarad_sigma1[Nind])**2/c(t(us.center)%*%us.center))
      pp_kp<-2*pnorm(-abs(TEST/sqrt(TEST_var))) #算p值
      p.adj_kp<-p.adjust(pp_kp,"BH") #BH procedure
      
      ind_result[Nind]=pp_kp #Nind 非0部分记录原始p值,用于估计集合中原假设大小
      ind_result[Nind[which(p.adj_kp<=0.05)]]=1+ind_result[Nind[which(p.adj_kp<=0.05)]] # 被调整后的p值 拒绝的部分在原来的基础上+1
    }
    
    ind_result=c(length(Nind),ind_result)
    
    
    return(ind_result)
  }
  
  
  
  
  
  ##################################################################
  ######################## foreach结束啦 ###########################
  ##################################################################
  
  
  stopCluster(cl)
  result=as.matrix(result[,which(result[1,]!=0)]) #只挑那些有Nind的部分来看，去掉啥也没选出来的结果
  
  if(length(result[1,])==0){
    rej=NULL
    return(list(rej=NULL,rej_freq=NULL,rej_index=NULL,all_freq=NULL))
  }else{
    
    pp.value=as.matrix(result[-1,]) #第一行是每一次 N_can的大小
    
    Nind.mat = matrix(0,ncol=ncol(pp.value),nrow=nrow(pp.value)) #一个pp.value一样大的0矩阵
    Nind.mat[(which(result[-1,]>=1&result[-1,]<2.5))]=1 #记录被拒绝的玩意
    
    pp.value[which( Nind.mat==1)]=pp.value[which( Nind.mat==1)]-1 #被拒绝的部分p值还原回去 这些都是真正的原始p值了
    
    
    
    len.pval = apply(pp.value,2,function(x){length(which(x<=(1-lmdpi)))}) #原始p值与0.5做比较 来估计集合大小了
    
    qq.value = pp.value #复制一下原始p值矩阵
    
    
    
    ############# 要开始 aggregation 了 #############3
    
    
    
    testfun=function(x){
      x[which(x<2.5)] = p.adjust(x[which(x<2.5)],"BH")
      return(x)
    }
    
    qq.value = apply(pp.value,2,testfun) #对原始p值重新校正，qq.value里都是校正后的p值和3
    
    qq.value[which(qq.value>=0.05)]=NA #接受原假设的部分被赋为NA
    
    Nind =  which(rowSums(!is.na(qq.value))>0) #至少被拒绝过一次的合集
    
    choose.matrix = matrix(0,nrow=nrow(result)-1,ncol=ncol(result)) ######## 行taxon 列 DF次数
    choose.matrix[ which(qq.value<0.05) ]=1 ######## 拒绝原假设的部分赋为1
    
    choose.matrix=as.matrix(choose.matrix)
    
    #### 得到inclution rate矩阵 ###
    inclu.rate = t(t(choose.matrix)/result[1,]) 
    inclu.rate = rowMeans(inclu.rate)
    tmp.ind = which(sort(inclu.rate)>0)[1] #从此处开始，inclu.rate非0
    
    
    
    if(is.na(tmp.ind)){ #没有非0的部分 全部被接受了
      rej=NULL
    }else{
      
      tmp.sum = rep(0,m)
      tmp.sum2 = 0
      for(flag in tmp.ind:m){ #从非0部分开始加
        tmp.sum2 = tmp.sum2+sort(inclu.rate)[flag]
        tmp.sum[flag] = tmp.sum2 #记录每一次的加和
      }
      tmp=min((1/(lmdpi)-mean( (len.pval)/(result[1,]*lmdpi))),1) #校正因子
      # print(tmp)
      tmp.ind = sort(which(tmp.sum<=(tmp)*fdrnomial),decreasing=T)[1]
      tmp.sum2 = sort(inclu.rate)[tmp.ind] #找出满足条件最小sort
      
      if(tmp.sum2==0){ #如果是0，那就all，否则
        rej = which(inclu.rate>tmp.sum2)
        rej = as.numeric(gsub('[taxon]',' ',rownames(Y)[rej]))
        print("all")
      }else{
        rej = which(inclu.rate>tmp.sum2)
        rej = as.numeric(gsub('[taxon]',' ',rownames(Y)[rej]))
        print(c("length",length(which(inclu.rate==tmp.sum2))))
  
      }
    }
    rejtaxon = as.matrix(rowSums(choose.matrix)[Nind])
    rej_index=Nind
    all_freq=rowSums(choose.matrix)
    rejtaxon.name=original.name[rej]
    #### 可能要补充一下和原来的核对一下 输出taxon名和指标
    
    return(list(rej=rej,rejtaxon.name=rejtaxon.name,rej_freq=rejtaxon,rej_index=rej_index,all_freq=all_freq)) #暂时先只返回结果，之后要什么再补再优化吧,缺一个配对的结果，晚点补，今天先把真实数据写掉
  }
  
  # return(rej)
 
}