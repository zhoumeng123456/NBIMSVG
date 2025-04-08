
library(truncnorm)
library(nlme)
library(spatstat.data)
library(spatstat.geom)
library(spatstat.random)
library(spatstat.explore)
library(spatstat.model)
library(spatstat.linnet)
library(spatstat)
library(splines)
library(mvtnorm)
library(Rcpp)
library(statmod)
library(R.utils)
library(MCMCpack)
library(Matrix)
library(MASS)
library(parallel)

#library(devtools)
sim_create <- function(gene_size =100,svgene_size=0.1,sv_mark=c(1,1),no_sv_mark = c(0,0),inf_size=0,phi=100,etamean=2,
                       xspace="linear",yspace="linear",seed=1,cell_dist=rep(1,6),domainnum=1){
  ###
  ###
  set.seed(seed)
  x_coords <- rep(0:31, each = 32)
  y_coords <- rep(0:31, times = 32)
  if(domainnum==1){
    z_type <- c(rep(rep(1:2,each=16),9),rep(1,96),rep(c(rep(3,12),rep(4,20)),10),rep(3,128),rep(rep(3:4,each=16),6))#横平竖直的分成了四个区域
  }else if(domainnum==2){
    z_type <- rep(1,1024)
  }else if(domainnum==3){
    z_type <- c(rep(1,512),rep(2,512))
  }else if(domainnum==4){
    z_type <- c(rep(rep(1:2,each=16),16),rep(3,512))
  }else if(domainnum==5){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==6){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(3,320))
  }else if(domainnum==7){
    #判断是否是3号的偏差影响太大
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(5,320))
  }else if(domainnum==8){
    z_type <- c(rep(c(rep(2,12),rep(1,20)),16),rep(c(rep(2,12),rep(3,20)),6),rep(4,320))
  }else if(domainnum==9){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,12),rep(1,8),rep(4,12)),8),rep(rep(3:4,each=16),8))
  }else if(domainnum==10){
    z_type <- c(rep(c(rep(1,20),rep(2,12)),16),rep(c(rep(3,4),rep(1,16),rep(4,12)),4),rep(rep(3:4,each=16),12))
  }else if(domainnum==11){
    #
    z_type <- c(rep(rep(2:3,each=16),16),rep(1,512))
  }
  
  
  location <- as.matrix(data.frame(x = x_coords, y = y_coords,z=z_type))
  
  x <- location[,1]
  y <- location[,2]
  
  x <- x-mean(x)
  x <- x/sd(x)
  y <- y-mean(y)
  y <- y/sd(y)
  
  x <- switch(xspace, "focal" = exp(-x^2/2),
              "period" = cos(2*pi*x),
              "linear" = x,
              stop("Invalid xspace!"))
  
  y <- switch(yspace,
              "focal" = exp(-y^2/2),
              "period" = cos(2*pi*y),
              "linear" = y,
              stop("Invalid yspace!"))
  
  kern_coord<-cbind(x,y)
  
  npoints <- nrow(location)
  rownames(location) = paste('spot', 1:npoints, sep = '')
  expres_marx = as.data.frame(matrix(NA, nrow = npoints, ncol = gene_size))
  rownames(expres_marx) = paste('spot', 1:npoints, sep = '')
  colnames(expres_marx) = paste('gene', 1:gene_size, sep = '')
  
  sv_points=svgene_size*gene_size
  sv_gene <- c(1:sv_points)##设定前多少个为SV基因
  no_sv_gene <- setdiff(1:gene_size, sv_gene)##设定其余的基因为非SV基因
  
  eta <- rnorm(gene_size,mean = 2,sd = 0.5)
  
  cell_matrix <- matrix(NA,nrow = npoints,ncol = 6)
  
  for(i in 1:npoints){
    if(z_type[i]==1){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,1,1,1,1,1))
    }else if(z_type[i]==2){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,3,5,7,9,11))
    }else if(z_type[i]==3){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(14,12,10,8,6,4))
    }else if(z_type[i]==4){
      cell_matrix[i,] <- rdirichlet(1,alpha = c(1,4,4,4,4,1))
    }else{
      cell_matrix[i,] <- rdirichlet(1,alpha = c(18,16,14,12,10,8))
    }
  }
  
  cell_mark <- rnorm(6,mean=0,sd=1)
  
  for(i in sv_gene){
    for(t in 1:npoints){
      
      mu = exp(eta[i]+sum(kern_coord[t,]*sv_mark)+sum(cell_matrix[t,]*cell_mark))
      
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  for(i in no_sv_gene){
    for(t in 1:npoints){
      
      mu= exp(eta[i]+sum(kern_coord[t,]*no_sv_mark)+sum(cell_matrix[t,]*cell_mark))
      expres_marx[t,i] <- rnbinom(1,mu=mu,size=phi)
    }
  }
  
  ##添加零膨胀率。
  expres_marx <- as.matrix(expres_marx)
  total.size <- npoints*gene_size
  zero_size<-floor(npoints*gene_size*inf_size)
  zeroin<-sample(c(1:total.size),zero_size)
  expres_marx[zeroin]<-0
  location <- as.matrix(location)
  spe <- list(expres_marx,location[,1:2],cell_matrix)
  
  return(spe)
}




CTIG<-function(spe.list){
  m=length(spe.list)##获取数据集的数量
  Y.list<-list()#numpt list
  coord.list<-list()
  samplenumber<-c()
  
  for(i in 1:m){
    Y.list[[i]]<-spe.list[[i]][[1]]
    samplenumber[i]<-nrow(Y.list[[i]])## get samplenember of each datasets
    coord<-spe.list[[i]][[2]]
    coord <-coord-colMeans(coord)#Center  coordinates of spots 
    coord <- coord / apply(coord,2,sd)# normalize coordinates of spots
    coord.list[[i]]<-coord
  }
  
  G <- ncol(Y.list[[1]])## gene numbers##这里直接以第一个数据集为例，后部分需要修改
  
  ##对于coord.list进行样条转化,由此我们获得了
  co_spline_list_list <- list()##用于容纳
  for(j in c(1:4)){
    co_spline_list <- list()
    for(i in 1:m){
      coord_spline1 <- bs(x=coord.list[[i]][,1],df = j,degree=j)
      coord_spline2 <- bs(x = coord.list[[i]][,2],df = j,degree=j)
      co_spline_list[[i]] <- cbind(coord_spline1, coord_spline2)##每一个列表内的元素都是一个矩阵，分别为x和y的基函数值，列合并
    }
    co_spline_list_list[[j]] <- co_spline_list
  }
  
  result <- list(Y.list,co_spline_list_list,samplenumber,G)
  return(result)
}

f <- function(x, p, q, r, s, t,h=0) {
  result <- ((log(1+r*x))^q)*exp(p*log(x)+s*x*log(x)-s*lgamma(x)-t*x)
  #result <-  (log(1 + r * x/t))^q * exp(p*log(x)-x+s*log(x/t)*x/t-s*lgamma(x/t)+h)/(t^(p+1))
  return(result)
}

gausslq <- function(n=100,p,q,r,s,t,h=0){
  out <- gauss.quad.prob(n,"gamma",alpha =1,beta = 1 )
  return(sum(f(out[[1]],p,q,r,s,t,h)*out[[2]]))
}

compu_phi <- function(n=200,s,t,ap){
  h=0
  result=0
  while(is.nan(result)||result==0||result==Inf||result==-Inf){
    result <- gausslq(n,p=ap,q=0,r=1,s=s,t=t-1,h)/gausslq(n,p=ap-1,q=0,r=1,s=s,t=t-1,h)
    h=h+100
  }
  return(result)
}

safe_compu_phi <- function(n=200,s,t,ap,timeout=0.05,old=1){
  result <- old
  tryCatch(
    {result <- withTimeout(compu_phi(n=200,s=s,t=t,ap=ap),timeout = timeout)},
    TimeoutException = function(ex){result <- old}
  )
  return(result)
}

log_compu_phi <- function(s,t,ap){
  j=0#用来计数迭代了多少次
  result=0
  h=0
  while(is.nan(result) ||result==0||result==Inf||result==-Inf){
    result <- log(gausslq(n=100, p=ap-1, q=0, r=1, s=s, t=t-1, h=h))
    h=h+100
    j=j+1
    
  }
  return(result-100*j)
}

safe_log_compu_phi <- function(s,t,ap,timeout=0.05,old=0){
  result <- old
  tryCatch(
    {result <- withTimeout(log_compu_phi(s=s,t=t,ap=ap),timeout = timeout)},
    TimeoutException = function(ex){result <- old}
  )
  return(result)
}


result_built2 <- function(result,genenumber=100,ak=1,acr=0.0001,dp=20,cp=1,thrs=0.5,c_alpha){
  G <- result[[4]]
  total_result <- matrix(NA,nrow = genenumber,ncol = 2)
  for(g in c(1:genenumber)){
    begin=Sys.time()
    Ysingle.list<-list()
    for(i in 1:length(result[[3]])){
      Ysingle.list[[i]]<-as.vector(result[[1]][[i]][,g]) #先提取第一个基因
    }
    ite_result <- VIZINB(y = Ysingle.list,splinevalue = result[[2]],c_alpha = c_alpha,samplenumber = result[[3]],splinelevel = result[[5]],ak=ak,acrate = acr,dp=dp,cp=cp,thrs=thrs)
    total_result[g,1] <- ite_result[1]
    total_result[g,2] <- ite_result[2]
    end=Sys.time()
    print(paste0("第",g,"个基因用时:",end-begin))
  }
  total_result <- apply(total_result,1,max)
  total_result[total_result>0.5] <- 1
  total_result[total_result<=0.5] <- 0
  aaa <- rep(1,2)
  aaa[1] <- sum(total_result[1:30])/30
  aaa[2] <- sum(total_result[31:genenumber])/70
  return(aaa)
}

#生成门槛值
BayFDR <- function(PPI, alpha){
  PPI_sorted = sort(PPI,decreasing = TRUE)
  k = 1
  fdr = 0
  while(fdr < alpha){
    fdr = mean(1 - PPI_sorted[1:k])
    k = k+1
    if(k > length(PPI_sorted)){
      k = length(PPI_sorted);
      break;
    }
  }
  return.value = PPI_sorted[k]
  return.value = ifelse(is.na(return.value), 0, return.value)
  return(return.value)
}

compute_diag_ABC_corrected <- function(A, B) {
  n <- nrow(A)
  diag_elements <- numeric(n)
  
  for (i in 1:n) {
    diag_elements[i] <- A[i, ] %*% B %*% t(A[i, , drop=FALSE])
  }
  
  return(diag_elements)
}


##c_alpha是一个list，安装了细胞的组成矩阵
VIZINB <- function(y,splinevalue,samplenumber,splinelevel,acrate=0.01,ak=1,dp=1.8,cp=0.2,thrs=0.5,c_alpha,knum=6){
  ##处理输入
  M=length(samplenumber)
  G=splinelevel
  ##定义参数(初始化)
  
  u_phi <- rep(100,M)#初始化shape parameter
  phiold <- rep(15,M)
  logphiold <- rep(1,M)
  lcp <- rep(1,M)
  
  u_pi <- rep(list(),M)#初始化零膨胀
  for(i in 1:M){
    u_pi[[i]] <- rep(0.5,samplenumber[i])
  }
  #对u_pi进行初始化
  for(i in 1:M){
    for(j in 1:samplenumber[i]){
      if(y[[i]][j]==0){
        u_pi[[i]][j]=0.5
      }else{
        u_pi[[i]][j]=0
      }
    }
  }
  
  a_pi=1
  b_pi=1
  
  u_g <- rep(list(),M)#初始化泊松参数gi
  for(i in 1:M){
    u_g[[i]] <- rep(3,samplenumber[i])
  }
  
  u_log_g <- rep(list(),M)#初始化泊松参数gi的函数log
  for(i in 1:M){
    u_log_g[[i]] <- rep(2,samplenumber[i])
  }
  
  u_exp_g <- rep(list(),M)
  for(i in 1:M){
    u_exp_g[[i]] <- rep(5,samplenumber[i])
  }
  
  a_phi=rep(0.001,M)
  b_phi=rep(0.001,M)
  
  c_m <- rep(list(),M)##library size n^m*(2G+1)
  for(i in 1:M){
    c_m[[i]] <- cbind(rep(1,samplenumber[i]),splinevalue[[i]],c_alpha[[i]])
  }
  
  u_theta <- matrix(rep(0.2,(2*G+1+knum)*M),nrow=(2*G+1+knum),ncol=M)#初始化theta的均值
  u_theta[1,] <- 2
  sigma_theta <- rep(list(),M)#初始化theta的方差
  for(i in 1:M){
    sigma_theta[[i]] <- diag(rep(1,(2*G+1+knum)))
  }
  
  eta_sigma = rep(1,M)
  alpha_sigma=rep(1,M)
  
  u_alpha=matrix(rep(0,2*M),nrow = M,ncol = 2)#第一层indicator
  u_sigmaT_beta=matrix(rep(1,2*M),nrow = M,ncol = 2)#beta的方差
  
  
  Gamma1=0.01
  Gamma2=0.01
  
  
  u_abetaT=matrix(rep(1,2*M),nrow = M,ncol = 2)#beta方差的参数
  
  A_k=rep(ak,2)##尚未定义
  
  u_qk <- rep(0.5,2)
  u_pk <- rep(0.5,2)
  
  u_uk <- rep(0.5,2)
  
  c_q=1
  d_q=1
  c_p=cp
  d_p=dp
  
  ELBO=0
  ELBO_OLD=1000
  signo=100
  ##定义迭代
  num=0
  u_uk_old <- u_uk
  stability_counter <- 0##用于计数u_uk的稳定状态

  while(signo>acrate){
    ##更新u_uk
    for(k in 1:2){
      upart1=1
      upart2=1
      for(i in 1:M){
        upart1=upart1*(u_alpha[i,k]*u_qk[k]-u_alpha[i,k]+1)*(u_alpha[i,k]*u_qk[k]-u_qk[k]+1)
        upart2=upart2*(u_alpha[i,k]*Gamma2-u_alpha[i,k]+1)*(u_alpha[i,k]*Gamma2-Gamma2+1)
      }
      upart1=upart1*u_pk[k]
      upart2=upart2*(1-u_pk[k])
      u_uk[k]=upart1/(upart2+upart1)
    }
    
    #更新u_qk
    for(k in 1:2){
      q1=u_uk[k]*sum(u_alpha[,k])+c_q
      q2=u_uk[k]*M+d_q-u_uk[k]*sum(u_alpha[,k])
      u_qk[k]=q1/(q1+q2)
    } 
    
    #更新u_pk
    for(k in 1:2){
      u_pk[k]=(c_p+u_uk[k])/(1+d_p+c_p)
    }

    
    ##更新均值u_theta and sigma_theta（eta,beta1,beta2) and u_alpha
    
    for(i in 1:M){
      # ##计算M_q_sigma
      
      M_q_1=u_alpha[i,1]*u_sigmaT_beta[i,1]+(1-u_alpha[i,1])/(Gamma1^2)
      M_q_2=u_alpha[i,2]*u_sigmaT_beta[i,2]+(1-u_alpha[i,2])/(Gamma1^2)
      M_q = diag(c(1/((eta_sigma[i])^2),rep(M_q_1,G),rep(M_q_2,G),rep(1/(alpha_sigma[i]^2),knum)))
      #计算均值
      theta1=as.vector(exp(-c_m[[i]]%*%u_theta[,i]+unname(compute_diag_ABC_corrected(c_m[[i]],sigma_theta[[i]]))/2))
      temp1=(1-u_pi[[i]])*u_g[[i]]*theta1
      D_u_theta= u_phi[i]*(t(c_m[[i]])%*%(temp1-(1-u_pi[[i]])))-(M_q%*%u_theta[,i])
      
      #t(c_m[[i]])%*%diag(as.vector(temp1))
      
      #计算方差
      vec_D_sigma_theta=u_phi[i]*(t(c_m[[i]]*as.vector(temp1))%*%c_m[[i]])+M_q
      
      #迭代theta的q密度函数
      # if(abs(det(vec_D_sigma_theta))<1e-15){
      #   sigma_theta[[i]] = ginv(vec_D_sigma_theta)
      # }else{
      #   sigma_theta[[i]] = solve(vec_D_sigma_theta)##使用逆矩阵
      # }
      
      rigd=1e-8
      while(abs(det(vec_D_sigma_theta))<.Machine$double.eps){
        vec_D_sigma_theta=vec_D_sigma_theta+rigd*diag(1+2*G+knum)
        rigd <- rigd*10
      }
      
      sigma_theta[[i]] = solve(vec_D_sigma_theta)
      u_theta[,i] = u_theta[,i]+sigma_theta[[i]]%*%D_u_theta

      if(num<=2){
        theta1=as.vector(exp(-c_m[[i]]%*%u_theta[,i]+unname(compute_diag_ABC_corrected(c_m[[i]],sigma_theta[[i]]))/2))
      }##校准初始值

      #更新u_g，一个list，包含多个向量
      g1=(1-u_pi[[i]])+u_phi[i]*((1-u_pi[[i]])*theta1)
      g2=(y[[i]]+u_phi[i]-1)*(1-u_pi[[i]])+1
      u_g[[i]] <- as.vector(g2/g1)
      u_log_g[[i]] <- as.vector(digamma(g2)-log(g1))
      u_exp_g[[i]] <- as.vector((g1/(g1+1))^g2)
      
      
      #更新ri(m),更新u_pi,这是一个列表
      for(j in 1:samplenumber[i]){
        if(y[[i]][j]==0){
          u_pi[[i]][j]=beta(a_pi+1,b_pi)/(beta(a_pi+1,b_pi)+beta(a_pi,b_pi+1)*u_exp_g[[i]][j])
        }
      }
      

      ##更新u_phi,这是一个向量
      #计算c1(是否直接使用向量相乘)
      c1=as.numeric(t(1-u_pi[[i]])%*%(c_m[[i]]%*%u_theta[,i]-u_log_g[[i]])+b_phi[i]+
                      t((1-u_pi[[i]])*u_g[[i]])%*%theta1)
      #计算N_pi
      N_pi = sum(1-u_pi[[i]])
      #可以计算后验的分布
      u_phi[i]=safe_compu_phi(n=200,s=N_pi,t=c1,ap=a_phi[i],old = u_phi[i])
      phiold[i] <- u_phi[i]
      
      #更新u_sigmaT_beta M*2的矩阵
      ##当k=1时
      sigma1=2*u_abetaT[i,1]+u_alpha[i,1]*(sum(u_theta[2:(G+1),i]^2)+sum(unname(diag(sigma_theta[[i]]))[2:(G+1)]))
      u_sigmaT_beta[i,1]= (G*u_alpha[i,1]+1)/sigma1
      ##当k=2时
      sigma2=2*u_abetaT[i,2]+u_alpha[i,2]*(sum(u_theta[(G+2):(2*G+1),i]^2)+sum(unname(diag(sigma_theta[[i]]))[(G+2):(2*G+1)]))
      u_sigmaT_beta[i,2]= (G*u_alpha[i,2]+1)/sigma2
      
      
      #更新u_aphiT,u_abetaT
      #u_aphiT[i]=1/(u_phi[i]+M_phi)
      u_abetaT[i,1]=1/(1/(A_k[1]^2)+u_sigmaT_beta[i,1])
      u_abetaT[i,2]=1/(1/(A_k[2]^2)+u_sigmaT_beta[i,2])
      
      
      #更新u_alpha M*2的矩阵
      for(k in 1:2){
        part1=(u_qk[k]*u_uk[k]+1-u_uk[k])*(Gamma2-Gamma2*u_uk[k]+u_uk[k])
        part2=((1-u_qk[k])*u_uk[k]+1-u_uk[k])*((1-Gamma2)*(1-u_uk[k])+u_uk[k])
        if(k==1){
          part5=exp(lgamma((G+G*u_alpha[i,k]+1)/2)-lgamma((G*u_alpha[i,k]+1)/2)-G/2*log(u_alpha[i,k]*(sum(u_theta[2:(G+1),i]^2)+sum(unname(diag(sigma_theta[[i]]))[2:(G+1)]))/2+u_abetaT[i,k]))
        }else{
          part5=exp(lgamma((G+G*u_alpha[i,k]+1)/2)-lgamma((G*u_alpha[i,k]+1)/2)-G/2*log(u_alpha[i,k]*(sum(u_theta[(G+2):(2*G+1),i]^2)+sum(unname(diag(sigma_theta[[i]]))[(G+2):(2*G+1)]))/2+u_abetaT[i,k]))
        }

        if(k==1){
              theta_sample <- mvrnorm(n=1,mu=u_theta[2:(G+1),i],Sigma = sigma_theta[[i]][2:(G+1),2:(G+1)])
        }else{
              theta_sample <- mvrnorm(n=1,mu=u_theta[(G+2):(2*G+1),i],Sigma = sigma_theta[[i]][(G+2):(2*G+1),(G+2):(2*G+1)])
        }
        #分子那项
        part3=((u_alpha[i,k]*sum(theta_sample^2)+2*u_abetaT[i,k])/((u_alpha[i,k]+1)*sum(theta_sample^2)+2*u_abetaT[i,k]))^(G*u_alpha[i,k]/2+0.5)
        #分母那项
        part4=exp(log(part2)-G*log(Gamma1)-log(part1)-log(part5)-sum(theta_sample^2)/2/Gamma1^2)
        temp_alpha=part3/(part3+part4)
        u_alpha[i,k]=temp_alpha
      }

      #计算ELBO
      ##第一部分
      ELBO=t(1-u_pi[[i]])%*%(-lgamma(y[[i]]+1)-u_phi[i]*(c_m[[i]]%*%u_theta[,i]))+
        -0.5*sum(diag(M_q%*%(u_theta[,i]%*%t(u_theta[,i])+sigma_theta[[i]])))
      -(1-u_alpha[i,1])*G*log(Gamma1)-(1-u_alpha[i,2])*G*log(Gamma1)
      +(1-u_uk[1])*(u_alpha[i,1]*log(Gamma2)+(1-u_alpha[i,1])*log(1-Gamma2))
      +(1-u_uk[2])*(u_alpha[i,2]*log(Gamma2)+(1-u_alpha[i,2])*log(1-Gamma2))
      -b_phi[i]*u_phi[i]
      
      for(j in 1:samplenumber[i]){
        ELBO=ELBO+u_pi[[i]][j]*log(beta(a_pi+1,b_pi))-log(u_pi[[i]][j]^(u_pi[[i]][j]))+(1-u_pi[[i]][j])*log(beta(a_pi,b_pi+1))-log((1-u_pi[[i]][j])^(1-u_pi[[i]][j]))
      }
      
      ##第二部份
      ELBO=ELBO+sum(-u_log_g[[i]]+lgamma(g2)
                    -g2*log((1-u_pi[[i]])*(u_phi[i]*theta1+1)))
      +c1*u_phi[i]+safe_log_compu_phi(s=N_pi,t=c1,ap=a_phi[i])+0.5*log(det(sigma_theta[[i]]))
      ##加入带有k=1/2部分
      for(k in 1:2){
        if(k==1){
          e_betak=sum(u_theta[2:(G+1),i]^2)+sum(unname(diag(sigma_theta[[i]])[2:(G+1)]))
        }else{
          e_betak=sum(u_theta[(G+2):(2*G+1),i]^2)+sum(unname(diag(sigma_theta[[i]])[(G+2):(2*G+1)]))
        }
        ELBO=ELBO-(u_alpha[i,k]*G+1)/2*log(u_abetaT[i,k]+u_alpha[i,k]*e_betak/2)+lgamma((u_alpha[i,k]*G+1)/2)
        +(u_alpha[i,k]*e_betak/2+u_abetaT[i,k])*u_sigmaT_beta[i,k] 
        -log(u_sigmaT_beta[i,k]+1/A_k[k]^2)
        -log(u_alpha[i,k]^(u_alpha[i,k]))-log((1-u_alpha[i,k])^((1-u_alpha[i,k])))
      }

    }
    ELBO=ELBO-log(u_uk[1]^u_uk[1])-log((1-u_uk[1])^((1-u_uk[1])))-log(u_uk[2]^u_uk[2])-log((1-u_uk[2])^((1-u_uk[2])))+
      log(beta(u_uk[1]*sum(u_alpha[,1])+c_q,u_uk[1]*(M-sum(u_alpha[,1]))+d_q))+log(beta(u_uk[1]+c_p,d_p-u_uk[1]+1))+
      log(beta(u_uk[2]*sum(u_alpha[,2])+c_q,u_uk[2]*(M-sum(u_alpha[,2]))+d_q))+log(beta(u_uk[2]+c_p,d_p-u_uk[2]+1))
    signo=abs(ELBO-ELBO_OLD)
    ELBO_OLD=ELBO
    
    
    
    if (all(abs(u_uk - u_uk_old) < 1e-9)) {  # 设定一个容忍度阈值
      stability_counter <- stability_counter + 1
    } else {
      stability_counter <- 0  # 如果 u_uk 发生变化，重置计数器
    }
    
    u_uk_old <- u_uk  # 更新旧的 u_uk 值
    
    if (stability_counter > 8 && signo < 1) {  
      break
    }
    
    num=num+1
    if(num>200){
      break
    }
    
    print(paste0(num," u1: ",u_uk[1],"u2:",u_uk[2]," error: ",round(signo,4)))
  }
  
  ##检测函数
  return(u_uk)
}

citgtest <- function(g,result,c_alpha,spe.list){
  Ysingle.list<-list()
  for(i in 1:length(result[[3]])){
    Ysingle.list[[i]]<-as.vector(result[[1]][[i]][,g]) #先提取第一个基因
  }
  begin <- Sys.time()
  tun_choose <- tun_spl(spe.list = spe.list,calpha.list = c_alpha,g=g)##得到对应的样条基阶数，这是一个数字
  
  zero_ratios <- sapply(Ysingle.list, function(x) mean(x == 0))
  # 检查零值比例的均值是否大于0.8
  #获得对应的b样条展开
  if(tun_choose==1){
    ak1=0.08
  }else if(tun_choose==2){
    ak1=0.05
  }else if(tun_choose==3){
    ak1=0.04
  }else{
    ak1=0.03
  }
  
  ##选择对应的ak1
  knum1=dim(c_alpha[[1]])[2]
  bbb <- VIZINB(Ysingle.list,result[[2]][[tun_choose]],result[[3]],tun_choose,acrate=0.01,ak=ak1,dp=1.8,cp=0.2,thrs=0.5,c_alpha=c_alpha,knum=knum1)
  end <- Sys.time()
  print(paste0("the",g,"gene use the time:",end-begin))
  return(bbb)
}

tun_spl <- function(spe.list,calpha.list,g=1){
  ##检查长度
  n1=length(spe.list)
  n2=length(calpha.list)
  if(n1==n2){
  }else{
    print("warning!the length of expression matrix do not match the calpha's")
  }
  tunning_choose=rep(1,n1)
  for(iii in c(1:n1)){
    aic_value <- c(1:4)
    y = unname(spe.list[[iii]][[1]][, g])
    
    # 计算零值比例
    zero_count <- sum(y == 0)
    samplenum=length(y)
    print(zero_count)
    total_count <- length(y)
    zero_ratio <- zero_count / total_count
    zero_ratio
    # 如果零值比例超过 90%，直接设置 min_position 为 1
    if (zero_ratio > 0.9) {
      min_position <- 1
    }else if(zero_ratio <0.7 && samplenum > 500 ){
      for (kkk in c(1:4)){
        splinelevel = kkk
        degreelevel = kkk
        coord_spline1 <- bs(x = spe.list[[iii]][[2]][, 1], df = splinelevel, degree = degreelevel)
        coord_spline2 <- bs(x = spe.list[[iii]][[2]][, 2], df = splinelevel, degree = degreelevel)
        x <- cbind(coord_spline1, coord_spline2)
        
        rows_with_sum_1 <- rowSums(calpha.list[[iii]]) == 1
        calpha.list[[iii]][rows_with_sum_1, 1] <- calpha.list[[iii]][rows_with_sum_1, 1] - 0.01
        calpha.list[[iii]][calpha.list[[iii]] < 0] <- 0
        
        data <- data.frame(x = x, calpha = calpha.list[[iii]], y = y)
        x_names = paste0("x", c(1:(2 * splinelevel)))
        celltypenum = dim(calpha.list[[iii]])[2]
        calpha_names = paste0("calpha", c(1:celltypenum))
        data_name = c(x_names, calpha_names, "y")
        names(data) = data_name
        formula_str <- paste(data_name[1:(2*splinelevel+celltypenum)], collapse = " + ")
        formula_count <- paste("y ~", formula_str)
        formula_zero <- paste("|", formula_str)
        formula_full <- as.formula(paste(formula_count, formula_zero))
        model_zinb <- zeroinfl(formula_full, data = data, dist = "negbin")
        summary(model_zinb)
        aic_value[kkk] <- AIC(model_zinb)
      }
      min_position <- which.min(aic_value)
    }else{
      for (kkk in c(1:4)) {
        splinelevel = kkk
        degreelevel = kkk
        coord_spline1 <- bs(x = spe.list[[iii]][[2]][, 1], df = splinelevel, degree = degreelevel)
        coord_spline2 <- bs(x = spe.list[[iii]][[2]][, 2], df = splinelevel, degree = degreelevel)
        x <- cbind(coord_spline1, coord_spline2)
        
        rows_with_sum_1 <- rowSums(calpha.list[[iii]]) == 1
        calpha.list[[iii]][rows_with_sum_1, 1] <- calpha.list[[iii]][rows_with_sum_1, 1] - 0.01
        calpha.list[[iii]][calpha.list[[iii]] < 0] <- 0
        
        data <- data.frame(x = x, calpha = calpha.list[[iii]], y = y)
        x_names = paste0("x", c(1:(2 * splinelevel)))
        celltypenum = dim(calpha.list[[iii]])[2]
        calpha_names = paste0("calpha", c(1:celltypenum))
        data_name = c(x_names, calpha_names, "y")
        names(data) = data_name
        # 使用 paste 函数生成字符串
        formula_str <- paste(data_name[1:(2*splinelevel+celltypenum)], collapse = " + ")
        # 零膨胀模型拟合
        formula_count <- paste("y ~", formula_str)
        formula_zero <- paste("|", formula_str)
        formula_full <- as.formula(paste(formula_count, formula_zero))
        model_zinb <- hurdle(formula_full, data = data, dist = "negbin")
        summary(model_zinb)
        aic_value[kkk] <- AIC(model_zinb)
      }
      min_position <- which.min(aic_value)
    }
    
    tunning_choose[iii]=min_position
  }
  tunning_final=max(tunning_choose)##选取最大值
  print(paste0("本list选择的b样条阶数为",tunning_final))
  # print(aic_value)
  # print(bic_value)
  # print(paste("aic选择的阶数是:", min_position))
  # print(paste("bic选择的阶数是:", min_position1))
  return(tunning_final)
}

CTIG<-function(spe.list){
  m=length(spe.list)##获取数据集的数量
  Y.list<-list()#numpt list
  coord.list<-list()
  samplenumber<-c()
  for(i in 1:m){
    Y.list[[i]]<-spe.list[[i]][[1]]#点位✖️基因
    samplenumber[i]<-nrow(Y.list[[i]])## get samplenember of each datasets
    coord<-spe.list[[i]][[2]]
    coord <-coord-colMeans(coord)#Center  coordinates of spots 
    coord <- coord / apply(coord,2,sd)# normalize coordinates of spots
    coord.list[[i]]<-coord##收集到了
  }
  
  G <- ncol(Y.list[[1]])## gene numbers##这里直接以第一个数据集为例，后部分需要修改
  
  ##对于coord.list进行样条转化,由此我们获得了
  co_spline_list_list <- list()##用于容纳
  for(j in c(1:4)){
    co_spline_list <- list()
    for(i in 1:m){
      coord_spline1 <- bs(x=coord.list[[i]][,1],df = j,degree=j)
      coord_spline2 <- bs(x = coord.list[[i]][,2],df = j,degree=j)
      co_spline_list[[i]] <- cbind(coord_spline1, coord_spline2)##每一个列表内的元素都是一个矩阵，分别为x和y的基函数值，列合并
    }
    co_spline_list_list[[j]] <- co_spline_list
  }
  
  result <- list(Y.list,co_spline_list_list,samplenumber,G)
  return(result)
}

ourmethod <- function(spelist,c_alpha,num_cores=23){
  cl <- makeCluster(num_cores)
  
  # 在每个节点加载所需的库
  clusterEvalQ(cl, {
    library(truncnorm)
    library(nlme)
    library(spatstat.data)
    library(spatstat.geom)
    library(spatstat.random)
    library(spatstat.explore)
    library(spatstat.model)
    library(spatstat.linnet)
    library(splines)
    library(mvtnorm)
    library(Rcpp)
    library(statmod)
    library(R.utils)
    library(MCMCpack)
    library(Matrix)
    library(MASS)
    library(pscl)
  })
  
  # 导出必要的函数和对象到并行环境
  clusterExport(cl, c("sim_create", "CTIG", "f", "gausslq", "compu_phi", 
                      "safe_compu_phi", "log_compu_phi", "safe_log_compu_phi", 
                      "result_built2", "BayFDR", "compute_diag_ABC_corrected", 
                      "VIZINB", "citgtest", "c_alpha", "spelist","tun_spl"))
  
  # 检查基因数量是否一致
  genenum <- length(colnames(spelist[[1]][[1]]))
  genename <- colnames(spelist[[1]][[1]])
  
  # 初始化 result 对象
  result_ctig <- CTIG(spe.list = spelist)
  
  # 并行处理的输入数据
  input_data <- 1:genenum
  chunk_size <- 100  # 每次并行处理100个基因
  total_result <- matrix(NA, nrow = genenum, ncol = 2)  # 用于存储结果的矩阵
  progress_intervals <- seq(0.1, 1, by = 0.1)  # 定义进度的间隔
  current_progress <- 0  # 当前进度
  
  # 并行分块处理
  for (i in seq(1, length(input_data), by = chunk_size)) {
    chunk_data <- input_data[i:min(i + chunk_size - 1, length(input_data))]
    
    # 使用 parLapply 并行处理每个块
    result_chunk <- parLapply(cl, chunk_data, function(x) {
      result <- tryCatch({
        withTimeout({
          citgtest(x, result = result_ctig, c_alpha = c_alpha, spe.list = spelist)
        }, timeout = 600, onTimeout = "error")
      }, error = function(e) {
        return(c(NA, NA))  # 错误处理，返回 NA
      })
      return(result)
    })
    
    # 将结果写入总结果矩阵
    for (j in seq_along(chunk_data)) {
      total_result[chunk_data[j], ] <- result_chunk[[j]]
    }
    
    # 计算并打印进度
    completed_fraction <- i / length(input_data)
    if (completed_fraction >= progress_intervals[current_progress + 1]) {
      current_progress <- current_progress + 1
      print(paste0(current_progress * 10, "% completed"))
    }
  }
  
  alpha=0.05/2
  gamma0=1/21
  out1=total_result
  total_result <- apply(total_result, 1, max)
  out2 <- total_result
  temp1=total_result
  temp=1-temp1
  temp=temp[temp>=gamma0]
  sorted_temp <- sort(temp)
  threshold_index <- ceiling(alpha/(1-alpha)/20*length(temp))
  threshold_value <- sorted_temp[threshold_index]
  out3=1-threshold_value
  temp1[temp1<=(1-threshold_value)]=0
  temp1[temp1>(1-threshold_value)]=1
  out4 = temp1
  svgenename <- as.vector(na.omit(genename[temp1 == 1]))
  out5=svgenename
  result_total <- list(out1, out2, out3, out4, out5)
  # 关闭集群
  stopCluster(cl)
  return(result_total)
}






