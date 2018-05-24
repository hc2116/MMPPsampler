utils::globalVariables(c("x", "perc", "alpha", "group"))

# Accumulation ######################################################

Forward_accumulation <- function(y_0T,Inter,Q,Lambda,messages=TRUE){
  Rho <- MASS::Null(Q)/sum(MASS::Null(Q))
  intermed <- Forward_Accumulation_cpp(y_0Tprime=y_0T,Inter=Inter,
                                       Q=Q,Lambda=Lambda,Rho=Rho,
                                       messages=messages,pp=12)
  return(list(x=intermed$x,LLH=intermed$LLH))
}

Interval_sampling_accumulation <- function(x_0T,y_0T.prime,
                                               Inter,Q,Lambda,
                                               G=NULL,messages=TRUE){
  if(messages==TRUE){cat("\r Interval sampling")}
  Tempres <- Interval_sampling_accumulation_cpp(x_0T=x_0T,
                                                y_0Tprime=y_0T.prime,
                                                Inter=Inter,
                                                Q=Q,
                                                Lambda=Lambda,messages=messages,pp=12)
  
  return(list(x=c(Tempres$x),t=c(Tempres$t),Count=c(Tempres$Count),R_n=Tempres$R_n))
}

Z_Sampling <- function(y_0T,path_x,path_t,
                             Lambda,Lambda_Y,Inter,messages=TRUE){
  if(messages==TRUE){cat("\r Latent-variable sampling")}
  intermed <- Z_Sampling_cpp(y_0T,path_x,path_t,Inter,
                             Lambda_Y,Lambda,
                             messages=messages)
  return(t(intermed$z))
}

GibbsSampler_hierarchical <- function(y_0T,M,Inter,
                                          alpha_Gamma_rate,
                                          beta_Gamma_rate,
                                          alpha_Gamma_Q,
                                          beta_Gamma_Q,
                                          alpha_Gamma_Y,
                                          beta_Gamma_Y,
                                          alpha_D = NULL,
                                          B=20,N=100,
                                          messages=TRUE){
  {
    if(class(M)!="numeric"|M<1|M != round(M)){
      cat("M must be integer! breaking")
      stop()}
    if(Inter<1|class(Inter)!="numeric"){
      cat("Value of Inter smaller than 1 or not numeric! breaking \n")
      stop()
    }
    if(length(y_0T)<2|class(y_0T)!="numeric"|min(y_0T)<0){
      cat("y_0T too short or not numeric! breaking \n")
      stop()
    }
    if(length(alpha_Gamma_rate)<M|class(alpha_Gamma_rate)!="numeric"){
      cat("alpha_Gamma_rate !=M or not numeric --> breaking \n")
      stop()
    }
    if(length(alpha_Gamma_Q)<M|class(alpha_Gamma_Q)!="numeric"){
      cat("alpha_Gamma_Q!=M or not numeric --> breaking \n")
      stop()
    }
    if(length(beta_Gamma_rate)<M|class(beta_Gamma_rate)!="numeric"){
      cat("beta_Gamma_rate!=M or not numeric --> breaking \n")
      stop()
    }
    if(length(beta_Gamma_Q)<M|class(alpha_Gamma_Q)!="numeric"){
      cat("beta_Gamma_Q!=M or not numeric --> breaking \n")
      stop()
    }
    if(length(alpha_Gamma_Y)!=1|class(alpha_Gamma_Y)!="numeric"){
      cat("alpha_Gamma_Y!=1 or not numeric --> breaking \n")
      stop()
    }
    if(length(beta_Gamma_Y)!=1|class(beta_Gamma_Y)!="numeric"){
      cat("beta_Gamma_Y!=1 or not numeric --> breaking \n")
      stop()
    }
    
    if(is.null(alpha_D)){
      alpha_D <- rep(1,M-1)
    }
    if(length(alpha_D)<(M-1)|class(alpha_Gamma_rate)!="numeric"){
      cat("alpha_D!=M-1 or not numeric--> breaking \n")
      stop()
    }
    
  }
  
  ptm_total <- proc.time()
  cpu_time <- c(0,0,0,0)
  t <- length(y_0T)
  if(alpha_Gamma_rate[1]==0){
    Lambda_S <- 0.00000001
  }else{
    Lambda_S <- stats::rgamma(1,shape=alpha_Gamma_rate[1],rate=beta_Gamma_rate[1])
  }
  Lambda_S <- matrix(sort(c(Lambda_S,stats::rgamma(M-1,shape=alpha_Gamma_rate[2:M],rate=beta_Gamma_rate[2:M]))),nrow=1)
  Lambda_Y_S <- stats::rgamma(1,shape=alpha_Gamma_Y,rate=beta_Gamma_Y)
  
  Q_S_r <- matrix(stats::rgamma(M,shape=alpha_Gamma_Q[1:M],rate=beta_Gamma_Q[1:M]),nrow=1)
  Q_S_d <- gtools::rdirichlet(M,alpha=alpha_D)*Q_S_r[1,]
  
  LLH=NULL
  x_samples=NULL
  Z_S <- NULL
  
  Q_S <- matrix(0,nrow=M,ncol=M)
  for(jj in 1:M){
    Q_S[jj,-jj] <- Q_S_d[jj+(1-1)*M,]
    Q_S[jj,jj] <- -Q_S_r[1,jj]
  }
  
  FW <- Forward_accumulation(ceiling(y_0T/10),Inter, Q=Q_S,
                             Lambda=Lambda_S[1,],messages=messages)
  path <- Interval_sampling_accumulation(FW$x,ceiling(y_0T/10),Inter = Inter,Q=Q_S,
                                             Lambda = Lambda_S[1,],
                                             messages=messages)
  ptm <- proc.time()
  for(j in 1:(N+B)){
    if(messages==TRUE){cat(paste("\r Sample trajectory",j,":                            \n"))}
    Z_intermed <- Z_Sampling(y_0T=y_0T,path_x=path$x,path_t=path$t,
                             Lambda=Lambda_S[j,],
                             Lambda_Y=Lambda_Y_S[j],
                             Inter=Inter,
                             messages=messages)
    Z_S <- rbind(Z_S,Z_intermed)
    cpu_time[1] <- cpu_time[1]+(proc.time()-ptm)[3]
    ptm <- proc.time()
    
    FW <- Forward_accumulation(Z_S[j,],Inter, Q=Q_S,
                               Lambda=Lambda_S[j,],
                               messages=messages)
    x_0T.s <- FW$x
    cpu_time[2] <- cpu_time[2]+(proc.time()-ptm)[3]
    ptm <- proc.time()
    path <- Interval_sampling_accumulation(x_0T.s,Z_S[j,],Inter = Inter,Q=Q_S,
                                               Lambda = Lambda_S[j,],
                                               messages=messages)
    #if(messages==TRUE){cat("\r Interval sampling done                      ")}
    cpu_time[3] <- cpu_time[3]+(proc.time()-ptm)[3]
    x_0T.star <- path$x
    times <- diff(path$t)
    x_samples <- cbind(x_samples, x_0T.s)
    
    if(j>B){
      LLH <- c(LLH,sum(stats::dpois(y_0T-Z_S[j,],lambda=Z_S[j,]*Lambda_Y_S[j],log=TRUE)))
    }
    State_time <- NULL
    for(jj in 1:M){
      State_time <- c(State_time,sum(times[x_0T.star[-length(x_0T.star)]==jj]))      
    }
    Changes <- matrix(0,nrow=M,ncol=M)
    for(ll in 1:M){
      for(mm in 1:M){
        Changes[ll,mm] <- sum(x_0T.star[-length(x_0T.star)]==ll&x_0T.star[-1]==mm)
      }
    }
    if(sum(path$Count<0)|sum(State_time<0)){
      if(messages==TRUE){cat("Counts or State_time below zero\n")}
      path$Count[path$Count<0]=0
      State_time[State_time<0]=0
    }
    if(alpha_Gamma_rate[1]==0){
      Lambda_S1 <- 0.00000001
    }else{
      Lambda_S1 <- stats::rgamma(1,shape=alpha_Gamma_rate[1]+path$Count[1],rate=beta_Gamma_rate[1]+State_time[1])
    }
    Lambda_S <- rbind(Lambda_S,
                      c(Lambda_S1,
                        stats::rgamma(M-1,shape=alpha_Gamma_rate[2:M]+path$Count[2:M],rate=beta_Gamma_rate[2:M]+State_time[2:M])))
    
    
    Q_S_r <- rbind(Q_S_r,stats::rgamma(M,shape=alpha_Gamma_Q[1:M]+rowSums(Changes)-diag(Changes),rate=beta_Gamma_Q[1:M]+State_time[1:M]))
    
    for(jj in 1:M){
      Q_S_d <- rbind(Q_S_d,gtools::rdirichlet(1,alpha=alpha_D+Changes[jj,-jj])*Q_S_r[j+1,jj])
    }
    Lambda_Y_S <- c(Lambda_Y_S,stats::rgamma(1,shape=alpha_Gamma_Y+sum(y_0T-Z_S[j,]),rate=beta_Gamma_Y+sum(Z_S[j,])))
  }
  if(messages==TRUE){cat(paste("\r Done                             "))}
  cpu_time[3] <- (proc.time()-ptm_total)[3]
  return(list(x=x_samples[,-(1:B)],
              Inter=Inter,
              y_0T=y_0T,
              Y=Z_S[,-(1:B)],
              cpu_time=cpu_time,
              Lambda_S=Lambda_S[-(1:B),],
              Q_r=Q_S_r[-(1:B),],
              Q_d=Q_S_d[-(1:M*B),],
              Lambda_Y=Lambda_Y_S[-(1:B)],
              LLH=LLH))
}

GibbsSampler <- function(y_0T,M,Inter,
                                          alpha_Gamma_rate,
                                          beta_Gamma_rate,
                                          alpha_Gamma_Q,
                                          beta_Gamma_Q,
                                          alpha_D = NULL,
                                          B=20,N=100,
                                          messages=TRUE){
  if(B<=0){
    B=1
  }
  if(N<=0){
    N=1
  }
  
  {
    if(class(M)!="numeric"|M<1|M != round(M)){
      cat("M must be integer! breaking")
      stop()}
    if(Inter<1|class(Inter)!="numeric"){
      cat("Value of Inter smaller than 1 or not numeric! breaking \n")
      stop()
    }
    if(length(y_0T)<2|class(y_0T)!="numeric"|min(y_0T)<0){
      cat("y_0T too short or not numeric! breaking \n")
      stop()
    }
    if(length(alpha_Gamma_rate)<M|class(alpha_Gamma_rate)!="numeric"){
      cat("alpha_Gamma_rate !=M or not numeric --> breaking \n")
      stop()
    }
    if(length(alpha_Gamma_Q)<M|class(alpha_Gamma_Q)!="numeric"){
      cat("alpha_Gamma_Q!=M or not numeric --> breaking \n")
      stop()
    }
    if(length(beta_Gamma_rate)<M|class(beta_Gamma_rate)!="numeric"){
      cat("beta_Gamma_rate!=M or not numeric --> breaking \n")
      stop()
    }
    if(length(beta_Gamma_Q)<M|class(alpha_Gamma_Q)!="numeric"){
      cat("beta_Gamma_Q!=M or not numeric --> breaking \n")
      stop()
    }
    
    if(is.null(alpha_D)){
      alpha_D <- rep(1,M-1)
    }
    if(length(alpha_D)<(M-1)|class(alpha_Gamma_rate)!="numeric"){
      cat("alpha_D!=M-1 or not numeric--> breaking \n")
      stop()
    }
  }
  
  ptm_total <- proc.time()
  cpu_time <- c(0,0,0,0)
  t <- length(y_0T)
  Lambda_S <- matrix(sort(stats::rgamma(M,shape=alpha_Gamma_rate[1:M],rate=beta_Gamma_rate[1:M])),nrow=1)
  Q_S_r <- matrix(stats::rgamma(M,shape=alpha_Gamma_Q[1:M],rate=beta_Gamma_Q[1:M]),nrow=1)
  Q_S_d <- gtools::rdirichlet(M,alpha=alpha_D)*Q_S_r[1,]
  
  LLH=NULL
  x_samples=NULL

  Q_S <- matrix(0,nrow=M,ncol=M)
  for(jj in 1:M){
    Q_S[jj,-jj] <- Q_S_d[jj+(1-1)*M,]
    Q_S[jj,jj] <- -Q_S_r[1,jj]
  }
  
  ptm <- proc.time()
  for(j in 1:(N+B)){
    if(messages==TRUE){cat(paste("\r Sample trajectory",j,":                            \n"))}
    cpu_time[1] <- cpu_time[1]+(proc.time()-ptm)[3]
    ptm <- proc.time()
    
    FW <- Forward_accumulation(y_0T,Inter, Q=Q_S,
                               Lambda=Lambda_S[j,],
                               messages=messages)
    x_0T.s <- FW$x
    cpu_time[2] <- cpu_time[2]+(proc.time()-ptm)[3]
    ptm <- proc.time()
    path <- Interval_sampling_accumulation(x_0T.s,y_0T,Inter = Inter,Q=Q_S,
                                               Lambda = Lambda_S[j,],
                                               messages=messages)
    cpu_time[3] <- cpu_time[3]+(proc.time()-ptm)[3]
    x_0T.star <- path$x
    times <- diff(path$t)
    x_samples <- cbind(x_samples, x_0T.s)
    
    State_time <- NULL
    for(jj in 1:M){
      State_time <- c(State_time,sum(times[x_0T.star[-length(x_0T.star)]==jj]))      
    }
    Changes <- matrix(0,nrow=M,ncol=M)
    for(ll in 1:M){
      for(mm in 1:M){
        Changes[ll,mm] <- sum(x_0T.star[-length(x_0T.star)]==ll&x_0T.star[-1]==mm)
      }
    }
    
    if(j>B){
      LLH <- c(LLH,sum(stats::dpois(path$Count[1:M],lambda=State_time*Lambda_S[j,],log=TRUE)))
    }
    if(sum(path$Count<0)|sum(State_time<0)){
      if(messages==TRUE){cat("Counts or State_time below zero\n")}
      path$Count[path$Count<0]=0
      State_time[State_time<0]=0
    }
    Lambda_S <- rbind(Lambda_S,
                        stats::rgamma(M,shape=alpha_Gamma_rate[1:M]+path$Count[1:M],rate=beta_Gamma_rate[1:M]+State_time[1:M]))
    
    Q_S_r <- rbind(Q_S_r,stats::rgamma(M,shape=alpha_Gamma_Q[1:M]+rowSums(Changes)-diag(Changes),rate=beta_Gamma_Q[1:M]+State_time[1:M]))
    
    for(jj in 1:M){
      Q_S_d <- rbind(Q_S_d,gtools::rdirichlet(1,alpha=alpha_D+Changes[jj,-jj])*Q_S_r[j+1,jj])
    }
  }
  if(messages==TRUE){cat(paste("\r Done                             "))}
  cpu_time[3] <- (proc.time()-ptm_total)[3]
  return(list(x=x_samples[,-(1:B)],
              y_0T=y_0T,
              Inter=Inter,
              cpu_time=cpu_time,
              Lambda_S=Lambda_S[-(1:B),],
              Q_r=Q_S_r[-(1:B),],
              Q_d=Q_S_d[-(1:M*B),],
              LLH=LLH))
}

.onUnload <- function (libpath) {
  library.dynam.unload("MMPPsampler", libpath)
}



# Plotting ##########################################################
modeesthenry <- function(x){
  return(as.numeric(rownames(table(x))[which.max(table(x))]))
}

modepercthenry <- function(x){
  return(1-max(table(x))/length(x))
}
stateperchenry <- function(x,M){
  return(tabulate(x,nbins=M)/length(x))
}

MMPPplot <- function(Sampler_Output=NULL,
                     title=" ",xaxis=" ",breaks=NULL,
                     colour=NULL){
  
  if(!is.null(Sampler_Output)){
    y_0T=Sampler_Output$y_0T
    M=dim(Sampler_Output$Q_r)[2]
    traj_samples=Sampler_Output$x
    Inter=Sampler_Output$Inter
  }
  Intervals <- seq(0,length(y_0T)*Inter,Inter)
  
  stateperchenry_temp <- function(x){
    return(stateperchenry(x,M=M))
  }
  c25 <- c("dodgerblue2","#E31A1C", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
           "black","gold1",
           "skyblue2","#FB9A99", # lt pink
           "palegreen2",
           "#CAB2D6", # lt purple
           "#FDBF6F", # lt orange
           "gray70", "khaki2",
           "maroon","orchid1","deeppink1","blue1","steelblue4",
           "darkturquoise","green1","yellow4","yellow3",
           "darkorange4","brown","dodgerblue2","#E31A1C", # red
           "green4",
           "#6A3D9A", # purple
           "#FF7F00", # orange
           "black","gold1",
           "skyblue2","#FB9A99", # lt pink
           "palegreen2",
           "#CAB2D6", # lt purple
           "#FDBF6F", # lt orange
           "gray70", "khaki2",
           "maroon","orchid1","deeppink1","blue1","steelblue4",
           "darkturquoise","green1","yellow4","yellow3",
           "darkorange4","brown")
  
  
  
  df1 = data.frame(t=Intervals[-1],
                   x=y_0T)
  df2 = data.frame(t=Intervals,
                   x=apply(traj_samples,1,modeesthenry),
                   x2=rowMeans(traj_samples),
                   sd=apply(traj_samples,1,stats::sd),
                   perc=apply(traj_samples,1,modepercthenry))
  df3 <- data.frame(t=rep(Intervals,M),
                    x=rep.int(1:M,rep(length(Intervals),M)),
                    alpha=c(t(apply(traj_samples,1,stateperchenry_temp))),
                    group=as.character(rep.int(1:M,rep(length(Intervals),M))))
  if(!is.null(colour)){
    Exp_times <- colour[colour>min(Intervals)&colour<max(Intervals)]
    Exp_indices <- 1
    Phases <- c(" ")
    for(i in 1:length(Exp_times)){
      Exp_indices <- c(Exp_indices,which(Intervals==min(Intervals[Intervals[-1]>Exp_times[i]])))
      Phases <- c(Phases,paste("Phase, ",i))
    }
    Exp_indices <- c(Exp_indices,length(Intervals))
    df1$group <- rep(Phases,times=diff(Exp_indices))
    
    p1 <- ggplot2::ggplot(df1,ggplot2::aes(x=t,y=x,group=group))+
      ggplot2::geom_ribbon(ggplot2::aes(x=t,ymax=max(x),ymin=x,fill=group),alpha=0.3)+
      ggplot2::scale_fill_manual(values= c25[1:length(Phases)])+
      ggplot2::guides(fill=FALSE)+
      ggplot2::geom_line()+
      ggplot2::theme_bw()+
      ggplot2::labs(x="time",y="Events",title=title)+
      ggplot2::theme(legend.title = ggplot2::element_blank(),
            axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank())
  }else{
    p1 <- ggplot2::ggplot(df1,ggplot2::aes(x=t,y=x))+
      ggplot2::geom_line()+
      ggplot2::theme_bw()+
      ggplot2::labs(x="time",y="Events",title=title)+
      ggplot2::theme(legend.title = ggplot2::element_blank(),
            axis.title.x=ggplot2::element_blank(),
            axis.text.x=ggplot2::element_blank(),
            axis.ticks.x=ggplot2::element_blank())
  }
  
  p2 = ggplot2::ggplot(df2,ggplot2::aes(x=t,y=x))+
    ggplot2::geom_line(data = df3,ggplot2::aes(x=t,y=x,group=group,size=2*alpha,alpha=0<alpha), color="firebrick1")+
    ggplot2::geom_line(color="#619CFF",size=1.2)+ #"firebrick1 or 2"
    ggplot2::theme_bw()+
    ggplot2::labs(x="time",y="Process state",title="Estimated Markov-process")+
    ggplot2::theme(legend.title = ggplot2::element_blank(),
          axis.title.x=ggplot2::element_blank(),
          axis.text.x=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank())+
    ggplot2::theme(legend.position="None")
  
  p3 = ggplot2::ggplot(df2,ggplot2::aes(x=t,y=perc))+
    ggplot2::geom_segment(ggplot2::aes(yend=perc,xend=t),color="coral2",yend=0)+
    ggplot2::theme_bw()+
    ggplot2::labs(x=xaxis,y="Uncertainty")+
    ggplot2::theme(legend.title = ggplot2::element_blank())+
    ggplot2::theme(legend.position="None")
  
  if(!is.null(breaks)){
    p1 <- p1+ ggplot2::scale_x_continuous(breaks=breaks[,1],
                                 labels=breaks[,2])
    p2 <- p2+ ggplot2::scale_x_continuous(breaks=breaks[,1],
                                 labels=breaks[,2])
    p3 <- p3+ ggplot2::scale_x_continuous(breaks=breaks[,1],
                                 labels=breaks[,2])
  }
  
  cowplot::plot_grid(p1, p2, p3, align = "v", nrow = 3, rel_heights = c(8/19, 8/19, 5/19))
}



