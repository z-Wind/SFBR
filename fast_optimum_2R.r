#需設定電壓公式，電阻％數與電阻範圍
#輸入理想電壓，差異值，樣本數

options(keep.source = FALSE)

#計算電壓值
voltage <- function(r1,r2)
{
  ref = 1.25
  if(length(r1) > 1)
  {
      ref <- rnorm(length(r1),ref,0.007/3)
  }
  return (ref*(1+r1/r2)) #更改成正確電壓公式
}


#計算合適值
fitness <- function(mean,sd,target_v,bond)
{
    return ((bond-abs(target_v-mean))/(sd))
}

#回傳電阻值
resistors <- function(lowerbond,upperbond,n=50000,type="1%")
{
	R_data <- read.csv (paste(getwd(),"/電阻表.csv",sep=""), head=TRUE, sep=",")
	
	if(type == "0.5%")
	{
		R_0.5_data <- na.omit(R_data[R_data[,2]>=lowerbond & R_data[,2]<=upperbond,2])
		return (R_0.5_data)
	}
	else if(type == "1%")
	{
		R_1_data <- na.omit(R_data[R_data[,1]>=lowerbond & R_data[,1]<=upperbond,1])
		return (R_1_data)
	}
}

#畫常態分佈圖
plot_norm <- function ( mean,sd,target_v,bond) 
{
    curve(dnorm(x,mean=mean,sd=sd)
		 ,xlim = c(target_v-bond,target_v+bond)
		 ,ylab = ""
         ,xlab = paste("3sd = ",3*sd))

    title(main=
        substitute(
        paste("Normal Distribution: ",bar(x),"=",mean," , ","s","=",sd),
        list(mean=mean,sd=sd)))
	legend("topleft",paste("Out Spec: ",
    						signif(pnorm(target_v-bond,mean=mean,sd=sd)+1-pnorm(target_v+bond,mean=mean,sd=sd),digits=4)*100,
							"%",sep=""),inset=.01,bty="n")
 
	x_upper <- mean - 3*sd
	x_lower <- target_v - bond
	x <- c(x_lower, seq(x_lower, x_upper, (x_upper-x_lower)/100), x_upper)
	y <- c(0, dnorm(seq(x_lower, x_upper, (x_upper-x_lower)/100),mean=mean,sd=sd), 0)
	polygon(x, y, col = "chartreuse4")
	
	x_upper <- target_v + bond
	x_lower <- mean + 3*sd
	x <- c(x_lower, seq(x_lower, x_upper, (x_upper-x_lower)/100), x_upper)
	y <- c(0, dnorm(seq(x_lower, x_upper, (x_upper-x_lower)/100),mean=mean,sd=sd), 0)
	polygon(x, y, col = "chartreuse4")
}

#給電阻畫出常態分佈圖
plot_norm_byR <- function(target_v,bond,r1,r2,r3=0,n=20000)
{
	r1_n <- rnorm(n,r1,r1*0.01/3)
	r2_n <- rnorm(n,r2,r2*0.01/3)
	r3_n <- rnorm(n,r3,r3*0.01/3)
	v <- voltage(r1_n,r2_n,r3_n)
	plot_norm(voltage(r1,r2,r3),sd(v),target_v,bond)
	legend("topright",paste(c("R1","R2","R3","fit"),signif(c(r1,r2,r3,fitness(voltage(r1,r2,r3),sd(v),target_v,bond)),digits=3)),inset=.02)
}


#main
optimum_R <- function(target_v,bond,n)
{	
    fit <- -Inf
	
	#設定電阻
	type1 = "1%"
	factor1 = 0.01
	type2 = "1%"
	factor2 = 0.01
	sample1 <- resistors(1e3,1e6,n,type=type1)
	sample2 <- resistors(1e3,1e6,n,type=type2)
	L2 <- length(sample2)
	
	record <- matrix(,,ncol=6)
	colnames(c("R1","R2","mean","sd","fit","3*sd"))
	
	ideal_i <- which(abs(target_v-
						(voltage(rep(sample1,each=L2)
								,sample2)))
						<= bond)
	for(i in ideal_i)
	{	
		j <- ifelse((i%%L2)==0,L2,i%%L2)
		i <- ceiling(i/L2)
		#增加電阻要修改之
		ideal_v <- voltage(sample1[i]
						  ,sample2[j])
		
		#增加電阻要修改之
		v <- voltage(rnorm(n,sample1[i],sample1[i]*factor1/3)
					,rnorm(n,sample2[j],sample2[j]*factor2/3))

		sd <- sd(v)
		fit_t <- fitness(ideal_v,sd,target_v,bond)
		
		if(fit_t > fit)
		{
			fit <- fit_t
			#增加電阻要修改之
			assm <- c(sample1[i]
					 ,sample2[j])
			names(assm) <- c(paste("R1",type1,sep="_")
							,paste("R2",type2,sep="_"))

			info <- c(ideal_v,sd,fit,3*sd)
			names(info) <- c("mean","sd","fit","3*sd")
			cat("---------------------------------------------------\n")
			print(assm)
			print(info)
			cat("---------------------------------------------------\n")				
			record <- rbind(record,c(assm,info))					
		}
		else if(info[1] == ideal_v)
		{
			#增加電阻要修改之
			assm_t <- c(sample1[i]
					   ,sample2[j])
			info_t <- c(ideal_v,sd,fit_t,3*sd)
			record <- rbind(record,c(assm_t,info_t))
		}
	}
	#畫圖
	plot_norm(mean=voltage(as.numeric(assm[1]),as.numeric(assm[2])),sd=as.numeric(info[2]),target_v,bond)
	legend("topright",paste(names(c(assm,info[3])),signif(c(as.numeric(assm),info[3]),digits=3)),inset=.02)

	return (na.omit(record))
}

target_v=17.7;bond=0.3;n=20000;
system.time(record <- optimum_R(target_v,bond,n))
filename <- paste(target_v,bond,format(Sys.time(), "%m%d%Y.csv"),sep="-")
write.csv(record,file = paste(getwd(),filename,sep="/"))