#需設定電壓公式，電阻％數與電阻範圍跟電阻數目
#輸入理想電壓，差異值，樣本數

options(keep.source = FALSE)
#計算電壓值
voltage <- function(r1,r2,r3=0,r4=NULL)
{
    ref <- 0.925
    if(length(r1) > 1)
    {
        ref <- rnorm(length(r1),ref,0.009/3)
    }
    if(length(r4)==0)
    {
        return (ref*(1+r1/(r2+r3)))#更改成正確電壓公式
    }
    else
    {
        return (ref*(1+r1/(r2+r3+r4))) #更改成正確電壓公式
    }
}


#計算合適值
fitness <- function(mean,sd,target_v,bond)
{
    return ((bond-abs(target_v-mean))/sd)
}

#回傳電阻值
resistors <- function(lowerbond,upperbond,n=50000,type="1%")
{
	R_data <- read.csv (paste(getwd(),"/電阻表.csv",sep=""), head=TRUE, sep=",")
	
	if(type == "0.5%")
	{
		R_0.5_data <- na.omit(R_data[R_data[,2]>=lowerbond & R_data[,2]<=upperbond,2])
		R_0.5 <- matrix(data=NA,nrow=n,ncol=length(R_0.5_data))
		for(i in 1:length(R_0.5_data))
		{
			R_0.5[,i] <- rnorm(n,R_0.5_data[i],R_0.5_data[i]*0.005/3)
		}
		colnames(R_0.5) <- R_0.5_data
		return (R_0.5)
	}
	else if(type == "1%")
	{
		R_1_data <- R_data[R_data[,1]>=lowerbond & R_data[,1]<=upperbond,1]
		R_1 <- matrix(data=NA,nrow=n,ncol=length(R_1_data))
		for(i in 1:length(R_1_data))
		{
			R_1[,i] <- rnorm(n,R_1_data[i],R_1_data[i]*0.01/3)
		}
		colnames(R_1) <- R_1_data
		return (R_1)
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
plot_norm_byR <- function(target_v,bond,r1,r2,r3=c(0,0),r4=c(0,0),n=20000)
{
	r1_n <- rnorm(n,r1[1],r1[1]*r1[2]/3)
	r2_n <- rnorm(n,r2[1],r2[1]*r2[2]/3)
	r3_n <- rnorm(n,r3[1],r3[1]*r3[2]/3)
    r4_n <- ifelse(r4[1]==0,0,rnorm(n,r4[1],r4[1]*r4[2]/3))
	v <- voltage(r1_n,r2_n,r3_n,r4_n)
	plot_norm(voltage(r1[1],r2[1],r3[1],r4[1]),sd(v),target_v,bond)
	legend("topright",paste(c("R1","R2","R3","fit"),signif(c(r1[1],r2[1],r3[1],fitness(voltage(r1[1],r2[1],r3[1]),sd(v),target_v,bond)),digits=3)),inset=.02)
}


#main
optimum_R <- function(target_v,bond,n,mode=0)
{	
    fit <- -Inf

	if(mode == 1)
	{
		#設定電阻
		type1 = "1%"
		type2 = "1%"
		type3 = "1%"
		sample1 <- resistors(1e3,1e6,n,type=type1)
		sample2 <- resistors(1e3,1e6,n,type=type2)
		sample3 <- resistors(1e3,1e6,n,type=type3)
		l3 = length(sample3[1,])
		
		record <- matrix(,,ncol=7)
		colnames(c("R1","R2","R3","mean","sd","fit","3*sd"))
		
		for(i in 1:length(sample1[1,]))
		{
			ideal_j <- which(abs(target_v-
    							(voltage(as.numeric(colnames(sample1)[i])
									,rep(as.numeric(colnames(sample2)),each=l3)
										,as.numeric(colnames(sample3)   ))))
								<= bond)

			for(j in ideal_j)
			{
				k <- ifelse((j%%l3)==0,l3,j%%l3)
				j <- ceiling(j/l3)

				#增加電阻要修改之
				ideal_v = voltage(as.numeric(colnames(sample1)[i])
								 ,as.numeric(colnames(sample2)[j])
								 ,as.numeric(colnames(sample3)[k]))
		
				#增加電阻要修改之
				v <- voltage(sample1[,i]
							,sample2[,j]
							,sample3[,k])

				sd = sd(v)
				fit_t <- fitness(ideal_v,sd,target_v,bond)

				if(fit_t > fit)
				{
					fit <- fit_t
					#增加電阻要修改之
					assm <- c(colnames(sample1)[i]
							 ,colnames(sample2)[j]
							 ,colnames(sample3)[k])
					names(assm) <- c(paste("R1",type1,sep="_")
									,paste("R2",type2,sep="_")
									,paste("R3",type3,sep="_"))
  
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
					assm_t <- c(colnames(sample1)[i]
							   ,colnames(sample2)[j]
							   ,colnames(sample3)[k])
					info_t <- c(ideal_v,sd,fit_t,3*sd)
					record <- rbind(record,c(assm_t,info_t))
				}
			}
		}
		cat("---------------------------------------------------\n")
		cat("Size\n")
		print(length(na.omit(record)))
		cat("---------------------------------------------------\n")
  
		#畫圖
		plot_norm(voltage(as.numeric(assm[1]),as.numeric(assm[2]),as.numeric(assm[3])),sd=as.numeric(info[2]),target_v,bond)
		legend("topright",paste(names(c(assm,info[3])),signif(c(as.numeric(assm),info[3]),digits=3)),inset=.02)
	}
	else
	{
		#設定電阻
		type1 = "1%"
		type2 = "1%"
		type3 = "1%"
		type4 = "1%"
 
 		sample1 <- resistors(1e0,1e6,n,type=type1)
 		sample2 <- resistors(1e0,1e6,n,type=type2)
		sample3 <- resistors(1e0,1e6,n,type=type3)
		sample4 <- resistors(1e0,1e6,n,type=type4)
		L1 <- length(sample1[1,])
		L2 <- length(sample2[1,])
		L3 <- length(sample3[1,])
		L4 <- length(sample4[1,])
		
		record <- matrix(,,ncol=8)
		colnames(c("R1","R2","R3","R4","mean","sd","fit","3*sd"))
		
		for(i in 1:length(sample1[1,]))
		{
			for(j in 1:length(sample2[1,]))
			{
				#增加電阻要修改之
				r1 = as.numeric(colnames(sample1)[i])
				r2 = as.numeric(colnames(sample2)[j])
				r3 = as.numeric(colnames(sample3))
				r4 = rep(as.numeric(colnames(sample4)),each = L3)
				ideal_t <- which(abs(target_v - voltage(r1,r2,r3,r4)) <= bond)

				for(index in ideal_t)
				{
					k <- (((index-1)%%L3) + 1)
					l <- (((ceiling(index/L3)-1)%%L4) + 1)
		
					#增加電阻要修改之
					ideal_v = voltage(as.numeric(colnames(sample1)[i])
									 ,as.numeric(colnames(sample2)[j])
									 ,as.numeric(colnames(sample3)[k])
									 ,as.numeric(colnames(sample4)[l]))
			
					#增加電阻要修改之
					v <- voltage(sample1[,i]
								,sample2[,j]
								,sample3[,k]
								,sample4[,l])

					sd = sd(v)
					fit_t <- fitness(ideal_v,sd,target_v,bond)

					if(fit_t > fit)
					{
						fit <- fit_t
						#增加電阻要修改之
						assm <- c(colnames(sample1)[i]
								 ,colnames(sample2)[j]
								 ,colnames(sample3)[k]
								 ,colnames(sample4)[l])
						names(assm) <- c(paste("R1",type1,sep="_")
										,paste("R2",type2,sep="_")
										,paste("R3",type3,sep="_")
										,paste("R4",type4,sep="_"))
	  
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
						assm_t <- c(colnames(sample1)[i]
								   ,colnames(sample2)[j]
								   ,colnames(sample3)[k]
								   ,colnames(sample3)[l])
						info_t <- c(ideal_v,sd,fit_t,3*sd)
						record <- rbind(record,c(assm_t,info_t))
					}
				}
			}
			#畫圖
			cat("---------------------------------------------------\n")
			cat("Size\n")
			print(length(na.omit(record)))
			cat("---------------------------------------------------\n")
			
			plot_norm(mean=voltage(as.numeric(assm[1]),as.numeric(assm[2]),as.numeric(assm[3]),as.numeric(assm[4])),sd=as.numeric(info[2]),target_v,bond)
			legend("topright",paste(names(c(assm,info[3])),signif(c(as.numeric(assm),info[3]),digits=3)),inset=.02)
		}
	}
	return (na.omit(record))
}

target_v=0.9;bond=0.02;n=10000;
mode=1;# 1 表示算 3 個電阻, 0 表示算 4 個電阻
system.time(record <- optimum_R(target_v,bond,n,mode))
filename <- paste(target_v,bond,format(Sys.time(), "%m%d%Y.csv"),sep="-")
write.csv(record,file = paste(getwd(),filename,sep="/"))