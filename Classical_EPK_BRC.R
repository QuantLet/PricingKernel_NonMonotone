# Clear variables and close windows
rm(list = ls(all = TRUE))
graphics.off()

# Set working directory
# setwd('~/...')    		# linux/mac os
# setwd('/Users/...') 		# windows

## Install packages
# install.packages(KernSmooth)
# install.packages(rugarch)
# install.packages(fGarch)
# install.packages(forecast)
# install.packages(TSA)
# install.packages(ks)
# install.packages(matlab)
# install.packages("optimx")

# Load packages
library(KernSmooth)
library(rugarch)
# library(fGarch)
library(forecast)
library(TSA)
library(ks)
library(matlab)
library(splines) # bs()
library(pracma)
library(optimx)
library(data.table) 
library(dplyr)
library(ggplot2)


# set directory
# path <- "~/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/PAPER/Pricing/code/test/"
path <- "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/main_code/"
setwd(path)


#--------------------------------------------------------#
# (0) data prepare, Q_file prepare
#--------------------------------------------------------#

# option <- readRDS(paste0(path,'20172022_processed_1_3_4.Rds'))
# Idata		= paste0(path,"BTC_USD_Quandl_2022.csv")		# index data
Idata		= paste0(path,"BTC_USD_yf_2017-2022.csv")		# index data
BTC_price		=  read.csv(Idata, header=TRUE)

# Qfile of zijin
path_Qfile <- paste0(path,"Q_density_files_new/")
directories <- list.files(path = path_Qfile, full.names = FALSE, recursive = FALSE)


# classical part

HDEgarch 		= function(garchfit, data, grid, maturity, N, start){
  
  # data = sp1.retts
  # grid = data_q$Strike
  # maturity = tau
  # N = 5000
  # start = SpotPrice 
  
  # Simulation
  garchsim 		= ugarchsim(garchfit, n.sim = round(maturity), 
                         n.start = 0, m.sim=N, startMethod=("sample"), 
                         mexsimdata=TRUE)

  
  # head(g11sim@simulation)[2]= extracting simulated return data
  returnsim		= as.vector(head(garchsim@simulation)[2])
  
  # Calculating spot prices from return data
  value			= matrix(0,ncol=N)
  for(i in 1:N) {
    value[,i]  	= start*exp(sum(returnsim$seriesSim[,i]))		
  }
  
  # Computing density on given grid: Either by using build in function (kde(...)):
  # kde() computes bandwith h wih function hpi() which uses Wand&Jones (1994) estimator
  # grid: strike
  ValDens			= kde(value[1,],eval.points=grid, gridsize=length(grid),
                  xmin=min(grid), xmax=max(grid))
  return(ValDens$estimate)
}


sp1 <- BTC_price
sp1.retts = ts(log(sp1$Adj.Close)[2:nrow(sp1)]-log(sp1$Adj.Close)[1:(nrow(sp1)-1)])

# Parameters for calculation/simulation
numbapprox  	= 2000			# fineness of the grid
N		= 5000			# No. of Simulations
# Check return series for ARMA effects, e.g. with the following function
# auto.arima(dax.retts, max.p=10, max.q=10, max.P=5, max.Q=5, 
# start.p=1, start.q=1,start.P=1, start.Q=1, stationary=T, seasonal=F)
p		= 0
q		= 0
arma		= c(p,q)
# specify garch order (need to be checked)
m		= 1
s		= 1
garch		= c(m,s)
# Specify GARCH model (default is standard GARCH)
# for changing GARCH-model + submodel, please refer to 
# rugarch package for further information
garchmodel	= "fGARCH"
submodel	= "GARCH"
# underlying distribution (default: "sstd" - skewed stundent t's)
# (alternatives: "norm" - normal, "ghyp"- generalized hyperbolic)
udist		= "sstd"
# set archm=T for ARCH in mean model (archpow specifies the power)
archm		= F
archpow		= 1
# set include.mean = F if you don't want to include a mean in the mean model
include.mean 	= T  

spec			= ugarchspec(variance.model = list(model = garchmodel, 
                                          garchOrder = garch, submodel = submodel), mean.model = 
                      list(armaOrder = arma, archm=archm,archpow=archpow,
                           include.mean=include.mean), distribution.model = udist)



# 4y 6y 8y 10y
# rolling = c(1, 2, 3, 4)
rolling = 4
# tau = 30
# tau_files <- grep( paste0("raw_Q_density.*tau",tau,"\\.csv$"), directories, value = TRUE)
# dates <- gsub("raw_Q_density_([0-9]{4}-[0-9]{2}-[0-9]{2})_tau.*", "\\1", tau_files)
# dates_Q <- dates[1:(length(dates) - 1)]

# date_tau = data.frame(dates_Q = c("2022-12-27","2022-10-04","2022-10-28"),tau = c(17,24,28))
# date_tau = data.frame(dates_Q = c("2022-11-17","2022-10-29","2022-07-22","2022-08-31"),tau = c(8,13,7,30))
date_tau = data.frame(dates_Q = c("2022-11-17","2022-10-29"),tau = c(8,13))

for (i in c(1:nrow(date_tau))){
  setwd(paste0(path, 'Q_density_files_new'))
  tau = date_tau$tau[i]
  figure_save = paste0(path,"EPKPlot/tau_",tau)
  dir.create(figure_save, showWarnings = FALSE)
 
  Q_band <- paste0(path_Qfile,"raw_Q_density_", date_tau$dates_Q[i], "_tau", date_tau$tau[i], "_band.csv")
  # 检查文件是否存在
  if (!file.exists(Q_band)) {
    print(paste("no files for", date_tau$dates_Q[i], "tau", date_tau$tau[i]))
    next  # 如果文件不存在，跳过这个循环
  }
  
  data_q = read.csv(Q_band)
  SpotPrice = mean(data_q$Strike * data_q$m)
  
  for (iRoll in rolling){
    sp1_h <- sp1[as.Date(sp1$Date) < as.Date(date_tau$dates_Q[i]), ]
    sp1_h <- sp1[max(1,(nrow(sp1_h)-iRoll*360)):nrow(sp1_h),]
    sp1.retts = ts(log(sp1_h$Adj.Close)[2:nrow(sp1_h)]-log(sp1_h$Adj.Close)[1:(nrow(sp1_h)-1)])
    
    garchfit 		= ugarchfit(data=sp1.retts, spec=spec, solver = "hybrid")
    
    HDE = HDEgarch(garchfit, sp1.retts, data_q$Strike, 
                   tau, N = 5000, SpotPrice)
    
    EPK = data.frame(M = data_q$m,EPK = data_q$SPD  / HDE, EPK_lo = data_q$SPD_lo/ HDE, EPK_up = data_q$SPD_up/ HDE, 
                     Q = data_q$SPD, Q_lo = data_q$SPD_lo,Q_up = data_q$SPD_up, P = HDE,Strike = data_q$Strike ) 
    
    EPK = EPK[EPK$EPK>0&EPK$M<1.07&EPK$M>0.8,]
    
    # Save the plot to a PNG file
    png(paste0(figure_save, "/BTC_", date_tau$dates_Q[i], "_tau_", tau, "_", iRoll, "YearRoll_Classical.png"), width = 1000, height = 800, bg = "transparent")
    
    # Ensure there are no missing or infinite values in EPK$M and EPK$EPK
    EPK <- EPK[!is.na(EPK$M) & !is.na(EPK$EPK) & !is.infinite(EPK$M) & !is.infinite(EPK$EPK), ]
    
    
    # Check if the data frame is not empty after removing NA and Inf values
    if (nrow(EPK) > 0) {
      # Calculate xlim and ylim ensuring finite values
      xlim_vals <- range(EPK$M, na.rm = TRUE, finite = TRUE)
      ylim_vals <- quantile(c(EPK$EPK_lo,EPK$EPK,EPK$EPK_up), probs = c(0.05, 0.95), na.rm = TRUE, finite = TRUE)
      
      # Debugging: Print the calculated xlim and ylim values
      print(paste("xlim:", xlim_vals))
      print(paste("ylim:", ylim_vals))
      
      # Plot only if xlim and ylim values are finite
      if (all(is.finite(xlim_vals)) && all(is.finite(ylim_vals))) {
        par(mar = c(5, 6, 4, 2) + 0.1)  # 增大左侧和底部边距
        plot(EPK$M, EPK$EPK, type = 'l', col = 'black', lwd = 6,
             xlab = "Moneyness", ylab = "EPK",
             cex.lab = 2,                # 坐标轴标签字体大小
             cex.axis = 2,
             xlim = xlim_vals, ylim = ylim_vals,   # x 和 y 轴的范围
             cex.main = 2)
        
        lines(EPK$M, EPK$EPK_lo, type = 'l', col = 'blue', lwd = 6, lty = 2)
        
        lines(EPK$M, EPK$EPK_up, type = 'l', col = 'blue', lwd = 6, lty = 2)
        
      } else {
        warning("Non-finite xlim or ylim values, plot not created.")
      }
    } else {
      warning("Data frame is empty after removing NA and Inf values, plot not created.")
    }
    
    # Close the graphics device
    dev.off()
    
    # plot(EPK$Strike, EPK$Q, type = 'l', col = 'red', 
    #      xlab = "Strike", ylab = "Q density",
    #      xlim = c(min(EPK$Strike), max(EPK$Strike)), ylim = c(quantile(EPK$Q, probs = c(0.01)),quantile(EPK$Q, probs = c(0.99))),
    #      main = paste(dates_Q[i], "Q density")
    # )
    # lines(EPK$Strike, EPK$P, col = 'blue')
    # 
    # 
    # plot(EPK$Strike, EPK$P, type = 'l', col = 'black', 
    #      xlab = "Strike", ylab = "P density",
    #      xlim = c(min(EPK$Strike), max(EPK$Strike)), ylim = c(quantile(EPK$P, probs = c(0.01)),quantile(EPK$P, probs = c(0.99))),
    #      main = paste(dates_Q[i], "Classical EPK")
    # )
    # 
    Plot_Out = data.frame(
      x = rep(EPK$M, 4),
      density = c(EPK$Q, EPK$Q_lo, EPK$Q_up, EPK$P),
      line = factor(rep(c("Q_density", "Q_density_Lowbound", "Q_density_Upbound", "P_density"), each = nrow(EPK)))
    )
    my_plot <-ggplot(Plot_Out, aes(x = x, y = density, color = line)) +
      geom_line(size = 1) +  
      labs(x=NULL,y = "Density")+
      theme(legend.key = element_rect(fill = "transparent"),
            axis.ticks.length = unit(-0.2, "cm"),
            panel.grid = element_blank(),
            axis.line = element_line(color = "black", linewidth = 0.5),
            plot.background = element_rect(fill = "transparent"),
            legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
            panel.background = element_rect(fill = "transparent"),
      )
    plot(my_plot)
    ggsave(paste0(figure_save,"/BTC_Density_",date_tau$dates_Q[i],"_tau_",tau,"_",iRoll,"YearRoll.png"), my_plot, width = 3200, height = 1600, units = "px")
    
  }
  
}

# compare Rookley Q with Fel Q

figure_save = paste0(path,"Q_compare/20220722_tau_7")
dir.create(paste0(path,"Q_compare"))
dir.create(figure_save)
# Rookley Q
data_q =  read.csv(paste0(path_Qfile,"raw_Q_density_2022-07-22_tau7.csv"))
data_q = data_q[order(data_q$m),]
data_q$Q_density = data_q$y


xlim_vals <- range(data_q$m, na.rm = TRUE, finite = TRUE)
ylim_vals <- quantile(data_q$Q_density, probs = c(0.05, 0.95), na.rm = TRUE, finite = TRUE)

# Save the plot to a PNG file
png(paste0(figure_save, "/BTC_Rookley_Q.png"), width = 800, height = 600, bg = "transparent")

plot(data_q$m, data_q$Q_density, type = 'l', col = 'black', lwd = 2,
     xlab = "Moneyness", ylab = "Q density",
     cex.lab = 1.25,                # 坐标轴标签字体大小
     cex.axis = 1.25,
     xlim = xlim_vals, ylim = ylim_vals,   # x 和 y 轴的范围
     cex.main = 1.25)

dev.off()


# Figlewski Q

data_q_F =  read.csv('Compare_Q_Figlewski/BRC/tau7/RND_Figlewski_7days_2022-07-22.csv')
data_q_F = data_q_F[order(data_q_F$Moneyness),]
data_q_F$Q_density = data_q_F$RND_K
data_q_F = data_q_F[data_q_F$Moneyness<=1.2 & data_q_F$Moneyness>=0.8,]

xlim_vals <- range(data_q_F$Moneyness, na.rm = TRUE, finite = TRUE)
ylim_vals <- quantile(data_q_F$Q_density, probs = c(0.05, 0.95), na.rm = TRUE, finite = TRUE)
png(paste0(figure_save, "/BTC_Figlewski_Q.png"), width = 800, height = 600, bg = "transparent")

plot(data_q_F$Moneyness, data_q_F$Q_density, type = 'l', col = 'black', lwd = 2,
     xlab = "Moneyness", ylab = "Q density",
     cex.lab = 1.25,                # 坐标轴标签字体大小
     cex.axis = 1.25,
     xlim = xlim_vals, ylim = ylim_vals,   # x 和 y 轴的范围
     cex.main = 1.25)
dev.off()




