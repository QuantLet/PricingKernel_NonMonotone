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

# Load packages
library(KernSmooth)
library(rugarch)
# library(fGarch)
library(forecast)
library(TSA)
library(ks)
library(matlab)
library(plotly)
library(htmlwidgets)
library(ggplot2)

#--------------------------------------------------------#
# (0) Adjust the following Parameters for different data #
#--------------------------------------------------------#

# path = "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/DAXdata/"
# # File Names 
# # 1.spot, 2.strike, 3.rate, 4.maturity, 5.oprice, 6.type, 7.ivola, 8.date
# Odata		= paste0(path,"odata20142212.txt")	# option data
# Idata		= paste0(path,"DAX.txt")		# dax index data
# 
# # Reading Option Data
# OdataAll 	= read.table(Odata, header=TRUE)
# Odata1		= subset(OdataAll,OdataAll[,7]<1 & OdataAll[,4]<0.9 
#                  & OdataAll[,6]==1)
# 
# # Reading Index Data
# dax1		= read.table(Idata, header=TRUE)[,2]
# periods		= length(dax1)
# dax		= c(1:periods)
# # Apply the foolowing loop if the data are given in reversed order
# # if not then set dax=dax1 
# dax		= flipud(dax1)
# dax.ts 		= ts(dax)
# dax.ret		= log(dax.ts)[2:periods]-log(dax.ts)[1:(periods-1)]
# dax.retts	= ts(dax.ret)

# BRC data
path = "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/Keynote_Figure/"
path = "/Users/rutingwang/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/Keynote_Figure/"

Odata		= paste0(path,"final_btc_option_2018-04-01-2018-06-30.csv")	# option data

Idata		= paste0(path,"BTC_USD_Quandl_2022.csv")		# index data
# 
OdataAll 	= read.csv(Odata, header=TRUE)

OdataAll <- readRDS(paste0(path,'final_btc_option_2018-04-01-2018-06-30.Rds'))

# Odata_sub = OdataAll[OdataAll$date == "2021-03-10",]
sub_date = "2018-04-01"
# sub_date ="2021-03-10"
Odata_sub = OdataAll[OdataAll$date == sub_date,]

Odata_sub$type = 0
Odata_sub$type[Odata_sub$option_type == "C"] = 1

Odata_sub$rate = 0
Odata_sub = Odata_sub[,c("indexprice","strike_price","interest_rate","maturity","markprice","type","iv","date")]
# 1.spot, 2.strike, 3.rate, 4.maturity, 5.oprice, 6.type, 7.ivola, 8.date
colnames(Odata_sub)[c(1,2,5,7)] = c("spot","strike","oprice","ivola")
Odata_sub$ivola = Odata_sub$ivola / 100

# Odata1		= subset(Odata_sub,Odata_sub[,4]>0 & Odata_sub[,6]==1)

Odata1		= subset(Odata_sub,Odata_sub[,4] >= 27/365& Odata_sub[,6]==1)

# Odata1		= subset(Odata_sub, Odata_sub[,6]==1)

# figure fo IV, call price, moneyness
# Example data
# moneyness <- mean(Odata1$spot,na.rm = TRUE)/Odata1$strike
moneyness <- Odata1$spot/Odata1$strike
IV<- Odata1$ivola
callprice <- Odata1$oprice
tau = Odata1$maturity

Odata1$Moneyness = Odata1$spot/Odata1$strike

Plot_Out = Odata1[Odata1$maturity == unique(Odata1$maturity),c('Moneyness','oprice')]

my_plot <-ggplot(Plot_Out, aes(x = Moneyness, y = oprice)) +
  geom_point(size = 2) +  # Scatter plot
  labs(x="Moneyness",y = "Option price")+
  theme(legend.key = element_rect(fill = "transparent"),
        axis.text = element_text(colour = 'black', size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(-0.2, "cm"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
        panel.background = element_rect(fill = "transparent"),
  )
my_plot
ggsave(paste0(path,"/tau_",unique(Odata1$maturity)*365,"_M_Op.png"), my_plot, width = 3200, height = 1600, units = "px")


my_plot <-ggplot(Plot_Out, aes(x = Moneyness, y = IV)) +
  geom_point(size = 2) +  # Scatter plot
  labs(x="Moneyness",y = "IV")+
  theme(legend.key = element_rect(fill = "transparent"),
        axis.text = element_text(colour = 'black', size = 18),
        axis.line = element_line(color = "black", linewidth = 1),
        axis.ticks.length = unit(-0.2, "cm"),
        panel.grid = element_blank(),
        plot.background = element_rect(fill = "transparent"),
        legend.background = element_rect(fill = "transparent"),  # Set legend background to transparent
        panel.background = element_rect(fill = "transparent"),
  )
my_plot
ggsave(paste0(path,"/tau_",unique(Odata1$maturity)*365,"_M_IV.png"), my_plot, width = 3200, height = 1600, units = "px")
