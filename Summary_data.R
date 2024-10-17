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
library(dplyr)
library(readxl)
# install.packages("reticulate")
library(reticulate)
library(xlsx)

# set directory
# path <- "~/Library/CloudStorage/OneDrive-NationalUniversityofSingapore/PAPER/Pricing/code/test/"
Path_raw =  "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/"
path <- "/Users/ruting/Library/Mobile Documents/com~apple~CloudDocs/PK_BTC/code_data/main_code/"
setwd(path)

raw_data_sp = read.csv('SP500Raw_data.csv')
raw_data_sp$strike = raw_data_sp$strike/1000
raw_data_sp$moneyness = raw_data_sp$index_price/raw_data_sp$strike

raw_data_cry = read.csv('CryptoRaw_data.csv')
raw_data_cry$date = as.POSIXct(raw_data_cry$date, tz = "UTC")
raw_data_cry$date = format(raw_data_cry$date, format = "%Y-%m-%d")

# SP 500 data
# date, tau, original obs, type, IV is not empty, moneyness between 0.8 and 1.2, 
date_tau =  read_excel("Sum_date.xlsx") %>% data.frame()
date_tau$dates_Q =  format(date_tau$dates_Q, format = "%Y-%m-%d")
date_tau$obs_raw = NA
date_tau$IV_Empty = NA
date_tau$IV_Zero= NA
date_tau$Moneyness_inRange = NA

for (i in c(1:nrow(date_tau))){
  Type = date_tau$Type[i]
  date = date_tau$dates_Q[i]
  tau = date_tau$tau[i]
  
  if (Type == "SP500"){
    temp_data = raw_data_sp[raw_data_sp$date == date & raw_data_sp$tau == tau,]
  }else{ 
    temp_data = raw_data_cry[raw_data_cry$date == date & raw_data_cry$tau == tau,]
    
  }

  date_tau$obs_raw[i] = nrow(temp_data)
  date_tau$IV_Empty[i] = length(which(is.na(temp_data$mark_iv)))
  date_tau$IV_Zero[i] = length(which(temp_data$mark_iv==0))
  date_tau$Moneyness_inRange[i] = length(which(temp_data$moneyness >= 0.8 & temp_data$moneyness <1.2
                                               & temp_data$mark_iv>0))
}
write.xlsx(date_tau,file = 'Data_summary.xlsx')

