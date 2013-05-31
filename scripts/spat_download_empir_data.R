setwd('~/maxent/spat')

dir.create('./data/raw_data')

library(RCurl)

source('./scripts/spat_functions.R')

urls = read.csv('./public_data_urls.csv', colClasses='character')

print('Downloading datasets, ...')

download_data(urls$url, urls$delim, './data/raw_data')

print('Downloading datasets, complete!')
