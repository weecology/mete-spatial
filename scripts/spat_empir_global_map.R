library(maps)
library(rgdal)

setwd('~/GIS/country')
country = readOGR('country.shp','country')

setwd('~/maxent/spat')

dat = read.csv('./data/empir_geo_coords.csv')
head(dat)

apply(dat[,-1], 2, range)

map('world', xlim=c(-125, -64), ylim=c(7, 45))
map('state', add=TRUE)
#axis(side=1)
#axis(side=2)
points(dat[,3], dat[,2], col='dodgerblue', pch=19, cex=1.25)

plot(country, xlim=c(-125, -64), ylim=c(7, 45))
points(dat[,3], dat[,2], col='dodgerblue', pch=19, cex=1.25)

## blow up on NC
data(us.cities)
map('county', c('north carolina,orange', 'north carolina,durham'))
points(dat[,3], dat[,2], col='dodgerblue', pch=19, cex=1.25)
map.cities(us.cities, country="NC")

## blow up Panama
map('world', 'panama')
axis(side=1)
axis(side=2)
plot(country, xlim=c(-80, -79), ylim=c(7, 10))
points(dat[,3], dat[,2], col='dodgerblue', pch=19, cex=1.25)

## To do
## output an OGR layer for google earth this way we can get quick satelitte
## imagry