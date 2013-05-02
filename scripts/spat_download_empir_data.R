
setwd('~/maxent/spat')

library(RCurl)

urls = read.csv('./public_data_urls.csv', colClasses='character')

download_data = function(urls, output_path) {
  for (i in seq_along(urls)) {
    temp = getURL(urls[i])
    dat = read.csv(textConnection(temp))
    namesplit = unlist(strsplit(urls[i], split='/'))
    filename = namesplit[length(namesplit)]
    write.csv(dat, file = file.path(path, filename), row.names=FALSE)
    print(paste(filename, ' was downloaded and written to ',
                 file.path(path, filename), sep=''))
  }
}

download_data(urls$url, './data/raw_data')


