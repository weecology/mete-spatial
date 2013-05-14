
setwd('~/maxent/spat')

library(RCurl)

urls = read.csv('./public_data_urls.csv', colClasses='character')

download_data = function(urls, delim, output_path) {
  for (i in seq_along(urls)) {
    temp = getURL(urls[i])
    if (delim[i] == 'comma')
      dat = read.csv(textConnection(temp))
    if (delim[i] == 'tab')
      dat = read.delim(textConnection(temp))
    namesplit = unlist(strsplit(urls[i], split='/'))
    filename = namesplit[length(namesplit)]
    if (delim[i] == 'comma')
      write.csv(dat, file = file.path(output_path, filename), row.names=FALSE)
    if (delim[i] == 'tab')
      write.table(dat, file = file.path(output_path, filename), sep='\t',
                  row.names=FALSE)
    print(paste(filename, ' was downloaded and written to ',
                 file.path(output_path, filename), sep=''))
  }
}

print('Downloading datasets, ...')

download_data(urls$url, urls$delim, './data/raw_data')

print('Downloading datasets, complete!')
