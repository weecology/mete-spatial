
print('Downloading datasets, ...')

dir.create('./data/raw_data')

library(RCurl)

clArgs = commandArgs(trailingOnly=T)
if (length(clArgs) > 0)
  sites = clArgs
if (length(clArgs) == 0)
  sites = c('ucsc', 'luquillo', 'oosting')

urls = read.csv('./public_data_urls.csv', colClasses='character')
urls = urls[urls$sitename %in% sites, ]

download_data = function(urls, delim, output_path) {
  ## Function that downloads flat data files from the web
  ## Arguments:
  ## urls: the web address of the flat data file
  ## delim: the column delimintors for each dataset, either
  ##   'comma' or 'tab'
  ## output_path: the directory that the data will be written to
  require(RCurl)
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

download_data(urls$url, urls$delim, './data/raw_data')

print('Downloading datasets, complete!')
