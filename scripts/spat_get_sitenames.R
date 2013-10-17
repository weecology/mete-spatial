
print('Finding and exporting site names, ...')

file_names = dir('./data/raw_data/')

raw_names = c("FERP07data.csv",          "species.csv",
              "LFDP_Census4-Part1.csv",  "LFDP_Census4-Part2.csv",
              "LFDP_spp.txt",            "Oosting_Trees_1998.txt",
              "bci_census7.txt",         "cocoli.txt",
              "cocolisp.txt",            "crosstimbers1998.csv",
              "m04_graveyard.csv",       "m07_landsend.csv",
              "m12_rocky.csv",           "m13_bormann.csv",
              "m91_woodbridge.csv",      "m92_baldmnt.csv",
              "m93_bryan.csv",           "m94_bigoak.csv",
              "serpentine_data.csv",     "sherman.txt",
              "shermansp.txt")

shrtnames = c('ucsc', NA, 'luquillo', 'luquillo', 'luquillo',
              'oosting', 'bci', 'coocli', 'cocoli', 'cross',
              'graveyard', 'landsend', 'rocky', 'borman',
              'woodbridge', 'baldmnt', 'bryan', 'bigoak', 'serp',
              'sherman', 'sherman')

data_key = cbind(shrtnames, raw_names)

sitenames = unique(data_key[match(file_names, data_key[ , 2]), 1])

write.table(matrix(sitenames, nrow=1), file='./data/shrtnames.txt',
            sep=' ', row.names=FALSE, col.names=FALSE, quote=FALSE)

print('Finding and exporting site names, complete!')
