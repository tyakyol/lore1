# This script deals with the rows with multiple exons. It converts them to one per line.

library(purrr)

dk01 = read.delim('data/DK01-07_insertion_table_sorted.txt.exons.ids.txt',
                  stringsAsFactors = FALSE,
                  header = FALSE)
dk08 = read.delim('data/DK08-10_insertion_table_sorted.txt.exons.ids.txt',
                  stringsAsFactors = FALSE,
                  header = FALSE)
dk11 = read.delim('data/DK11-13_insertion_table_sorted.txt.exons.ids.txt',
                  stringsAsFactors = FALSE,
                  header = FALSE)
dk14 = read.delim('data/DK14-16_insertion_table_sorted.txt.exons.ids.txt',
                  stringsAsFactors = FALSE,
                  header = FALSE)
dk17 = read.delim('data/DK17-19_insertion_table_sorted.txt.exons.ids.txt',
                  stringsAsFactors = FALSE,
                  header = FALSE)

data = rbind(dk01, dk08, dk11, dk14, dk17)
rm(dk01, dk08, dk11, dk14, dk17)
colnames(data) = c('Line', 'Chr', 'Pos', 'Dir', 'Exon')
data = data[data$Exon != '.', ]

foo = vector(mode = 'list', length = nrow(data))
for(i in seq_len(nrow(data))) {
  if(grepl(pattern =';', x = data$Exon[i])) {
    x = unlist(strsplit(x = data$Exon[i], split = ';'))
    for(j in seq_len(length(x))) {
      foo[[i]][[j]] = unlist(list(data[i, 1], data[i, 2], data[i, 3], data[i, 4], x[j]))
    }
  } else {
    foo[[i]] = NA
  }
}
foo = foo[!is.na(foo)]

bar = vector(mode = 'list', length = length(foo))
for(i in seq_len(length(foo))) {
  bar[[i]] = lapply(foo[[i]], function(x) t(as.data.frame(x)))
}

baz = vector(mode = 'list', length = length(bar))
for(i in seq_len(length(bar))) {
  baz[[i]] = reduce(bar[[i]], function(x, y) rbind(x, y))
}

df = as.data.frame(reduce(baz, function(x, y) rbind(x, y)))
colnames(df) = c('Line', 'Chr', 'Pos', 'Dir', 'Exon')
write.csv(df, 'data/multiple_exons.csv', row.names = FALSE)
