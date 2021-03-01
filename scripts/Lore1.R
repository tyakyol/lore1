library(magrittr)
library(dplyr)
library(RVenn)
library(doBy)
library(stringr)
library(purrr)

# Load the data & data wrangling ----------------
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
df = read.delim('data/multiple_exons.csv',
                stringsAsFactors = FALSE,
                sep = ',')

dk = paste0('DK', str_pad(1:19, width = 2, pad = 0))

# Combine all the data
data = rbind(dk01, dk08, dk11, dk14, dk17)
rm(dk01, dk08, dk11, dk14, dk17)
colnames(data) = c('Line', 'Chr', 'Pos', 'Dir', 'Exon')
data = data[data$Exon != '.', ]
data = data[!grepl(pattern = ';', x = data$Exon), ]
data = rbind(data, df)
rm(df)
data$Genes = sapply(strsplit(x = data$Exon, split = '\\.'), function(x) x[[1]])

# Load gff3 file
gff = read.delim('data/genes.txt',
                 stringsAsFactors = FALSE,
                 header = FALSE)
gff$Genes = sapply(strsplit(x = gff$V9, split = '='), function(x) x[[2]])
colnames(gff)[4] = 'Start'
colnames(gff)[5] = 'End'

# Combine insertion with gff3
comb = left_join(data, gff, by = 'Genes')

# Calculate distance to the gene start
comb = comb %>% 
  mutate('Distance' = ifelse(V7 == '+', Pos - Start, End - Pos))

# Aggregate
agg_d = aggregate(Distance ~ Exon, data = comb, c)
agg_n = aggregate(Line ~ Exon, data = comb, c)
agg_c = aggregate(Chr ~ Exon, data = comb, c)
agg_p = aggregate(Pos ~ Exon, data = comb, c)
agg_s = aggregate(Dir ~ Exon, data = comb, c)

# Remove duplicates
f = lapply(agg_p$Pos, function(x) which(!duplicated(x)))
agg_d$Distance = map2(agg_d$Distance, f, function(x, y) x[y])
agg_n$Line = map2(agg_n$Line, f, function(x, y) x[y])
agg_c$Chr = map2(agg_c$Chr, f, function(x, y) x[y])
agg_p$Pos = map2(agg_p$Pos, f, function(x, y) x[y])
agg_s$Dir = map2(agg_s$Dir, f, function(x, y) x[y])

# Insertion selection ---------------------------
# Line
y = list()
for(i in seq_len(nrow(agg_d))) {
  if(length(agg_d$Distance[[i]]) > 2) {
    if(all(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])) {
      y[[i]] = agg_n$Line[[i]][which.minn(agg_d$Distance[[i]], 2)]
    } else if(sum(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19]) >= 2) {
      r = which(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])
      e = r[which.minn(agg_d$Distance[[i]][r], 2)]
      y[[i]] = agg_n$Line[[i]][e]
    } else {
      y[[i]] = agg_n$Line[[i]][which.minn(agg_d$Distance[[i]], 2)]
    }
  } else {
    y[[i]] = agg_n$Line[[i]]
  }
}
unite(Venn(y))
table(unlist(lapply(y), function(x) substr(x, 1, 4)))

# Chromosome
a = list()
for(i in seq_len(nrow(agg_d))) {
  if(length(agg_d$Distance[[i]]) > 2) {
    if(all(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])) {
      a[[i]] = agg_c$Chr[[i]][which.minn(agg_d$Distance[[i]], 2)]
    } else if(sum(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19]) >= 2) {
      r = which(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])
      e = r[which.minn(agg_d$Distance[[i]][r], 2)]
      a[[i]] = agg_c$Chr[[i]][e]
    } else {
      a[[i]] = agg_c$Chr[[i]][which.minn(agg_d$Distance[[i]], 2)]
    }
  } else {
    a[[i]] = agg_c$Chr[[i]]
  }
}

# Position
b = list()
for(i in seq_len(nrow(agg_d))) {
  if(length(agg_d$Distance[[i]]) > 2) {
    if(all(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])) {
      b[[i]] = agg_p$Pos[[i]][which.minn(agg_d$Distance[[i]], 2)]
    } else if(sum(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19]) >= 2) {
      r = which(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])
      e = r[which.minn(agg_d$Distance[[i]][r], 2)]
      b[[i]] = agg_p$Pos[[i]][e]
    } else {
      b[[i]] = agg_p$Pos[[i]][which.minn(agg_d$Distance[[i]], 2)]
    }
  } else {
    b[[i]] = agg_p$Pos[[i]]
  }
}

# Direction
s = list()
for(i in seq_len(nrow(agg_d))) {
  if(length(agg_d$Distance[[i]]) > 2) {
    if(all(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])) {
      s[[i]] = agg_s$Dir[[i]][which.minn(agg_d$Distance[[i]], 2)]
    } else if(sum(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19]) >= 2) {
      r = which(substr(agg_n$Line[[i]], 1, 4) %in% dk[8:19])
      e = r[which.minn(agg_d$Distance[[i]][r], 2)]
      s[[i]] = agg_s$Dir[[i]][e]
    } else {
      s[[i]] = agg_s$Dir[[i]][which.minn(agg_d$Distance[[i]], 2)]
    }
  } else {
    s[[i]] = agg_s$Dir[[i]]
  }
}

# Write out the data tables as gff3 -------------
whole = data[, c(2, 3, 1)]
whole$Pos2 = data$Pos
whole$V2 = rep('PGSB', nrow(whole))
whole$V3 = rep('SNP', nrow(whole))
whole$V6 = rep('.', nrow(whole))
whole$V7 = ifelse(data$Dir == 'F', '+', '-')
whole$V8 = rep('.', nrow(whole))
whole = whole[, c('Chr', 'V2', 'V3', 'Pos', 'Pos2', 'V6', 'V7', 'V8', 'Line')]
whole$Line = paste0('ID=', whole$Line)
write.table(whole, 'results/whole.gff3',
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = FALSE)

subd = data.frame(matrix(NA, nrow = length(flatten_chr(a)), ncol = 3))
subd[, 1] = flatten_chr(a)
subd[, 2] = flatten_int(b)
subd[, 3] = flatten_chr(y)
colnames(subd) = c('Chr', 'Pos', 'Line')
subd$Pos2 = subd$Pos
subd$V2 = rep('PGSB', nrow(subd))
subd$V3 = rep('SNP', nrow(subd))
subd$V6 = rep('.', nrow(subd))
subd$V7 = ifelse(flatten_chr(s) == 'F', '+', '-')
subd$V8 = rep('.', nrow(subd))
subd = subd[, c('Chr', 'V2', 'V3', 'Pos', 'Pos2', 'V6', 'V7', 'V8', 'Line')]
subd$Line = paste0('ID=', subd$Line)
write.table(subd, 'results/subdata.gff3',
            quote = FALSE,
            sep = '\t',
            row.names = FALSE,
            col.names = FALSE)
