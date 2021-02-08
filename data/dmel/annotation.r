library(GenomicRanges)
library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
library(tidyverse)

keep_seqs <- c("chr2L", "chr2R", "chr3R", "chr3L")

tx <- TxDb.Dmelanogaster.UCSC.dm3.ensGene
reduce <- IRanges:::reduce

sq_filters <- function(x, chroms) x[seqnames(x) %in% chroms]

dm_exs <- exons(tx, filter=list(tx_chrom=keep_seqs))
dm_introns <- unlist(sq_filters(intronsByTranscript(tx), keep_seqs))
dm_5UTRs <- unlist(sq_filters(fiveUTRsByTranscript(tx), keep_seqs))
dm_3UTRs <- unlist(sq_filters(threeUTRsByTranscript(tx), keep_seqs))

seqlens <- seqlengths(tx)
seqlens <- seqlens[names(seqlens) %in% keep_seqs]
G <- sum(seqlens)

d <- tribble(~ feature, ~ total_bp, 
                 'exon', sum(width(reduce(dm_exs, ignore.strand=TRUE))),
                 'intron', sum(width(reduce(dm_introns, ignore.strand=TRUE))),
                 'utr', sum(width(reduce(dm_5UTRs))) + 
                        sum(width(reduce(dm_3UTRs)))) 

# compute non-coding (drop exon; we use CDS)
noncoding <- G - d %>% pull(total_bp) %>% sum()

d <- d %>% add_row(feature='noncoding', total_bp=noncoding) %>% 
  mutate(prop_genome = total_bp / G) 

# sanity check --- all the features 1?
stopifnot(sum(d$prop_genome) == 1)

write_tsv(d, 'dmel_annotation_summary.tsv')

sl <- tibble(chrom = names(seqlens), length=seqlens)
write_tsv(sl, 'dmel_seqlens.tsv')




# col_types <- "c_ciic___"
# col_names <- c('seqname', 'source', 'feature', 'start', 'end', 'score', 'filter', 'strand', 'group', 'attribute')

# gff <- read_tsv('./dmel-all-r5.33.gff.gz', comment="#", col_types=col_types, col_names=col_names)

# seqlens <- read_tsv('./dmel-5.33-seqlens.txt', col_names=c("seqname", "length"))

# keep_features <- c('exon', 'intron', 'gene', 'five_prime_UTR', 'three_prime_UTR')
# keep_seqs <- c("2L", "2R", "3R", "3L")

# # only keep autosomes
# seqlens2 <- seqlens %>% filter(seqname %in% keep_seqs)

# # keep only key features 
# features <- gff %>% filter(feature %in% keep_features)

# # create a GRanges object
# gr <- with(features, GRanges(seqnames = seqname, 
#                              ranges = IRanges(start, end), feature=feature))

# # collapse overlapping ranges, within a feature type
# out <- lapply(split(gr, seqnames(gr)), function(y) {
#        lapply(split(y, y$feature), function(x) sum(width(IRanges:::reduce(x))))
# })

# # total basepairs
# d <- as.data.frame(do.call(rbind, out)) %>% rownames_to_column('seqname') %>%
#   unnest(c(exon, five_prime_UTR, gene, intron, three_prime_UTR))

# # get the lengths of all features on subset of chroms
# d2 <- d %>% filter(seqname %in% keep_seqs) %>% select(-seqname) %>% 
#          summarize_all(sum)

# # get number of features
# dfeats <- tibble(seqname=as.vector(seqnames(gr)), 
#                  start=start(gr),
#                  end=end(gr),
#                  feature=gr$feature,
#                  width=width(gr)) 

# dfeats2 <- dfeats %>% filter(seqname %in% keep_seqs) %>%  ungroup() %>%
#             group_by(feature) %>% summarize(n = n(), length=mean(width))
                 
# # dfeats %>% group_by(seqnames, feature) %>% summarize(n = n()) 
# # feats %>% group_by(feature) %>% summarize(n = n(), length=mean(width))

# # merging
# tbps <- tibble(total_bps = as.vector(t(d2)), feature=names(d2))

# x <- dfeats2 %>% left_join(tbps) %>% 
#       mutate(genome_length = sum(seqlens2$length)) %>%
#       mutate(prop = total_bps / genome_length)
    

# write_tsv(x, 'dmel_annotation_summary.tsv')
