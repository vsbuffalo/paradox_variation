# stapley_et_al_2017_data.r -- Stapley et al (2017) data
library(tidyverse)


source('../R/utilities.r')

# Stapley et al (2017) data, with my hand-collected body length data
ds_raw <- read_csv('./stapley_et_al_2017.csv') 

ds <- ds_raw %>% 
        # convert the genome size to Gb like the Corbett-Detig data
        mutate(genome_size = genome_size_Mb * 1e6 / 1e9)  %>%
        # select out certain columns
        select(species, species_ott, group, genome_size, hcn,  
               map_length, sexual_system, size,
               num_markers, parpath) %>%
        # fix parpath 
        mutate(parpath=case_when(parpath == 'n' ~ FALSE,
                                 parpath == 'y' ~ TRUE,
                                 TRUE ~ NA)) %>%
        # add metadata
        mutate(data_source='Stapley et al (2017)') %>%
        mutate(group=tolower(group)) %>% 
        # rename group
        mutate(group=case_when(
                               group == 'animals' ~ 'animal',
                               group == 'plants' ~ 'plant',
                               TRUE ~ group
                               )) %>%
       mutate(genus = map_chr(strsplit(species, ' '), 1)) %>%
       select(-species_ott)

# Stapley data for Apis cerana is wrong...
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076459
# the largest linkage block, not total map length reported
ds$map_length[ds$species == 'Apis cerana'] <- 3942.7

# # grt the taxonomial datas from taxadb
# ds <- ds %>% mutate(id_raw=map(genus, taxadb::get_ids)) %>%
#               mutate(id = map_int(id_raw, ~ as.integer(sub('ITIS:', '', .))))


# ## Lookup and merge in taxonomy 
# taxa_file <- 'stapley_et_al_2017_data_hier.Rdata'
# if (!file.exists(taxa_file)) {
#   hier <- map(ds$id, possibly(hierarchy_full, NA))
#   save(hier, file=taxa_file)
# } else {
#   load(taxa_file)
# }
# ds$taxons <- NA
# ds$taxons <- hier

# # extract out relevant rank levels
# ds <- ds %>% mutate(kingdom = map_chr(taxons, extract_rank), 
#      phylum = map_chr(taxons, extract_rank, rank=c('Phylum', 'Division')),
#      subphylum = map_chr(taxons, extract_rank, rank=c('Subphylum', 'Subdivision'), default=NA),
#      family = map_chr(taxons, extract_rank, rank='Family', default=NA),
#      order = map_chr(taxons, extract_rank, rank='Order', default=NA),
#      superphylum = map_chr(taxons, extract_rank, rank=c('Superphylum', 'Superdivision'), 
#                            default=NA),
#      class = map_chr(taxons, extract_rank, rank='Class', default=NA))
# save(ds, file='stapley_et_al_2017.Rdata')

# final dataframe
ds_flat <- ds %>% mutate(log10_size=log10(size), 
                         # convert to Morgans
                         map_length=map_length/100) 

### Merge in the range size data
#ranges_df <- read_tsv('gbif_ranges.tsv')

#ds_flat %>% left_join(ranges_df) %>%
#  #filter(species %in% ld$species) %>%
#  filter(kingdom == 'Animalia') %>%
#  group_by(family) %>% summarize_at(c('range', 'map_length', 'size'), mean, na.rm=TRUE) %>%
#  ggplot(aes(map_length, -log10(range / size))) + 
#  geom_point() + geom_smooth()


#ds_flat %>% left_join(ranges_df) %>%
#  filter(kingdom == 'Animalia', -log10(range / size) < -4) %>%
#  ggplot(aes(map_length, -log10(range / size))) + 
#  geom_point() + geom_text_repel(aes(map_length, -log10(range/ size), label=species))

#ds_flat %>% left_join(ranges_df) %>%
#  filter(kingdom == 'Animalia') %>%
#  ggplot(aes(map_length, -log10(range / size))) + 
#  geom_point() #+ geom_text_repel(aes(map_length, log10(range/ size), label=species))


write_tsv(ds_flat, 'stapley_et_al_2017_updated.tsv')


