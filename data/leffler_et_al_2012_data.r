## leffler_et_al_2012_data.r -- Leffler et al (2012) data
library(tidyverse)
library(ritis)
library(ggrepel)
source('../R/utilities.r')

# Load Leffler et al (2012) diversity data, with 
# my hand-collected body length data
ld <- read_tsv('leffler_et_al_2012.tsv')
colnames(ld) <- c('species', 'common_name', 'phylum', 'class', 'total_loci',
                  'mean_diversity_pops', 'thetaW_pi', 'chr_type', 'site_type', 
                  'sampling', 'reference', 'mating_system', 'habitat', 
                  'range', 'included_study', 'included_sex_auto', 'size', 
                  'size_source', 'size_comments')

# # I did not enter body size data in the Leffler data that exists elsewhere; 
# # we merge that in now
# dsc <- read_tsv('combined_data.tsv')
# size_lookup <- as.list(c(dsc$size[!is.na(dsc$size)], ld$size[!is.na(ld$size)]))
# names(size_lookup) <- c(dsc$species[!is.na(dsc$size)], ld$species[!is.na(ld$size)])

# # before lookup missing size data
# table(is.na(ld$size))
# tmp <- sapply(size_lookup[ld$species], function(x) ifelse(is.null(x), NA, x))
# ld$size <- as.numeric(tmp)
# # after merge
# table(is.na(ld$size))


## A few subspecies with different names *should* be averaged
# so we change their names here
# Combine intestinalis 
ciona <- regexpr("Ciona intestinalis", ld$species) != -1
ld$species[ciona] <- "Ciona intestinalis"
# subtypes A & B are listed as cosmopolitan and broad endemic. I average
# and call ie cosmopolitan.
ld$range[ciona] <- "Cosmopolitan"

# Mus subspecies are averaged since the GBIF data doesn't seem to distinguish 
# between subspecies
ld[regexpr('Mus musculus', ld$species) != -1, ]$species <- 'Mus musculus'


# average the diversity values, looking only at π, autosomes included in the study
lds <- ld %>% mutate(diversity = strsplit(mean_diversity_pops, ';',)) %>%
  # diversity is in percent; thus we divide by 100
  unnest(diversity)  %>% mutate(diversity=as.numeric(diversity)/100) %>%
  select(-mean_diversity_pops) %>% 
  filter(chr_type == 'Autosome', included_study == 'Y', thetaW_pi == 'Pi') %>%
  # average over pop measures
  group_by(species, phylum, size, site_type) %>% 
  summarize(diversity=mean(diversity), 
            habitat=unique(habitat),
            mating_system=unique(mating_system),
            range_cat=unique(range)) %>%
  mutate(genus = map_chr(species, ~ strsplit(., ' ')[[1]][1]))

# We don't have drosophila size data for a lot of species. 
# However, within this genus, there isn't a ton of size variation, so we just
# impute using the average
ave_drosoph_size <- mean(lds$size[lds$genus == 'Drosophila'], na.rm=TRUE)
lds$size[lds$genus == 'Drosophila' & is.na(lds$size)] <- ave_drosoph_size


# # get taxonomical IDs
# lds <- lds %>% mutate(id_raw=map(genus, taxadb::get_ids)) %>%
#                mutate(id = map_int(id_raw, ~ as.integer(sub('ITIS:', '', .))))

# ## Lookup and merge in taxonomy 
# lds_taxa_file <- 'leffler_et_al_2012_hier.Rdata'
# if (!file.exists(lds_taxa_file)) {
#   hier <- map(lds$id, possibly(hierarchy_full, NA))
#   save(hier, file=lds_taxa_file)
# } else {
#   load(lds_taxa_file)
# }
# lds$taxons <- NA
# lds$taxons <- hier

# # extract out relevant rank levels
# lds <- lds %>% ungroup() %>% mutate(kingdom = map_chr(taxons, extract_rank), 
#      phylum = map_chr(taxons, extract_rank, rank=c('Phylum', 'Division')),
#      subphylum = map_chr(taxons, extract_rank, rank=c('Subphylum', 'Subdivision'), default=NA),
#      family = map_chr(taxons, extract_rank, rank='Family', default=NA),
#      order = map_chr(taxons, extract_rank, rank='Order', default=NA),
#      superphylum = map_chr(taxons, extract_rank, rank=c('Superphylum', 'Superdivision'), 
#                            default=NA),
#      class = map_chr(taxons, extract_rank, rank='Class', default=NA))

# # manually fix taxonomy
# lds$phylum[lds$species == "Aphanarthrum subglabrum/glabrum"] <- 'Arthropoda'
# lds$class[lds$species == "Aphanarthrum subglabrum/glabrum"] <- 'Insecta'
# lds$order[lds$species == "Aphanarthrum subglabrum/glabrum"] <- '	Coleoptera'

# lds$phylum[lds$species == "Cecidostiba fungosa"] <- 'Arthropoda'
# lds$class[lds$species == "Cecidostiba fungosa"] <- 'Insecta'
# lds$order[lds$species == "Cecidostiba fungosa"] <- 'Hymenoptera'

# lds$phylum[lds$species == "Mesobuthus cyprius/gibbosus"] <- 'Arthropoda'
# lds$class[lds$species == "Mesobuthus cyprius/gibbosus"] <- 'Arachnida'
# lds$order[lds$species == "Mesobuthus cyprius/gibbosus"] <- 'Scorpiones'

# lds$phylum[lds$species == "Nasonia vitripennis"] <- 'Arthropoda'
# lds$class[lds$species == "Nasonia vitripennis"] <- 'Insecta'
# lds$order[lds$species == "Nasonia vitripennis"] <- 'Hymenoptera'

# lds$phylum[lds$species == "Neurospora crassa"] <- 'Ascomycota'
# lds$class[lds$species == "Neurospora crassa"] <- 'Sordariomycetes'
# lds$order[lds$species == "Neurospora crassa"] <- 'Sordariales'

# lds$phylum[lds$species == "Paracoccidioides brasiliensis PS2/PS3/S1"] <- 'Ascomycota'
# lds$class[lds$species == "Paracoccidioides brasiliensis PS2/PS3/S1"] <- 'Eurotiomycetes'
# lds$order[lds$species == "Paracoccidioides brasiliensis PS2/PS3/S1"] <- 'Onygenales'

# lds$phylum[lds$species == "Paracoccidioides lutzii"] <- 'Ascomycota'
# lds$class[lds$species == "Paracoccidioides lutzii"] <- 'Eurotiomycetes'
# lds$order[lds$species == "Paracoccidioides lutzii"] <- 'Onygenales'

# lds$phylum[lds$species == "Phlebotomus ariasi"] <- 'Arthropoda'
# lds$class[lds$species == "Phlebotomus ariasi"] <- 'Insecta'
# lds$order[lds$species == "Phlebotomus ariasi"] <- 'Diptera'

# lds$phylum[lds$species == "Pieris rapae"] <- 'Arthropoda'
# lds$class[lds$species == "Pieris rapae"] <- 'Insecta'
# lds$order[lds$species == "Pieris rapae"] <- 'Lepidoptera'

# lds$phylum[lds$species == "Plasmodium falciparum"] <- 'Apicomplexa'
# lds$class[lds$species == "Plasmodium falciparum"] <- 'Aconoidasida'
# lds$order[lds$species == "Plasmodium falciparum"] <- 'Haemosporida'

# # save the full thing
# save(lds, file='leffler_et_al_2012_updated.Rdata')

# # ditch the list columns
write_tsv(lds, 'leffler_et_al_2012_updated.tsv')

# lds %>% ggplot(aes(range_cat, range)) + geom_bar(stat='identity')


# lds2 <- lds %>% filter(!is.na(range)) %>% arrange(range) %>% 
#            filter(range_cat != 'Unknown') %>%
#            mutate(range_cat=as.factor(tolower(range_cat)))
# cols <- RColorBrewer:::brewer.pal(4, 'Spectral')

# plot(seq_along(lds2$species), log10(lds2$range), col=cols[lds2$range_cat], pch=19)
# text(seq_along(lds2$species), log10(lds2$range), lds2$species, col=cols[lds2$range_cat])
# # legend()


# lds2 %>% 
#   group_by(species, range_cat) %>% summarize(range=mean(range, na.rm=TRUE)) %>%
#   ungroup() %>%
#   arrange(range) %>%
#   mutate(order = 1:nrow(.)) %>% 
#   ggplot(aes(order, log10(range), color=range_cat, label=species)) + 
#   geom_point() + geom_text_repel(size=1.8)

# lds2 %>% 
#   filter(kingdom == 'Animalia') %>%
#   group_by(species, range_cat) %>% 
#   summarize_at(c('range', 'size'), mean, na.rm=TRUE) %>%
#   ungroup() %>%
#   arrange(range/size) %>%
#   mutate(order = 1:nrow(.)) %>% 
#   ggplot(aes(order, log10(range/size), color=range_cat, label=species)) + 
#   geom_point() + geom_text_repel(size=1.8)

# lds2 %>% filter(genus == 'Pinus') %>% select(species, range, size)

# lds2 %>%
#   group_by(family) %>%
#   summarize_at(c('diversity', 'range', 'size'), mean, na.rm=TRUE) %>%
#   filter(range/size < 1e10) %>%
#   ggplot(aes(range/size, diversity, label=family)) +
#   geom_smooth(method='lm') +
#   geom_point() + geom_text_repel(size=1.8) + scale_x_log10() + 
#   scale_y_log10()

# d <- lds2 %>%
#   group_by(family) %>%
#   summarize_at(c('diversity', 'range', 'size'), mean, na.rm=TRUE)  %>%
#   filter(log10(range/size) < 1e10)

# summary(lm(log10(diversity) ~ log10(range / size), data=d))
 
