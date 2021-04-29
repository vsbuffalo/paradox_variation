# stapley_leffler_combined.r -- make main dataset, add in range data
## CACHE FILES
# combined_dataset_hier.Rdata: taxanomical hierarchy 
# redlist.Rdata: redlist species, downloaded from web
# combined_tree.Rdata: tree of life synthetic tree (no branch lengths)
# datalife.Rdata: synthetic calibrated tree
# popsize_chains.Rdata: Stan chains (these created by popsize_models.r, 
#                       which is sourecd here)

# NOTE: the Damuth data are not dry body mass -- I checked a few species.
# neither are Romiguier et al.

## OUTPUT DATA
# size_mass_density_lms.Rdata: linear models on Damuth 1987 data
# calibrated_phylogeny.Rdata: data for metazoan taxa in tree, and calibrated tree
# redlist_popsizes.tsv: population sizes from redlist data
# combined_data.tsv: the main dataset, for all metazona taxa (not parpath), even those not in phylogeny

library(ape)
library(phytools)
library(datelife)
library(rotl)
library(scales)
library(tidyverse)
library(rredlist)
library(ritis)
map <- purrr::map
source('../R/utilities.r')

opar <- par(no.readonly=TRUE)

#### Load all in Raw Datasets
dl <- read_tsv('leffler_et_al_2012_updated.tsv') %>% 
        # Leffler merged some species; we average them here.
        # we also merge all site types -- TODO, may want to do it both ways?
        group_by_at(vars(-size, -diversity, -site_type)) %>% 
        summarize_all(mean, na.rm=TRUE) %>%
        mutate(diversity_data_source = 'Leffler et al.')

ds <- read_tsv('stapley_et_al_2017_updated.tsv') %>% 
        select(-group, -log10_size, -data_source) %>%
        mutate(maplen_data_source = 'Stapley et al.') %>%
        mutate(diversity = NA)

dc <- read_tsv('corbett_detig_2015_updated.tsv')  %>%
        mutate(diversity_data_source = 'Corbett-Detig et al.',
               maplen_data_source =    'Corbett-Detig et al.',
               # rename the original Corbett-Detig range
               range_cd = range) %>%
        rename(diversity = obs_pi) %>%
        # 
        select(-range, -kingdom)


# Corbett-Detig's data has Dmel as 1.5mm which is 60% of what it should be.
# This was leading the average length across datasets to be too small, 
# and in turn leading to a much larger density, 10^9.4 rather than 10^9
# which is what it is for other species 
# I correct this here manually, based on Wikipedia's value (REVISIONS)
dc[dc$species == 'Drosophila melanogaster', ]$size <- 2.5e-3

# Corbett-Detig's data has przewalskii's horse's range as 10^7. These strange
# little horsies have a range that's more about ~300km² forcing it to be NA
# here should allow the correction in in species_range.r data to go through
dc[dc$species == 'Equus ferus przewalskii', ]$log10_range <- NA
dc[dc$species == 'Equus ferus przewalskii', ]$range_cd <- NA

# Romiguier et al data
dr <- read_tsv('romiguier_et_al_2014_updated.tsv') %>%
        mutate(diversity_data_source = 'Romiguier et al.') %>%
        select(species, genus, diversity = piS, fecundity, 
               propagule, size, diversity_data_source, bodymass_g)
#### Primary Data Joining
# This big ol' pipe monster is the main data cleaning/processing pipeline 
# as we merge in all the datasets. Many merges require that numeric 
# columns are averaged, which can leave NaNs when two NAs are averaged. 
# This requires lots of repetitive steps to mix, as all the datasets are
# combined.
data_source_merge <- function(x, y) {
  # first, anything not NA takes priority, then x
  out <- x
  out[is.na(out)] <- y
  out
}

d <- full_join(dl, ds, by=c('species', 'genus')) %>%
        # combine the data source columns
        rowwise() %>%
        mutate(diversity=mean(c(diversity.x, diversity.y), na.rm=TRUE),
               size=mean(c(size.x, size.y), na.rm=TRUE)) %>%
        # if two NAs are averaged, they lead to an NaN, so we fix this
        mutate(diversity=ifelse(is.nan(diversity), NA, diversity),
               size=ifelse(is.nan(size), NA, size)) %>%
        select(-size.x, -size.y, 
               -diversity.x, -diversity.y) %>%
        full_join(dc, by=c('species')) %>%
        # again, combine the data source columns
        mutate(diversity_data_source=data_source_merge(diversity_data_source.x, 
                                                       diversity_data_source.y)) %>%
        mutate(maplen_data_source=data_source_merge(maplen_data_source.x, 
                                                    maplen_data_source.y)) %>%
        select(-diversity_data_source.x, -diversity_data_source.y) %>%
        select(-maplen_data_source.y, -maplen_data_source.x)  %>%
        # merge in the Corbett-Detig diversities
        rowwise() %>%
        mutate(diversity=mean(c(diversity.x, diversity.y), na.rm=TRUE),
              size = mean(c(size.x, size.y), na.rm=TRUE),
              map_length=mean(c(map_length.x, map_length.y), na.rm=TRUE)) %>%
        # if two NAs are averaged, they lead to an NaN, so we fix this
        mutate(diversity=ifelse(is.nan(diversity), NA, diversity),
               size=ifelse(is.nan(size), NA, size),
               genome_size = mean(c(genome_size.x, genome_size.y), na.rm=TRUE),
               map_length=ifelse(is.nan(map_length), NA, map_length)) %>%
        select(-diversity.x, -diversity.y, 
               -size.x, -size.y,
               -genome_size.x, -genome_size.y,
               -map_length.x, -map_length.y) %>%
        full_join(dr, by=c('species', 'genus')) %>%
        mutate(diversity_data_source=data_source_merge(diversity_data_source.x, 
                                                       diversity_data_source.y)) %>%
        select(-diversity_data_source.x, -diversity_data_source.y) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(diversity = mean(c(diversity.x, diversity.y), na.rm=TRUE),
               size = mean(c(size.x, size.y), na.rm=TRUE)) %>%
        # again, if two NAs are averaged, they lead to an NaN, so we fix this
        mutate(diversity=ifelse(is.nan(diversity), NA, diversity),
               size=ifelse(is.nan(size), NA, size)) %>%
        select(-size.x, -size.y, -diversity.x, -diversity.y) %>%
        ungroup() %>%
        # now, deal with duplicates. We average numeric columns
        group_by_if(is.character) %>% summarize_all(mean, na.rm=TRUE) %>%
        mutate(diversity=ifelse(is.nan(diversity), NA, diversity),
               size=ifelse(is.nan(size), NA, size),
               map_length=ifelse(is.nan(map_length), NA, map_length)) %>%
        # fix the parpath column (only in Stapley data)
        mutate(parpath = as.logical(ifelse(is.nan(parpath), NA, parpath))) %>%
        # turn all numeric column's nans to NAs
        mutate_if(is.numeric, ~ ifelse(is.nan(.), NA, .)) %>%
        mutate(log10_diversity = log10(diversity))

#### Manual Curation of social/parpath species
# deal with remaining parpath entried -- looked these up on Wikipedia / Google
# External parasites (Ixodes), and parasitoid wasps are not counted.
parpath <- c("Phytophthora capsici", # an oomycete
             "Plasmodium vivax", #  maralia
             "Meloidogyne hapla", # vegetable nematode, a plant pathogen
             "Schistosoma mansoni") # blood fluke, a human parasite

d$parpath <- case_when(is.na(d$parpath) & d$species %in% parpath ~ TRUE,
                       is.na(d$parpath) & !(d$species %in% parpath) ~ FALSE,
                       TRUE ~ d$parpath)

# list of social insects, for recomb plot
social <- c("Acromyrmex echinatior",
            "Apis mellifera",
            "Apis cerana",
            "Bombus terrestris",
            # wasp
            "Vespula vulgaris",
            # ants
            "Pogonomyrmex rugosus",
            # termites 
            "Reticulitermes flavipes", 
            "Reticulitermes grassei", 
            "Reticulitermes lucifugus")

d$social <- FALSE
d$social[d$species %in% social] <- TRUE

stopifnot(nrow(d) == length(unique(d$species)))

#### Load in Range Data Estimated from Occurrence Data
## Load and clean up range size
ranges <- read_tsv('gbif_ranges.tsv')  %>%
              # post-hoc fixes
              # some Ciona entries -- GBIF data was downloaded with old name
              mutate(key=ifelse(key == 'Ciona intestinalis A', 
                                  'Ciona intestinalis', key)) %>%
              mutate(key=ifelse(key == 'Mus musculus castaneus', 
                                  'Mus musculus', key)) %>%
              select(-genus)

merge_ranges <- function(x) {
  if (nrow(x) == 0) return(tibble(is_terrestrial = NA, n_occ=NA, range=NA))
  out <- x %>% summarize(is_terrestrial = unique(is_terrestrial),
                         n_occ=mean(n_occ, na.rm=TRUE),
                         range=mean(range, na.rm=TRUE))
  stopifnot(nrow(out) == 1)
  return(out)
}

d <- d %>% 
         # before merge, we want to perserve the "good" species names 
         # not combined species keys as in original Leffler et al
         mutate(species_orig = species) %>%
         nest_join(ranges, by=c('species'='key')) %>%
         mutate(ranges = map(ranges, merge_ranges)) %>%
         unnest(ranges) %>%
         ungroup() %>%
         rowwise() %>%
         # this averages the ranges from Corbett-Detig and mine
         # ignoring NAs
         mutate(range = mean(c(range, range_cd), na.rm=TRUE),
                range = ifelse(is.nan(range), NA, range)) %>%
         select(-range_cd) 
stopifnot(nrow(d) == length(unique(d$species)))
stop()

#### Get Taxonomical Info
# get the taxonomial datas from taxadb (via Catalog of Life)
d <- d %>% mutate(id_raw=map_chr(clean_name(species), taxadb::get_ids)) 
# find all missing IDs and get genera-level IDs
missing_sps <- d$species[is.na(d$id_raw)]
genera <- map_chr(missing_sps, extract_genus, genus_only=TRUE)
genera_ids <- map(genera, taxadb::get_ids)
d$id_raw[is.na(d$id_raw)] <- genera_ids
d <- d %>% mutate(id = map_int(id_raw, ~ as.integer(sub('ITIS:', '', .))))
 
## Lookup and merge in taxonomy 
TAXA_FILE <- 'combined_dataset_hier.Rdata'
if (!file.exists(TAXA_FILE)) {
 # grab taxonomical info
  hier <- map(d$id, possibly(hierarchy_full, NA))
  save(hier, file=TAXA_FILE)
} else {
  load(TAXA_FILE)
}
d$taxons <- NA
d$taxons <- hier


## Extract taxonomical groups
d <- d %>% ungroup() %>% 
    # a little hack to move genus in the col order
    mutate(.genus = genus) %>% 
    select(-genus) %>% 
    mutate(kingdom = map_chr(taxons, extract_rank), 
           superphylum = map_chr(taxons, extract_rank, rank=c('Superphylum', 'Superdivision'), 
                                 default=NA),
           phylum = map_chr(taxons, extract_rank, rank=c('Phylum', 'Division')),
           subphylum = map_chr(taxons, extract_rank, rank=c('Subphylum', 'Subdivision'), 
                               default=NA),
           class = map_chr(taxons, extract_rank, rank='Class', default=NA),
           order = map_chr(taxons, extract_rank, rank='Order', default=NA),
           family = map_chr(taxons, extract_rank, rank='Family', default=NA),
           # extract genus if its NA
           genus = ifelse(is.na(.genus), sub('(\\w+) .*', '\\1', species), .genus),
           clean_name = map_chr(species, clean_name)) %>%
    select(-.genus)
           

## Manual taxonomy fixes
source('taxonomy_fixes.r')
for (fix in taxon_fixes) {
  sps <- fix$species
  for (level in taxon_levels) {
    stopifnot(any(d$species == sps)) # makes sure we find a species to change
    d[d$species == sps, level] <- fix[[level]]
  }
}


#### Primary Filtering
# This study uses only metazoan taxa, and we ignore the parpath species 
# which since they (1) have dymamics w/ hosts and (2) have low occurrence 
# data, have shoddy estimates of range sizes.
# Note that we have range data for plant taxa too, which we do not use here
# (and it too is generally less reliable). Still, it is useful for debugging
# the range estimation process.
d <- d %>% 
      filter(kingdom == 'Animalia') %>%
      filter(!parpath)


#### IUCN Red List
# IUCN Red List entries
# this is used for quality control in range estimation
REDLIST_FILE <- 'redlist.Rdata'
if (!file.exists(REDLIST_FILE)) {
  message('no redlist file found, downloading... ')
  # animals only
  da <- d %>% filter(kingdom == 'Animalia')
  safe_rl_search <- possibly(rl_search, NA)
  rdl_entries <- map(da$species, safe_rl_search)
  save(rdl_entries, file=REDLIST_FILE)
} else {
  load(REDLIST_FILE)
}

# get the occurrences / km2 out
species <- map(rdl_entries, 'name') 
species <- map_chr(species, ~ ifelse(is.null(.), NA, .))

get_col <- function(x, col='eoo_km2', raw=FALSE) {
  if (col %in% colnames(x)) {
    text <- x[, col]
    if (raw) return(text)
    if (is.na(text)) return(NA)
    if (regexpr('-', text, fixed=TRUE) != -1) {
      # deal with a range
      return(mean(as.numeric(unlist(strsplit(text, '-', fixed=TRUE)))))
    } 
    return(as.numeric(sub('>', '', text)))
  }
  return(NA)
}

redlist_cat <- map_chr(map(rdl_entries, 'result'), get_col, col='category',
                       raw=TRUE)
# the lc "indicates they have not been re-evaluated since 2000" - WP
# I just convert -- affects very few
redlist_cat[redlist_cat == 'LR/lc'] <- "LC"
redlist_cat[redlist_cat == 'LR/nt'] <- "NT" # upgrade to NT
redlist_cat[redlist_cat == 'LR/nt'] <- "LC"
redlist_df <- tibble(species, redlist_cat)

# merge in IUCN cats 
d <- d %>% left_join(redlist_df, by='species')  

eoo_km2 <- map_dbl(map(rdl_entries, 'result'), get_col) 
aoo_km2 <- map_dbl(map(rdl_entries, 'result'), get_col, col='aoo_km2') 
rdl_popsizes <- tibble(species = species, 
                       eoo_km2=eoo_km2) %>%
                 filter(!is.na(eoo_km2))

# build a small dataset of cases where we have real EOO / km2 
# data, to connect our RS proxy to real counts
rdl_popsizes  <- d %>% mutate(RS = range / size) %>% 
  select(species, kingdom, family, RS, range, size) %>%
  right_join(rdl_popsizes, by='species') %>%
  mutate(D_proxy = log10(1/size))

write_tsv(rdl_popsizes, 'redlist_popsizes.tsv')

# EOO km2 is a range estimate -- we can validate our range 
# measures this way.
ggplot(rdl_popsizes, aes(log10(eoo_km2), log10(range), color=family)) + geom_point()

# pred_eoo_km2 <- function(RS) {
#   if (is.na(RS)) return(NA)
#   predict(rdl_fit, newdata=tibble(RS = RS))
# }


#### Estimate the body mass / body size relationship from Romiguier et al data
# these are simple linear models predictions -- now this is handled all in Stan 
# in the full dataset; kept here for comparison
dr <- dr %>% mutate(log10_mass_g = log10(bodymass_g), 
                    log10_size = log10(size))
size_mass_fit <- lm(log10_mass_g ~ log10_size, data=dr)
summary(size_mass_fit)
plot(log10(bodymass_g) ~ log10(size), data=dr)
abline(size_mass_fit)

pred_log10_mass <- function(log10_size) {
  # predicts in grams
  stopifnot(length(log10_size) == 1)
  if (is.na(log10_size)) return(NA)
  predict(size_mass_fit, newdata=tibble(log10_size = log10_size))
}

#### Load in the body mass / density data from Damuth (1987)
# Again this is for some simplistic density estimates from 
# linear models; the final data uses Stan

# Damuth 1987 data
dd <- read_tsv('./damuth1987.tsv') 

# poikilotherms densities are divided by 30
poikilotherms <- c('pisces', 'reptilian', 'amphibia', "terrestrial arthropods", 
                   "other terrestrial invertebrates")
terrestrial_animals <- c("mammals", "amphibia", "reptilian", 
                         "terrestrial arthropods", "other terrestrial invertebrates")
dd <- dd %>% mutate(poikilothermic = group %in% poikilotherms) %>%
       mutate(density = density) %>%
       mutate(raw_density = density, 
              density = ifelse(poikilothermic, density / 30, density)) %>%
        mutate(log10_mass_g = log10(mass_g), log10_density = log10(density)) %>%
        filter(group %in% terrestrial_animals)

mass_density_fit <- lm(log10_density ~ log10_mass_g, data=dd)
summary(mass_density_fit)
predict(mass_density_fit, newdata=data.frame(log10_mass_g=-3.6))

dd_groups <- dd %>% group_by(group) %>% nest() %>% 
               mutate(fit = map(data, ~ lm(log10_density ~ log10_mass_g, data=.))) 

ggplot(dd) + geom_point(aes(log10_mass_g, log10_density, color=group)) + 
  geom_smooth(mapping=aes(log10_mass_g, log10_density), method='lm') +
  scale_x_continuous(breaks = seq(-5, 5, 2)) + 
  scale_y_continuous(breaks = seq(-3, 9, 2), limits=c(-3, 9)) 

# save the linear models
save(size_mass_fit, mass_density_fit, file='size_mass_density_lms.Rdata')

# for predicting density based on lm on Damuth's data
pred_log10_density <- function(log10_mass_g) {
  stopifnot(length(log10_mass_g) == 1)
  if (is.na(log10_mass_g)) return(NA)
  predict(mass_density_fit, newdata=tibble(log10_mass_g=log10_mass_g))
}


#### Stan Population Size Estimates
# Previous popsize measures use linaer models and don't carry forward uncertainty.
# I have also coded up a Stan model that estimates the size / body mass
# relationship from the Romiguier et al data, and the body mass / density
# relationship from the Damuth 1987 dataset, and then uses posterior draws to
# estimate the population sizes as range * density.

source('./popsize_models.r')

# how to the lm and Stan estimates compare? 
sapply(rstan:::extract(bayes_fit, pars), mean)

# augment data, estimate the popsize proxies
do <- d # for debugging
# d <- do


# Some of the naming is bad here. To clarify:
# 1. pred_log10_body_mass and pred_log10_density are lm predictions
# 2. log10_density_pred and log10_body_mass_pred are random samples from the
# the Stan model using posterior parameter samples

## Join in the Stan data
d <- d %>%
       left_join(bayes_popsizes) %>%
       mutate(log10_size = log10(size)) %>%
       # predicted body size from the lm -- in grams
       mutate(pred_log10_body_mass=map_dbl(log10_size, pred_log10_mass))  %>%
       # predicted density from the lm 
       mutate(pred_log10_density=map_dbl(pred_log10_body_mass, 
                                         pred_log10_density)) %>%
       # Nc proxy
       mutate(pred_log10_N = log10(10^pred_log10_density * range)) %>%
       mutate(RS = range / size, 
              log10_RS=log10(RS)) %>%
       mutate(log10_size = log10(size),
              log10_range = log10(range),
              log10_map_length=log10(map_length))

# these are diagnostics to ensure a close fit between
# the Stan lognormal model and the log-log regression models
# basically, a comparison to lm and non-regularized versions
da <- d %>% filter(kingdom == 'Animalia')
plot(log10_density_pred ~ pred_log10_density, da); abline(a=0, b=1)
plot(log10_body_mass_pred ~ pred_log10_body_mass, da); abline(a=0, b=1)
plot(log10_body_mass ~ pred_log10_body_mass, da); abline(a=0, b=1)
plot(log10_popsize ~ pred_log10_N, da); abline(a=0, b=1)

# popsize range posterior
hist(psr)

#### Tree of Life Phylogenetic 
# No branch lengths; those are added by datelife later.
# Lots of taxa matching first.
# how many still have missing taxonomical data? Two species  --- fine.
d %>% filter_at(vars(kingdom, phylum, class, order, family, genus), any_vars(is.na(.))) %>%
  select(species, kingdom, phylum, class, order, family, genus)


# get taxa
all_taxa_df <- tnrs_match_names(unique(d$clean_name)) %>% as_tibble() 

## Phylognetic Tree
# only for animal species -- since the plant size data is to arbitrary 

taxa <- tnrs_match_names(unique(sort(da$clean_name)), context_name = "Animals") 
taxa_df <- taxa %>% as_tibble()

# what significantly differs?
taxa_df %>% filter(tolower(unique_name) != search_string) %>%
  select(search_string, unique_name) %>% as.data.frame()

# inspect the number of multiple matches
table(taxa$number_matches)

# what has more than one match? do the approximate matches look reasonable?
taxa_df %>% filter(number_matches > 1) %>% 
  select(search_string, unique_name)

# these all seem reasonable -- bogotana is close enough to pseudoobscura

# get the tree of life tree
TREE_FILE <- 'combined_tree.Rdata'
if (!file.exists(TREE_FILE)) {
  tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])
  save(tr, file=TREE_FILE)
} else {
  load(file=TREE_FILE)
}
tro <- tr # make a copy


# now, match the original dataframe with the tree data
da <- da %>% mutate(species_orig = species) %>%
             mutate(species_lower = tolower(clean_name))

# removed some extraneous stuff
tips_simple <- sub('(species_in_domain_Eukaryota)_', '', tro$tip.label, fixed=TRUE)

# A lot of this is complicated, as we need to join the tree tips
# without taxa names, then with our original dataset.
# The taxa lookups mess the names up a bit, which is why this is so hard.
da_tr <- tibble(ott_name=tro$tip.label, tip_species = tips_simple) %>%
           tidyr:::extract(tip_species, into=c('genus', '.species', 'ott_id'), 
                   '([^_]+)_([\\w_]+)_ott(\\d+)',
                   remove=FALSE, convert=TRUE) %>%
           # deal with subspecies names
           separate(.species, into=c('.species', 'subspecies')) %>%
           unite(species, genus, .species, sep=' ', remove=FALSE)  %>%
           select(-.species) %>%
           # join in the taxa IDs by OTT
           left_join(taxa_df, by='ott_id') %>%
           # join in the main dataframe, by the taxa search string 
           # (lower case version of the species) 
           left_join(da, by=c(search_string='species_lower'))  %>%
           # we keep species, genus from the taxa_df (x) -- this should be more
           # accurate
           select(-genus.y, -species.y) %>% 
           rename(species = species.x, genus = genus.x)
 

# match trees with tips (a precaution; should be in same order)
idx <- match(tro$tip.label, da_tr$ott_name)
tr <- drop.tip(tro, tr$tip.label[is.na(idx)])
tr$tip.label <- da_tr$species[idx]
# order dataframe
da_tr <- da_tr[idx, ]

#### Build Datelife Calibrated Phylogeny
# Use datelife to calibrate the phylogeny
# note: the API for this is rapidly change; it's still under 
# rapid development as I'm using it. This is why I've still uesd tol_induced_subtree()
# and manually done some steps with other packages -- I've encountered some bugs
# due to beta nature of the software
DATELIFE_FILE <- 'datelife.Rdata'
if (!file.exists(DATELIFE_FILE)) {
  met_dq <- make_datelife_query(da_tr$clean_name)
  met_dr <- get_datelife_result(input=met_dq)
  met_phyloall <- summarize_datelife_result(datelife_result = met_dr)
  # plot_phylo_all(met_phyloall)
  met_syntree <- get_otol_synthetic_tree(input = met_dq)
  met_allcal <- get_all_calibrations(met_phyloall)
  met_cal <- use_calibrations(phy=met_syntree, met_allcal)
  save(met_dq, met_dr, met_syntree, met_allcal, met_cal, file=DATELIFE_FILE)
} else {
  load(DATELIFE_FILE)
}

tips_clean <- function(x) {
  # again, remove junk added to tip labels
  out <- trimws(gsub('_', ' ', sub('-species_in_domain_Eukaryota', '', x, fixed=TRUE)), 'both')
  gsub('(\\w+) (\\w+).*', '\\1 \\2', out)
}

clean_tip_labels <- tips_clean(met_cal$tip.label)
met_cal$tip.label <- clean_tip_labels

# ladderize phylogeny
met_cal_rooted <- ladderize(root(met_cal, 'Amphimedon queenslandica'))

# match up the tree data with the tip labels
idx <- match(met_cal$tip.label, da_tr$species)
# did we lose anything? ensure not.
stopifnot(length(met_cal$tip.label[is.na(idx)]) == 0)
# notice: we only keep data from *calibrated* phylogeny
da_tr <- da_tr[idx, ]
stopifnot(nrow(da_tr) == length(met_cal$tip.label))
# again, assert we have the same numbers of tips and entries in our subset of
# the Stapley data
stopifnot(nrow(da_tr) == length(tr$tip.label))
# rename things
otl_tr <- tr
tr <- met_cal



#### Save the data
# Save the combined data set
d_flat <- da %>% select(-id_raw, -taxons)
write_tsv(d_flat, 'combined_data.tsv')

# save the phylogenetic data
save(da_tr, tr, otl_tr, file='calibrated_phylogeny.Rdata')





