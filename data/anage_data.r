library(tidyverse)

ctypes <- cols(
  HAGRID = col_character(),
  Kingdom = col_character(),
  Phylum = col_character(),
  Class = col_character(),
  Order = col_character(),
  Family = col_character(),
  Genus = col_character(),
  Species = col_character(),
  `Common name` = col_character(),
  `Female maturity (days)` = col_double(),
  `Male maturity (days)` = col_double(),
  `Gestation/Incubation (days)` = col_double(),
  `Weaning (days)` = col_double(),
  `Litter/Clutch size` = col_double(),
  `Litters/Clutches per year` = col_double(),
  `Inter-litter/Interbirth interval` = col_double(),
  `Birth weight (g)` = col_double(),
  `Weaning weight (g)` = col_double(),
  `Adult weight (g)` = col_double(),
  `Growth rate (1/days)` = col_double(),
  `Maximum longevity (yrs)` = col_double(),
  Source = col_character(),
  `Specimen origin` = col_character(),
  `Sample size` = col_character(),
  `Data quality` = col_character(),
  `IMR (per yr)` = col_double(),
  `MRDT (yrs)` = col_double(),
  `Metabolic rate (W)` = col_double(),
  `Body mass (g)` = col_double(),
  `Temperature (K)` = col_double(),
  References = col_number()
)


cnames <- c(
  'HAGRID',
  'kingdom',
  'phylum',
  'class',
  'order',
  'family',
  'genus',
  'species',
  'common_name',
  'female_maturity_days',
  'male_maturity_days',
  'gestation_incubation_days',
  'weaning_days',
  'litter_clutch_size',
  'litters_clutches_per_year',
  'inter_litter_interbirth_interval',
  'birth_weight_g',
  'weaning_weight_g',
  'adult_weight_g',
  'growth_rate_recip_days',
  'maximum_longevity_yrs',
  'source',
  'specimen_origin',
  'sample_size',
  'data_quality',
  'imr_per_yr',
  'mrdt_yrs',
  'metabolic_rate_w',
  'body_mass_g',
  'temperature_k',
  'references'
)

da_raw <- read_tsv('./anage_data.tsv', col_types=ctypes)
colnames(da_raw) <- cnames

da <- da_raw %>% 
         mutate(species = map2_chr(genus, species, ~ paste(.x, .y, sep=" ")))
stop()
# write_tsv(da, '../data/anage_data_updated.tsv')

taxon_keys <- c('species', 'phylum', 'kingdom', 'order', 'class', 'genus', 'family')

ggplot(da, aes(body_mass_g, litter_clutch_size, color=phylum))  + geom_point() + scale_x_log10() + scale_y_log10()

ggplot(da, aes(body_mass_g, litters_clutches_per_year))  + geom_point() + scale_x_log10() + scale_y_log10() + geom_smooth(method='lm')
ggplot(da, aes(body_mass_g, female_maturity_days))  + geom_point() + scale_x_log10() + scale_y_log10() + geom_smooth(method='lm')

ggplot(da, aes(body_mass_g, litters_clutches_per_year))  + geom_point() + scale_x_log10() + scale_y_log10()



# dsc %>% inner_join(da, by=taxon_keys) %>% 
#   ggplot(aes(map_length, log10(litter_clutch_size * litters_clutches_per_year), color=family)) + geom_point()

# drs %>% inner_join(da, by='species') %>%
#   mutate(psize = (weaning_weight_g/adult_weight_g)^3) %>%
#   ggplot(aes(psize, propagule_size_cm)) + 
#   geom_point()


