# data.r -- join all the Corbet-Detig datasets
library(tidyverse)

maps_raw <- readxl:::read_xls('journal.pbio.1002112.s016.XLS')
maps <- maps_raw[, c(1:3, 6, 10)]
colnames(maps) <- c('species', 'num_chr', 'map_length', 'totalMb_covered', 'ave_rec')

div_cord <- readxl:::read_xls('S2_Data.xls')

colnames(div_cord) <- c("species", "abbrv", "kingdom", "size", 
                        "range", "pred_pi",
                        "obs_pi", "impact_of_sel", 
                        "log10_size", "log10_range", 
                        "neut_model_conf")

# genome size
gs <- readxl:::read_xls('journal.pbio.1002112.s018.XLS')[c(1, 2, 4, 10)]
colnames(gs) <- c('species', 'genome_species', 'assembly_size', 'genome_size')

# selfing
self <- readxl:::read_xls('journal.pbio.1002112.s017.XLS')[c(1, 3, 4)]
colnames(self) <- c('species', 'selfing', 'self_est')


# Below is the Table S5 of McGaugh et al. 2012
# for reference.
# This paper is cited as the map source for D. pseudoobscura 
# in Corbett-Detig et al (2015). The map length of 
# pseudoobscura in CD2015 is 0.535 (CD2015 Table S13)
# which seems very low. McGaugh et al say they only use 43%
# of the chromosome.
mcg12 <- tribble(~ Mb, ~ cMMb,
                 0.021916,10.38,
                 0.019957,4.67,
                 0.017055,3.64,
                 0.021615,2.87,
                 0.024667,6.29,
                 0.020844,21.25,
                 0.019964,6.9,
                 0.019212,3.59,
                 0.021546,4.11,
                 0.022377,3.52,
                 0.021324,13.85,
                 0.02736,4.31,
                 0.02314,8.07,
                 0.012614,2.34,
                 0.018506,1.59,
                 0.017728,6.1)

# In looking the D. pseudoobscura data, I cannot
# recreate CD's estimate and it's a very strong 
# outlier, so I am ignoring this
maps$map_length[maps$species ==  "Drosophila pseudoobscura"] <- NA

# main join
dc <- dplyr:::left_join(div_cord, maps) %>% 
  dplyr::left_join(gs, by='species') %>%
  dplyr::left_join(self, by='species')


# ## Manual map length additions:

# # For some missing values, I have found their total map lengths in the
# # literature; these are entered below

# # from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2728870/, table 1
# dc$map_length[dc$species == 'Mus musculus castaneus'] <- 1445.23	

# # assuming same genetic map
# # from https://www.frontiersin.org/articles/10.3389/fpls.2017.00706/full
# dc$map_length[dc$species == 'Zea mays ssp parviglumis'] <- 2236.66


# if we assume Equus ferus map length ~ E. 
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC310772/
dc$map_length[dc$species == 'Equus ferus przewalskii'] <- 2000

# from closely related B. mori https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-1-r21
dc$map_length[dc$species == 'Bombyx mandarina'] <- 1413

# https://genomebiology.biomedcentral.com/articles/10.1186/gb-2007-8-4-r66
dc$map_length[dc$species == 'Apis mellifera scutellata'] <- 4114.5

# # two related species; as a rough approx, we average them
# # https://www.cambridge.org/core/journals/fruits/article/integrated-genetic-map-of-citrus-based-on-rapd-markers/78D957230473160950D924DCD2A7A68A
# # https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-593
# dc$map_length[dc$species == 'Citrus reticulata'] <- (527 + 1084.1)/2

# # https://www.ncbi.nlm.nih.gov/pubmed/30347873
# dc$map_length[dc$species == 'Citrullus lanatus lanatus'] <- 1508.94

# # https://www.cell.com/molecular-plant/pdf/S1674-2052(15)00176-8.pdf
# dc$map_length[dc$species == 'Cucumis sativus'] <- 1384.4  

# # https://www.sciencedirect.com/science/article/pii/S0888754319303945?via%3Dihub
# # the cross between the two G. soja varities (pop II); could consider averaging in HI too
# dc$map_length[dc$species == 'Glycine soja'] <- 2622.8 

# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1456771/
# dc$map_length[dc$species == 'Heliconius melpomene melpomene'] <- 1616

# # https://www.ncbi.nlm.nih.gov/pubmed/20920197
# dc$map_length[dc$species == 'Ovis canadensis'] <- 3051

# # a relative, Prunus mume
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4463843/
# dc$map_length[dc$species == 'Prunus davidiana'] <- 1550.62

# # https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3430-7
# dc$map_length[dc$species == 'Sorghum bicolor subsp. verticilliflorum'] <- 2158.1

# # https://www.frontiersin.org/articles/10.3389/fpls.2016.00437/full
# # worth seeing Table 5; this map is a bit larger than other published records but looks
# # high quality
# dc$map_length[dc$species == "Cucumis sativus var. hardwickii"] <- 1046.97


# Load the Wilfert et al (2007) data
# I have used various sources to approximate the body size in meters. 
# For entries left blank, I either could not find a reliable size, or 
# the data already existed in Corbett-Detig et al.

## Commented out for now, since the Stapley et al data is much larger
# dwf <- read_csv('wilfert_et_al_2007.csv')
# dwfs <- dwf[!is.na(dwf$size), c('species', 'class', 'genetic_size_cM', 'size')]
# colnames(dwfs) <- c('species', 'class', 'map_length', 'size')

# dd <- full_join(d, dwfs)
# dd$log10_size <- log10(dd$size)


# we clean up some names; simplifying since we need to do automated taxonomical lookups
related <- c("Apis mellifera scutellata"="Apis mellifera",
	"Citrullus lanatus lanatus"="Citrullus lanatus",
	"Cucumis sativus var. hardwickii"="Cucumis sativus", 
	"Heliconius melpomene melpomene"="Heliconius melpomene",
	"Sorghum bicolor subsp. verticilliflorum"="Sorghum bicolor",
	"Zea mays ssp parviglumis"="Zea mays")

tmp <- related[dc$species]
tmp[is.na(tmp)] <- dc$species[is.na(tmp)]
dc$species <- tmp

# clean up selfing data
dc$selfing <- c('no'=FALSE, 'yes'=TRUE)[dc$selfing]

dc <- dc %>% 
        # rename kingdoms
        mutate(kingdom = 
               case_when(kingdom == 'animal' ~ 'Animalia',
                         kingdom == 'plant' ~ 'Plantae',
                         TRUE ~ NA_character_)) %>%
        # convert map lengths to Morgans
        mutate(map_length = map_length / 100)
              
write_tsv(dc, path='corbett_detig_2015_updated.tsv')


