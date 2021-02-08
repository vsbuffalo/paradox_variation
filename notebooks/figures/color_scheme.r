library(tidyverse)
library(scales)

# we want consistent phylum colors, which we set here

xx <- read_tsv('../../data/combined_data.tsv')  %>%
       filter(kingdom == 'Animalia') 


nphyla <- length(unique(xx$phylum))
phyla_cols <- hue_pal()(nphyla)
names(phyla_cols) <- sort(unique(xx$phylum))


 
