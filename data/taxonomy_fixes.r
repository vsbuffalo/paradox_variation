## Manual taxonomy fixes
taxon_levels <- c('kingdom', 'phylum', 'class', 'order', 'family')

taxon_fixes <- list(
     list(species = "Agaricus bisporus",
          kingdom = 'Fungi', 
          phylum='Basidiomycota', 
          class='Agaricomycetes', 
          order='Agaricales',
          family='Agaricaceae'),

     list(species = "Chimonanthus praecox",
          kingdom = 'Plantae', 
          phylum=NA, 
          class=NA, 
          order='Laurales',
          family='Calycanthaceae'),

     list(species = "Aphanarthrum subglabrum/glabrum",
          kingdom = 'Animalia', 
          phylum='Arthropoda', 
          class='Insecta', 
          order='Coleoptera',
          family='Curculionidae'),

     list(species="Cecidostiba fungosa",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Hymenoptera',
          family='Pteromalidae'),

     list(species="Mesobuthus cyprius/gibbosus",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Arachnida',
          order='Scorpiones',
          family='Buthidae'),

     list(species="Nasonia vitripennis",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Hymenoptera',
          family='Pteromalidae'),

     list(species="Neurospora crassa",
          kingdom='Animalia',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Sordariales',
          family='Sordariaceae'),

     list(species="Paracoccidioides brasiliensis PS2/PS3/S1",
          kingdom='Animalia',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Onygenales',
          family='Ajellomycetaceae'),

     list(species="Paracoccidioides lutzii",
          kingdom='Animalia',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Onygenales',
          family='Ajellomycetaceae'),

     list(species="Phlebotomus ariasi",
          kingdom='Fungi',
          phylum='Arthropoda',
          class='Insecta',
          order='Diptera',
          family='Psychodidae'),

     list(species="Pieris rapae",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Lepidoptera',
          family='Pieridae'),

     list(species="Plasmodium falciparum",
          kingdom='Chromista',
          phylum='Apicomplexa',
          class='Aconoidasida',
          order='Haemosporida',
          family='Plasmodiidae'),

     list(species="Abatus agassizi",
          kingdom='Animalia',
          phylum='Echinodermata',
          class='Echinoidea',
          order='Spatangoida',
          family='Schizasteridae'),

     list(species="Abatus cordatus",
          kingdom='Fungi',
          phylum='Basidiomycota',
          class='Agaricomycetes',
          order='Agaricales',
          family='Schizasteridae'),

     list(species="Aspergillus nidulans",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Eurotiales',
          family='Trichocomaceae'),

     list(species="Aspergillus nidulans",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Eurotiales',
          family='Trichocomaceae'),

     list(species="Bicyclus anynana",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Lepidoptera',
          family='Nymphalidae'),

     list(species="Bostrycapulus aculeatus",
          kingdom='Animalia',
          phylum='Mollusca',
          class='Gastropoda',
          order='Littorinimorpha',
          family='Calyptraeidae'),

     list(species="Cecidostiba fungosa",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Hymenoptera',
          family='Cynipidae'),

     list(species="Chimonanthus praecox",
          kingdom='Plantae',
          phylum=NA,
          class=NA,
          order=NA,
          family=NA),

     list(species="Coprinus cinereus",
          kingdom='Fungi',
          phylum='Basidiomycota',
          class='Agaricomycetes',
          order='Agaricales',
          family='Psathyrellaceae'),

     list(species="Eimeria maxima",
          kingdom='Chromista',
          phylum='Apicomplexa',
          class='Conoidasida',
          order='Eucoccidiorida',
          family='Eimeriidae'),

     list(species="Eimeria tenella",
          kingdom='Chromista',
          phylum='Apicomplexa',
          class='Conoidasida',
          order='Eucoccidiorida',
          family='Eimeriidae'),

     list(species="Eunicella cavolinii",
          kingdom='Animalia',
          phylum='Cnidaria',
          class='Anthozoa',
          order='Alcyonacea',
          family='Gorgoniidae'),

     list(species="Eunicella verrucosa",
          kingdom='Animalia',
          phylum='Cnidaria',
          class='Anthozoa',
          order='Alcyonacea',
          family='Gorgoniidae'),

     list(species="Fusarium circinatum",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Hypocreales',
          family='Nectriaceae'),

     list(species="Fusarium graminearum",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Hypocreales',
          family='Nectriaceae'),

     list(species="Fusarium oxysporum",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Hypocreales',
          family='Nectriaceae'),

     list(species="Fusarium verticillioides",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Hypocreales',
          family='Nectriaceae'),

     list(species="Globodera rostochiensis",
          kingdom='Animalia',
          phylum='Nematoda',
          class='Secernentea',
          order='Tylenchida',
          family='Heteroderidae'),

     list(species="Histoplasma capsulatum",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Onygenales',
          family='Ajellomycetaceae'),

     list(species="Kryptolebias marmoratus",
          kingdom='Animalia',
          phylum='Chordata',
          class='Actinopterygii',
          order='Cyprinodontiformes',
          family='Rivulidae'),

     list(species="Magnaporthe grisea",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Magnaporthales',
          family='Magnaporthaceae'),

     list(species="Melampsora lini",
          kingdom='Fungi',
          phylum='Basidiomycota',
          class='Pucciniomycetes',
          order='Pucciniales',
          family='Melampsoraceae'),

     list(species="Mellicta athalia",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Lepidoptera',
          family='Nymphalidae'),

     list(species="Mellicta parthenoides",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Lepidoptera',
          family='Nymphalidae'),

     list(species="Mesobuthus cyprius/gibbosus",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Arachnida',
          order='Scorpiones',
          family='Buthidae'),

     list(species="Nasonia giraulti",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Hymenoptera',
          family='Pteromalidae'),

     list(species="Nasonia giraulti",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Hymenoptera',
          family='Pteromalidae'),

     list(species="Nasonia vitripennis",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Hymenoptera',
          family='Pteromalidae'),

     list(species="Necora puber",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Malacostraca',
          order='Decapoda',
          family='Portunidae'),

     list(species="Neurospora crassa",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Sordariomycetes',
          order='Sordariales',
          family='Sordariaceae'),

     list(species="Paracoccidioides brasiliensis PS2/PS3/S1",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Onygenales',
          family='Ajellomycetaceae'),

     list(species="Paracoccidioides lutzii",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Eurotiomycetes',
          order='Onygenales',
          family='Ajellomycetaceae'),

     list(species="Phlebotomus ariasi",
          kingdom='Animalia',
          phylum='Arthropoda',
          class='Insecta',
          order='Diptera',
          family='Psychodidae'),

     list(species="Plasmodium chabaudi",
          kingdom='Chromista',
          phylum='Apicomplexa',
          class='Aconoidasida',
          order='Haemosporida',
          family='Plasmodiidae'),

     list(species="Plasmodium vivax",
          kingdom='Chromista',
          phylum='Apicomplexa',
          class='Aconoidasida',
          order='Haemosporida',
          family='Plasmodiidae'),

     list(species="Pseudocercospora fijiensis",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Dothideomycetes',
          order='Capnodiales',
          family='Mycosphaerellaceae'),

     list(species="Ruditapes philippinarum",
          kingdom='Animalia',
          phylum='Mollusca',
          class='Bivalvia',
          order='Venerida',
          family='Veneridae'),

     list(species="Saccharina japonica",
          kingdom='Chromista',
          phylum='Ochrophyta',
          class='Phaeophyceae',
          order='Laminariales',
          family='Laminariaceae'),

     list(species="Tripylus abatoides",
          kingdom='Animalia',
          phylum='Echinodermata',
          class='Echinoidea',
          order='Spatangoida',
          family='Prenasteridae'),

     list(species="Trypanosoma brucei",
          kingdom='Excavata',
          phylum='Euglenozoa',
          class='Kinetoplastea',
          order='Trypanosomatida',
          family='Trypanosomatidae'),

     list(species="Undaria pinnatifida",
          kingdom='Chromista',
          phylum='Ochrophyta',
          class='Phaeophyceae',
          order='Laminariales',
          family='Alariaceae'),

     list(species="Zymoseptoria tritici",
          kingdom='Fungi',
          phylum='Ascomycota',
          class='Dothideomycetes',
          order='Capnodiales',
          family='Mycosphaerellaceae'))

