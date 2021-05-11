# NOTE:
# There is one important dependency not in here: if the raw data 
# the body sizes changes, one needs to delete the cached popsize_chains.Rdata
# object. Remove combined_data.tsv and re-run make if the raw size data changes

all: combined_data.tsv post_pred_draws.Rdata stapley_et_al_2017_updated.tsv leffler_et_al_2012_updated.tsv stapley_et_al_2017_otl.Rdata romiguier_et_al_2014_updated.tsv

## Corbett-Detig et al (2015) data

# CD2015 diversity data 
CD2015_div=S2_Data.xls

# CD2015 map lengths 
CD2015_maplens=journal.pbio.1002112.s016.XLS

# CD2015 genome-size
CD2015_genomsize=journal.pbio.1002112.s018.XLS

# CD2015 selfing data
CD2015_selfing=journal.pbio.1002112.s017.XLS

CD2015_DATA=$(CD2015_div) $(CD2015_maplens) $(CD2015_genomsize) $(CD2015_selfing) 

corbett_detig_2015_updated.tsv: $(CD2015_DATA) corbett_detig_2015_data.r
	Rscript corbett_detig_2015_data.r

## Stapley et al (2017) data
S2017_DATA=stapley_et_al_2017.csv

stapley_et_al_2017_updated.tsv: $(S2017_DATA) stapley_et_al_2017_data.r
		Rscript stapley_et_al_2017_data.r

## Leffler et al (2012) Data
leffler_et_al_2012_updated.tsv leffler_et_al_2012_hier.Rdata: leffler_et_al_2012.tsv leffler_et_al_2012_data.r
	Rscript leffler_et_al_2012_data.r

## Romiguier et al (2014) Data
romiguier_et_al_2014_updated.tsv: romiguier_et_al_2014_supptable_2.xls romiguier_et_al_2014_data.r
	Rscript romiguier_et_al_2014_data.r

gbif_ranges.tsv: leffler_et_al_2012_updated.tsv stapley_et_al_2017_updated.tsv romiguier_et_al_2014_updated.tsv species_ranges.r ../R/range_funcs.r species_range_fixes.r 
	Rscript species_ranges.r

## Combined Dataset
combined_data.tsv combined_dataset_hier.Rdata calibrated_phylogeny.Rdata post_pred_draws.Rdata: stapley_et_al_2017_updated.tsv leffler_et_al_2012_updated.tsv corbett_detig_2015_updated.tsv gbif_ranges.tsv combined_dataset.r taxonomy_fixes.r
	Rscript combined_dataset.r


clean-cache:
	rm -f popsize_chains.Rdata post_pred_draws.Rdata gbif_ranges.tsv

clean-cache-all:
	rm -f popsize_chains.Rdata combined_tree.Rdata
	rm -f datelife.Rdata combined_data_hier.Rdata 
	rm -f taxadb_species_ids.Rdata
	rm -f taxadb_genera_ids.Rdata
	# rm -f redlist.Rdata  ## this almost never needs to be cleaned

clean:
	rm -f combined_data.tsv leffler_et_al_2012_updated.tsv
	rm -f combined_data_hier.Rdata leffler_et_al_2012_hier.Rdata
