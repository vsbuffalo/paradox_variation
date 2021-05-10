
main_analysis.html: main_analysis.Rmd ../R/utilities.r ../data/combined_data.tsv ../data/calibrated_phylogeny.Rdata ../data/post_pred_draws.Rdata
	Rscript render_notebook.r

clean-cache:
	rm -rf ../data/phylo_mm_lognormal_chains.Rdata
	rm -rf ../data/phylo_mm_regression_chains.Rdata
	rm -rf ../data/phylo_mm_regression_by_phyla_chainsdf.Rdata
	rm -rf ../data/phylo_mm_regression_phyla_alt_chains.Rdata
	rm -rf ../data/phylo_mm_regression_body_size_chains.Rdata
	rm -rf ../data/phylo_mm_regression_range_chains.Rdata
