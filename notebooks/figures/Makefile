all: figures supp_figures

MSFIGDIR = ../../manuscript/figures

figures: $(MSFIGDIR)/diversity_popsize_full.pdf $(MSFIGDIR)/diversity_pcm.pdf $(MSFIGDIR)/figure_3.pdf $(MSFIGDIR)/Nc_figure.pdf

supp_figures: $(MSFIGDIR)/range_categories.pdf $(MSFIGDIR)/iucn_redlist_ranges.pdf \
	            $(MSFIGDIR)/damuth.pdf $(MSFIGDIR)/body_mass_length.pdf \
							$(MSFIGDIR)/biconmap.pdf \
							$(MSFIGDIR)/diversity_popsize_averages.pdf \
							$(MSFIGDIR)/genome_size_popsize.pdf \
							$(MSFIGDIR)/genome_size_maplength.pdf \
              $(MSFIGDIR)/diversity_popsize_linear.pdf \
              $(MSFIGDIR)/corbett_detig.pdf \
              $(MSFIGDIR)/biconmap_diversity.pdf \
              $(MSFIGDIR)/div_popsize_post_phyla.pdf \
              $(MSFIGDIR)/diversity_range_bodymass.pdf \
              $(MSFIGDIR)/diversity_popsize_redlist.pdf \
              $(MSFIGDIR)/Nc_figure.pdf \
              $(MSFIGDIR)/figure_3_alt.pdf \
              $(MSFIGDIR)/linked_sel_params.pdf \
              $(MSFIGDIR)/species_counts.pdf \
              $(MSFIGDIR)/diversity_leffler_ranges.pdf

clean_repel:
	rm -f repel_genome_size_popsize.tsv repel_genome_size_maplength.tsv repel_full.tsv repel.tsv

clean: clean_repel
	rm -f range_categories.pdf iucn_redlist_ranges.pdf figure_3.pdf figure_3_alt.pdf damuth.pdf diversity_pcm.pdf
	rm -f body_mass_length.pdf biconmap.pdf diversity_popsize_full.pdf 
	rm -f diversity_popsize.pdf genome_size_popsize.pdf  genome_size_maplength.pdf

# components of diversity_pcm.pdf
# needs to be updated in Illustrator too
node_heights.pdf: ../../data/nodeheights.Rdata color_scheme.r node_heights.r  ../../data/main_datasets.Rdata
	Rscript node_heights.r

Nc_figure.pdf: ../../data/combined_data.tsv diversity_popsize_full.r color_scheme.r ../../data/diversity_popsize_chains.Rdata Nc_figure.r
	Rscript Nc_figure.r

popsize_tree.pdf: ../../data/nodeheights.Rdata color_scheme.r popsize_tree.r
	Rscript popsize_tree.r

diversity_pcm.pdf: node_heights.pdf popsize_tree.pdf div_popsize_post.pdf
	# dumby rule
	echo "You need to run Illustator!"

diversity_leffler_ranges.pdf: species_ranges.r ../../data/combined_data.tsv diversity_leffler_ranges.r
	Rscript diversity_leffler_ranges.r


range_categories.pdf: species_ranges.r ../../data/combined_data.tsv species_ranges.r
	Rscript species_ranges.r

iucn_redlist_ranges.pdf: redlist_ranges.r ../../data/redlist_popsizes.tsv
	Rscript redlist_ranges.r

figure_3_alt.pdf: figure_3_alt.r ../../data/linked_sel_div_changes.Rdata
	Rscript figure_3_alt.r


figure_3.pdf: figure_3.r ../../data/linked_sel_div_changes.Rdata
	Rscript figure_3.r

damuth.pdf: damuth.r ../../data/damuth1987.tsv
	Rscript damuth.r

body_mass_length.pdf: body_mass_length.r ../../data/romiguier_et_al_2014_updated.tsv ../../data/popsize_chains.Rdata
	Rscript body_mass_length.r

biconmap.pdf: biconmap.r ../../data/biconmap_data.Rdata
	Rscript biconmap.r

diversity_popsize_linear.pdf: ../../data/combined_data.tsv diversity_popsize.r color_scheme.r ../../data/diversity_popsize_chains.Rdata
	Rscript diversity_popsize_linear.r

diversity_popsize_redlist.pdf: ../../data/combined_data.tsv diversity_popsize_redlist.r color_scheme.r ../../data/diversity_popsize_chains.Rdata
	Rscript diversity_popsize_redlist.r

diversity_popsize_full.pdf: ../../data/combined_data.tsv diversity_popsize_full.r color_scheme.r ../../data/diversity_popsize_chains.Rdata
	Rscript diversity_popsize_full.r

diversity_popsize_averages.pdf: ../../data/combined_data.tsv diversity_popsize_averages.r color_scheme.r ../../data/diversity_popsize_chains.Rdata
	Rscript diversity_popsize_averages.r

genome_size_maplength.pdf: ../../data/combined_data.tsv genome_size_maplength.r color_scheme.r
	Rscript genome_size_maplength.r

genome_size_popsize.pdf: ../../data/combined_data.tsv genome_size_popsize.r color_scheme.r
	Rscript genome_size_popsize.r

corbett_detig.pdf: ../../data/corbett_detig_data.Rdata color_scheme.r corbett_detig_fig.r
	Rscript corbett_detig_fig.r

div_popsize_post.pdf: ../../data/diversity_popsize_chains.Rdata color_scheme.r div_popsize_post.r ../../data/combined_data.tsv
	Rscript div_popsize_post.r

div_popsize_post_phyla.pdf: ../../data/phyla_diversity_popsize_chains.Rdata color_scheme.r div_popsize_post_phyla.r 
	Rscript div_popsize_post_phyla.r


biconmap_diversity.pdf: ../../data/biconmap_data.Rdata color_scheme.r biconmap_diversity.r
	Rscript biconmap_diversity.r

diversity_range_bodymass.pdf: diversity_range_bodymass.r ../../data/main_datasets.Rdata ../../data/diversity_range_bodymass_chains.Rdata color_scheme.r
	Rscript diversity_range_bodymass.r

linked_sel_params.pdf: ../../data/sweep_param_ests.tsv ../../data/shapiro_et_al_2007.tsv linked_sel_params.r
	Rscript linked_sel_params.r

species_counts.pdf: ../../data/species_by_class.tsv species_counts.r
	Rscript species_counts.r



#../../manuscript/figures/diversity_body_size.pdf 
# ../../manuscript/figures/map_length_body_size.pdf ../../manuscript/figures/map_length_body_size_bayesian.pdf ../../manuscript/figures/fecundity_body_size.pdf 

# ../../manuscript/figures/diversity_body_size.pdf: diversity_body_size.pdf
	# cp diversity_body_size.pdf ../../manuscript/figures/

$(MSFIGDIR)/diversity_popsize_averages.pdf: diversity_popsize_averages.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/diversity_popsize_redlist.pdf: diversity_popsize_redlist.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/diversity_popsize_full.pdf: diversity_popsize_full.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/iucn_redlist_ranges.pdf: iucn_redlist_ranges.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/diversity_leffler_ranges.pdf: diversity_leffler_ranges.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/range_categories.pdf: range_categories.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/damuth.pdf: damuth.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/body_mass_length.pdf: body_mass_length.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/biconmap.pdf: biconmap.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/figure_3_alt.pdf: figure_3_alt.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/figure_3.pdf: figure_3.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/genome_size_popsize.pdf: genome_size_popsize.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/genome_size_maplength.pdf: genome_size_maplength.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/diversity_popsize_linear.pdf: diversity_popsize_linear.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/diversity_pcm.pdf: diversity_pcm.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/corbett_detig.pdf: corbett_detig.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/biconmap_diversity.pdf: biconmap_diversity.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/div_popsize_post_phyla.pdf: div_popsize_post_phyla.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/diversity_range_bodymass.pdf: diversity_range_bodymass.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/Nc_figure.pdf: Nc_figure.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/linked_sel_params.pdf: linked_sel_params.pdf
	cp $< $(MSFIGDIR)

$(MSFIGDIR)/species_counts.pdf: species_counts.pdf
	cp $< $(MSFIGDIR)



# ../../manuscript/figures/map_length_body_size.pdf: map_length_body_size.pdf ../../data/body_size_map_length.Rdata
# 	cp map_length_body_size.pdf ../../manuscript/figures/

# ../../manuscript/figures/map_length_body_size_bayesian.pdf: map_length_body_size_bayesian.pdf 

# ../../manuscript/figures/fecundity_body_size.pdf: fecundity_body_size.pdf
# 	cp fecundity_body_size.pdf ../../manuscript/figures/

# fecundity_body_size.pdf: fecundity_body_size.r ../../data/body_size_fecundity.Rdata ../../data/romiguier_et_al_2014_updated.tsv
# 	Rscript fecundity_body_size.r

# map_length_body_size.pdf: map_length_body_size.r plot_funs.r  ../../data/body_size_map_length.Rdata ../../data/diversity_size_theory.Rdata
# 	Rscript map_length_body_size.r

# map_length_body_size_bayesian.pdf: map_length_body_size_bayesian.r plot_funs.r ../../data/body_size_map_length_bayesian.Rdata
# 	Rscript map_length_body_size_bayesian.r

# diversity_body_size.pdf: diversity_body_size.r plot_funs.r ../../data/leffler_et_al_2012_updated.tsv 
# 	Rscript diversity_body_size.r
