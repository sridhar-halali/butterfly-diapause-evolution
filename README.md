Tempo and mode of diapause evolution in butterflies
Sridhar Halali, Etka Yapar, Christopher W Wheat, Niklas Wahlberg, Karl Gotthard, Nicolas Chazot, Sören Nylin, Philipp Lehmann


Dataset and R codes for modelling diapause evolution in butterflies

Summary of the contents

## Trait data
- `diapause_classification_data_RAW.xlsx`: This file contains raw data for five, three and binary diapause classification along with the data on tree id (as in the phylogeny), species, genus, and references used for diapause classification. Note that the full citations corresponding to these references are provided in the supplementary information. All the columns are self-explanatory and contain information on the tip labels in the phylogeny, genus, species and species family followed by five, three and two-state diapause classification. Note that there are ’NA’ rows in the ‘tree_species’ column which simply means that the species name for that specific genus is not known. 
- `diapause_classification_data_PROCESSED.xlsx`: This file contains the data on diapause classification but has been processed (e.g. by removing ‘skip’ and ‘skip_high_altitude’ rows from ‘diapause_classification_data_RAW.xlsx’ file) for phylogenetic comparative analyses. 


## Phylogenies
- `S_3b_core_analysis_median_ages.tre`: Global genus-level phylogeny obtained from Chazot et al. [1]
- `butterfly_phylogeny_pruned.trees`: This is a pruned phylogeny (from the global phylogeny from Chazot et al [1]) consisting of 952 tips ready for running the comparative analyses.
- `ButterflyPhyloMulti_pruned.trees`: This is a pruned phylogeny with 952 tips consisting of 200 posterior trees from Chazot et al. [1]
## R Scripts
- `Preprocessing data & some plots.Rmd`: Code for processing raw data (e.g. renaming some columns/diapause classification), pruning phylogeny and making some basic plots 
- `Fitting Mk models, simmap & ancestral state estimation_phytools.Rmd`: Code for fitting the Mk models, extracting transition rates, carrying out stochastic mapping and ancestral state estimation using both maximum likelihood and stochastic mapping in phytools R package
- `Counting transitions & rate through time.Rmd`: Code for counting number of state transition (using 1000 stochastic maps) and calculating rate through time 
- `Fitting hidden rate models_corHMM.Rmd`: Code for fitting several hidden rate models using the corHMM R package. 
- `Plotting uncertainty in ancestral estimation.Rmd`: Code for calculating the uncertainty in ancestral state estimation across methods and root priors  using the sum of squares metric
- `source_rate_through_time_Hughes_et_al_2021.R`: R code provided by Hughes et al. [2] (code available here: https://github.com/jakeberv/mammal_arboreality) for calculating rate through time for both single and using posterior trees.  
- `getRates_modified_Etka_Yapar.R`: modified function from Hughes et al. [2], mainly to plot 95% CI on points, when stochastic maps are generates for a single best (e.g. consensus) tree.  


[1] Chazot, N., Wahlberg, N., Freitas, A. V. L., Mitter, C., Labandeira, C., Sohn, J. C., ... & Heikkilä, M. (2019). Priors and posteriors in Bayesian timing of divergence analyses: the age of butterflies revisited. Systematic biology, 68(5), 797-813.

[2] Hughes, J. J., Berv, J. S., Chester, S. G., Sargis, E. J., & Field, D. J. (2021). Ecological selectivity and the evolution of mammalian substrate preference across the K–Pg boundary. Ecology and Evolution, 11(21), 14540-14554.
