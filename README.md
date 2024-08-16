Tempo and mode of diapause evolution in butterflies
Sridhar Halali, Etka Yapar, Christopher W Wheat, Niklas Wahlberg, Karl Gotthard, Nicolas Chazot, Sören Nylin, Philipp Lehmann


Dataset and R codes for modelling diapause evolution in butterflies

Summary of the contents

## Trait data
- `diapause_classification_data_RAW.xlsx`: This file contains contains raw data on five, three and binary diapause classification along with the data on tree id (as in the phylogeny), species genus and references used for diapause classification. Note that the full citations corresponding to these references are provided in the supplementary information.
- `diapause_classification_data_PROCESSED.xlsx`: This file contains the data on diapause classification but has been processed (e.g. by removing NA rows) for running phylogenetic comparative analyses. 


## Phylogenies
- `S_3b_core_analysis_median_ages.tre`: Global genus-level phylogeny obtained from from Chazot et al. (2019):10.1093/sysbio/syz002 
- `butterfly_phylogeny_pruned.trees`: This is a pruned phylogeny consisting of 952 tips ready for running the comparative analyses. 

## Code
- `Preprocessing data & some plots.Rmd`: Code for processing raw data (e.g. renaming some columns/diapause classification), pruning phylogeny and making some basic plots 
- `Fitting Mk models, simmap & ancestral state estimation_phytools.Rmd`: Code for fitting the Mk models, extracting transition rates, carrying out stochastic mapping and ancestral state estimation using both maximum likelihood and stochastic mapping in phytools R package
- `Counting transitions & rate through time.Rmd`: Code for counting number of state transition (using 1000 stochastic maps) and calculating rate through time 
- `Fitting hidden rate models_corHMM.Rmd`: Code for fitting several hidden rate models using the corHMM R package. 
- `Plotting uncertainty in ancestral estimation.Rmd`: Code for calculating the uncertainty in ancestral state estimation across methods and root priors  using the sum of squares metric


