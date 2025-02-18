Results from GBA model validation

figures/

* plot_allocation.html - predicted protein allocation from the model P17
* plot_proteomics_ecocyc.html - plots of 3 experimental proteomics datasets
* validation_plots.html - other validation datasets (fraction of degraded rRNA, RNAP allocation, biomass composition)
* GBA_Model_P17_1parameter_test.pdf - model sensitivity to different parameters
* GBA_Model_P17_1ribocomp_test.pdf - predictions of growth rate at different ribosome compositions and different parameters

code/

Scripts to produce figures in the folder figures
* plot_allocation.R - plot proteome allocation predicted by GBA
* plot_proteomics_ecocyc.R - plot experimental proteomics data
* validation_plots.R - plot experimental data

code/Results GBA
* GBA Model P17old, mean time (25.5s) results.csv - results from GBA

code/Models
* P17old - model with 17 reactions
* M19.ods - model with 19 reactions
* R2deg.ods - model with 7 reactions

data/

Experimental data from E. coli
* biomass.csv - biomass composition data from
* bremer_percent_transcription.csv - RNA allocation at different growth rates from Bremer & Dennis 2008
* gausing_RNA_deg.csv - fraction of degraded rRNA at different growth rates
* EV2-Samples-1.csv, EV3-Samples-2.csv, EV8-AbsoluteMassFractions-1.csv, EV9-AbsoluteMassFractions-2.csv - proteomics data from Mori 2021
* Schmidt_proteomics_table_s23.csv, Schmidt_proteomics_table_s6.csv - proteomics data from Schmidt 2016
* p17_groups_detailed.txt - EcoCyc annoation of proteome
