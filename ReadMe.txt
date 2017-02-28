A) GENERAL SCRIPTS

randNum_seeds.R
Sets up a list of random numbers to set the seed when parallelizing scripts.


B) DSPR SPECIFIC SCRIPTS

1) Simulate Phenotypes 
GetPhenotypes_DSPR.R and GetPhenotypes_multi_DSPR.R
Inputs:
pARILs.rda
snpfreqs_RILS_A_chromosome.txt
SNPtable_A_RILs_inferred_Release4_chromosome.txt
A_RILLIST_release4_chromosome.txt

Output:
pheno_biallelic_",100*eff,".rda
pheno_multi3_gene_",100*eff,".rda

2) Make Kinship Matricies for each chromosome
kinshipMatrix_Genetic_A_LOCO.R
Uses MySQL database 
output Kin_LOCO_A.rda

3) Set up haplotypes
makebiglist_DSPR.R

input: Data/DSPR/SimPhenos/pheno_biallelic_1.rda
../GenomeCache/Release2/A/

output:
Data/DSPR/pA_biglist.rda
Data/DSPR/pA_bigarray.rda

4) Do mapping
mapping_DSPR.R
mapping_DSPR_multi.R

input:
pA_bigarray.rda
Kin_LOCO_A.rda
pheno_biallelic_",100*eff,".rda
paste("Data/DSPR/SimPhenos/pheno_multi3_gene_",100*eff,".rda",sep="")

output: 
LODseed_biallelic_KLOCO_E",100*eff,"_S",samp.size,".rda
paste("Data/DSPR/Obs_LODs/LODs_multi3gene_KLOCO_E",100*eff,"_S",samp.size,".rda",sep="")

5) Collate results
combine_map_results_LOCOseed.R
  input:
  paste("Data/DSPR/SimPhenos/pheno_biallelic_",100*eff,".rda",sep="")
  paste("Data/DSPR/Obs_LODs/LODseed_biallelic_KLOCO_E",100*eff,"_S",samp.size,".rda",sep="")
  output: 
  Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda
  
combine_map_results_multi_LOCO.R
  input:
  paste("Data/DSPR/SimPhenos/pheno_multi3_gene_",100*eff,".rda",sep="")
  paste("Data/DSPR/Obs_LODs/LODs_multi3gene_KLOCO_E",100*eff,"_S",samp.size,".rda",sep="")
  output: 
  Data/DSPR/Obs_LODs/map_combined_multiLOCO.rda
  Data/DSPR/SimPhenos/dspr.multi.pos.rda

6) Random Phenotypes
mapping_DSPR_randomPhenos.R
input: 
Data/DSPR/perms_DSPR_info.rda
Data/DSPR/pA_bigarray.rda
Data/DSPR/Kin_LOCO_A.rda
output:
paste("Data/DSPR/Perm_LODs/LODSET_",samp.size,"_",ii,".rda",sep="")

7) Get False Positive Rate 
DSPR_FPR.R
input:
paste('Data/DSPR/Perm_LODs/',i,sep='')
output:
Data/DSPR/FPR_rates_ALL.rda
Data/DSPR/FWER.rda

C) DGRP SPECIFIC SCRIPTS

1) Get Kinship and eliminate lines
Kinship_AllLines_DGRP.R
input:
Data/DGRP/dgrp2.rda
output:
Data/DGRP/dropped_lines.rda

2) Impute genotypes
impute_missing_DGRP.R
input:
Data/DGRP/dgrp2.rda
Data/dropped_lines.rda
output:
Data/DGRP/imputed_genos.rda

3) Get Kinship with LOCO method
input:
Data/dgrp2.rda
Data/dropped_lines.rda
output:
Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda

4) Simulate phenotypes
Get_Phenotypes_DGRP.R
input:
Data/DGRP/imputed_genos.rda
output:
Data/DGRP/SimPhenos/allphenos_noinv.rda

Get_Phenotypes_DGRP_multi.R
Data/DGRP/imputed_genos.rda
Data/genes_Dmel_5.txt
output:
Data/DGRP/SimPhenos/allphenos_noinv_multi.rda

4) Do mapping
mapping_DGRP.R
mapping_DGRP_multi.R
input:
Data/randNum.rda
Data/DGRP/imputed_genos.rda
Data/DGRP/SimPhenos/allphenos_noinv_multi.rda
Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda
output:
paste("Data/DGRP/Obs_LODs/DGRPsetLOCO_",s,"_",st.ind,"_",en.ind,".rda",sep="")
paste("Data/DGRP/Obs_LODs/DGRPsetLOCOmulti_",samp.size,"_",st.ind,"_",en.ind,".rda",sep="")

5) Collate results
input:
paste("Data/DGRP/Obs_LODs/DGRPsetLOCO_",s,"_",st.ind,"_",en.ind,".rda",sep="")
paste("Data/DGRP/Obs_LODs/DGRPsetLOCOmulti_",s,"_",st.ind,"_",en.ind,".rda",sep="")
output:
Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda
Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda

6) Random phenotypes
mapping_DGRP_randomPhenos.R
input:
Data/randNum.rda
Data/DGRP/imputed_genos.rda
Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda
output:
should be: file=paste("Data/DGRP/Perm_LODs/All_lods_randLOCO_",samp.size,"_",i,".rda",sep="")
*large set of files. not included in the zip. re- run script to generate.

7)Get false positive rates etc.
FWER:
DGRP_FWER.R
input:
paste("Data/DGRP/Perm_LODs/",ffs[file.i],sep="")
output:
Data/DGRP/FWER.rda

False Positive Rate:
DGRP_FPR.R
input:
Data/geneticmapLIST.rda
Data/DGRP/imputed_genos.rda
paste("Data/DGRP/Perm_LODs/",i,sep="")
output:
paste('Data/DGRP/FPR_rates',tt,'.rda',sep='')

D) Analysis and common scripts

FWER_Region.R
Get family-wise error rate for mapped region (+/-5Mb). 
input:
Data/DGRP/imputed_genos.rda
paste("Data/DGRP/Perm_LODs/",i,sep="")
paste('Data/DSPR/Perm_LODs/',i,sep='')
output:
Data/DGRP/FWER_region_DGRP.rda
Data/DSPR/FWER_region_DSPR.rda

Compare_freqs_pops.R
compare frequencies in both populations and make histogram
input:
Data/DGRP/imputed_genos.rda
Data/DSPR_Release4/snpfreqs_RILS_A_'chr'.txt
output:
Plots/Freq_compare.pdf

DSPR_K_LOCO_PVE_Plot.R
makes plot for DSPR effects with and without LOCO method
input:
Data/DSPR/Obs_LODs/CompareK/map_combined_v1.rda
Data/DSPR/Obs_LODs/CompareK/map_combined_K.rda
Data/DSPR/Obs_LODs/CompareK/map_combined_KLOCO.rda
output:
Plots/DSPR_K_LOCO.pdf
Plots/DSPR_K_LOCO_AllPOS.pdf

FPR_power_beavis.R
Compares false positive rate, power, and PVE bias for various thresholds for DGRP and DSPR. Makes plot of these results
input:
Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda
Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda
Data/DGRP/FPR_rates_ALL.rda
Data/DSPR/FPR_rates_ALL.rda
output:
Plots/FPR_pow_beav.pdf

PVE_Plots.R
Makes plots for simulated QTL for DGRP and DSPR showing PVE for different effects & sample sizes
input:
Data/DGRP/Obs_LODs/map_combined_LOCO_biallelic.rda
Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda
Data/DSPR/Obs_LODs/map_combined_LOCOseed.rda
Data/DSPR/Obs_LODs/map_combined_multiLOCO.rda
output:
Data/DGRP/Obs_LODs/All_Obs_LODs.rda
Data/DSPR/Obs_LODs/All_Obs_LODs.rda
Plots/DSPR_DGRP_beavis.pdf
Plots/DGRP_ALLeff_AtQTL.pdf

E) VALIDATION 
Get_Phenotypes_Validate_within.R
input:
Data/DSPR/Obs_LODs/All_Obs_LODs.rda
Data/DSPR_Release4/pARILs.rda
paste('Data/DSPR_Release4/A_RILLIST_release4_',chromosomes[ii],'.txt',sep='')
paste('Data/DSPR_Release4/SNPtable_A_RILs_inferred_Release4_',chromosomes[ii],'.txt',sep='')
Data/DSPR/SimPhenos/dspr.multi.pos.rda
Data/DSPR/Obs_LODs/All_Obs_LODs.rda
Data/DGRP/Obs_LODs/All_Obs_LODs.rda
Data/DGRP/imputed_genos.rda
Data/DGRP/Obs_LODs/map_combined_LOCO_multi.rda
output:
Data/DSPR/SimPhenos/DSPR_validate_biallelic.rda
Data/DSPR/SimPhenos/DSPR_validate_multi.rda
Data/DGRP/SimPhenos/DGRP_validate_multi.rda
Data/DGRP/SimPhenos/DGRP_validate_biallelic.rda

Validation_btw_pops.R
Generate phenotypes to validate between populations
input:
Data/DSPR/Obs_LODs/All_Obs_LODs.rda
Data/DGRP/imputed_genos.rda
Data/DSPR_Release4/pARILs.rda
paste('Data/DSPR_Release4/A_RILLIST_release4_',chromosomes[ii],'.txt',sep='')
paste('Data/DSPR_Release4/SNPtable_A_RILs_inferred_Release4_',chromosomes[ii],'.txt',sep='')
Data/DGRP/Obs_LODs/All_Obs_LODs.rda
output:
Data/DGRP/SimPhenos/DGRP_validate_btwpop.rda
Data/DSPR/SimPhenos/DSPR_validate_btwpop.rda

validation_DGRP_biallelic.R
within population validation mapping (biallelic)

input:
Data/DGRP/imputed_genos.rda
Data/DGRP/SimPhenos/DGRP_validate_biallelic.rda
Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda
output:
Data/DGRP/Obs_LODs/Validate_LODs_biallelic_DGRP.rda

validation_DGRP_btw.R
mapping in DGRP for between population validation
input:
Data/DGRP/imputed_genos.rda
Data/DGRP/SimPhenos/DGRP_validate_btwpop.rda
Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda
output:
paste("Data/DGRP/Obs_LODs/Validate_LODs_btw_DGRP_",samps,".rda",sep="")

validation_DGRP_mulit.R
within population validation mapping (multiallelic)
input:
Data/DGRP/imputed_genos.rda
Data/DGRP/SimPhenos/DGRP_validate_multi.rda
Data/DGRP/DGRP_kinshipLOCO_R2_relCut.rda
output:
Data/DGRP/Obs_LODs/Validate_LODs_multi_DGRP.rda

validation_DSPR_biallelic.R
within population validation mapping (biallelic)
input:
Data/DSPR/pA_bigarray.rda
Data/DSPR/Kin_LOCO_A.rda
Data/DSPR/SimPhenos/DSPR_validate_biallelic.rda
output:
Data/DSPR/Obs_LODs/Validate_LODs_biallelic_DSPR.rda

validation_DSPR_btw.R
mapping in DSPR for between pop validation
input:
Data/DSPR/pA_bigarray.rda
Data/DSPR/Kin_LOCO_A.rda
Data/DSPR/SimPhenos/DSPR_validate_btwpop.rda
output:
paste("Data/DSPR/Obs_LODs/Validate_LODs_btw_DSPR_",ss,".rda",sep="")

validation_DSPR_multi.R
within population validation mapping (multiallelic)
input:
Data/DSPR/pA_bigarray.rda
Data/DSPR/Kin_LOCO_A.rda
Data/DSPR/SimPhenos/DSPR_validate_multi.rda
output:
Data/DSPR/Obs_LODs/Validate_LODs_multi_DSPR.rda

Validation_summary.R
Summarize validation results. Inputs are output files from validation scripts above.






