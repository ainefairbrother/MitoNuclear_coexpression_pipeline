################# pipeline to process GTEx data #################
## Author: Aine S Fairbrother-Browne                            #
## Date: 01/2020                                                #
## Description: pipleline to process GTEx TPM files and output  #
# nuc-mito correlation matrices                                 #
## Steps (functions):                                           #
# 1. Preprocess raw TPM file (preprocessCountTPMFiles.py)       #
# 2. Filtration of genes to remove unexpressed genes, and       #
# samples with no expression values (filterNullGenesAndSamps)   #
# 3. Log10 & median normalisation (log10MedNormalise)           #
# 4. Mask outliers (maskOutliers)                               #
# 5  Regress out covariates (gtex_regress_covariates)           #
# 6. Generate nuclear-mitochondrial correlation and             #
# corresponding pvalue matrices (genCorrs)                      #
# Run with: python3.7 implementing_whole_pipeline.py            #
#################################################################

################# add paths of fns to sys path ##################
import sys
sys.path.insert(0,'/home/abrowne/Scripts/data_processing/')
sys.path.insert(0,'/home/abrowne/projects/GTEx_6p/scripts/')
sys.path.insert(0,'/home/abrowne/projects/recount_PD_case_control/scripts/')
sys.path.insert(0,'/home/abrowne/projects/ROSMAP/scripts/')

################# def target dir containing files ###############

pwd = "/home/abrowne/projects/ROSMAP/data/celltype_prop_corrected/"

################# import pipeline functions #####################

from preprocessCountTPMFiles import *
from filterGenes import *
from RPKM_to_TPM import *
from fpkm_to_TPM import *
from counts2TPM import *
from splitCaseCtrl import *
from maskOutliers import *
from cov_regression_gtex import *
from cov_regression_rosmap import *
from genCorrs import *
from genCorrSummaryTable import *

################# pipeline ######################################

################# clean counts_all_tpm files ####################
# 
# print('Pre-processing count_all_TPM files')
# preprocessCountTPMFiles(
#   file_dir=pwd,
#   pattern=".txt",
#   outdir=None,
#   outlabel="",
#   sra_table_path="/home/abrowne/projects/GTEx_6p/gtex_metadata/SraRunTable_RNAseq.txt")

################# count/fpkm/rpkm --> tpm #######################

# # convert counts to TPMs
# print("converting counts to TPMs")
# counts2TPM(
#   file_dir=pwd,
#   pattern="_0filtered.csv",
#   filesep=",",
#   outdir=None,
#   outlabel="")
# # # out label = _TPM.csv

# # convert RPKMs to TPMs
# print("converting RPKMs to TPMs")
# rpkm_to_tpm(
#   rpkm_dir = pwd,
#   out_dir = None,
#   pattern = "_0filtered.csv")
# # out label = _TPM.csv

# # convert fpkms to TPMs
# print("converting FPKMs to TPMs")
# fpkm_to_tpm(
#   file_dir=pwd,
#   out_dir=None,
#   pattern="_preprocessed.csv",
#   outfile_label="")
# # out label = _TPM.csv

################# filter for genes with TPM > 0 #################

# # filter for set of genes with RPKM/TPM/count>0 in all samples
# print("filtering genes")
# filterNullGenesAndSamps(
#   file_dir=pwd,
#   pattern="_preprocessed.csv",
#   filesep=",",
#   outdir=None,
#   outlabel="")
# # out label = _0filtered.csv

##################### log10 median norm #########################

# #importing R script log10MedNormalise.R
# print("log10 median normalising")
# import rpy2.robjects as robjects # to allow calling R script
# r_source = robjects.r['source']
# r_source("""/home/abrowne/Scripts/data_processing/log10MedNormalise.R""")
# log10MedNormalise = robjects.globalenv["log10MedNormalise"]
# log10MedNormalise(
#   file_dir=pwd,
#   pattern="*_0filtered.csv",
#   filesep=",",
#   outlabel="",
#   med_norm="TRUE")
# # # out label = log10_mediannorm_TPM.csv

################# case/control split ############################

# # split cases and controls - recount2
# splitCaseCtrl(
#   file_dir=pwd,
#   metadata_dir= "/home/abrowne/recount_PD_case_control/metadata/SRP058181_sample_metadata.csv",
#   pattern="masked_outliers.csv",
#   filesep=",",
#   outdir=None)

########### filter for gene set in all files ####################

# # filter for same gene set in all regions
# filterGenesInAllFiles(
#   file_dir=pwd,
#   pattern="_masked_outliers.csv",
#   filesep=",",
#   outdir=None,
#   outlabel="",
#   med_fill_na=False)
# # out label = _genefiltered.csv

################# mask outliers #################################

# # mask outliers
# print("masking outliers")
# maskGeneOutliers(
#   file_dir=pwd,
#   pattern="log10_mediannorm_TPM.csv",
#   filesep=",",
#   outdir=None,
#   outlabel="")
# # out label = _masked_outliers.csv

################# regress out covs ##############################

# # regress out harccoded covariates ROSMAP
# rosmap_regress_covariates(
#   tpm_dir=pwd,
#   pattern=".csv",
#   meta_path="/home/abrowne/projects/ROSMAP/metadata/ROSMAP_meta_preprocessed.csv",
#   out_dir=None,
#   outlabel=""
#   )

# # regress out harccoded covariates GTEx
# print("regressing out covariates")
# gtex_regress_covariates(
#   tpm_dir = pwd,
#   pattern = "masked_outliers.csv",
#   meta_path = "/home/abrowne/projects/GTEx_6p/gtex_metadata/GTEx_Data_V6_Annotations_SampleAttributesDS.txt",
#   pheno_path = "/home/abrowne/projects/GTEx_6p/gtex_metadata/pheno6p.txt",
#   out_dir = None,
#   outlabel = ""
#   )
# # out label = _residuals.csv

################# mt-nuc corrs ##################################

# generate mt-nuc correlations
print("generating mt-nuc correlations")
genCorrs(
  file_dir=pwd,
  pattern="_residuals_celltypeprop_corrected.csv",
  outdir=None,
  method="spearman",
  outlabel="",
  random_shuffle_cols = False,
  all_corrs = False
  )
# out label = _spearman_corrs.csv, _spearman_pvals.csv

################# generate summary table ########################

# # generate summary_table for correlations + pvalues
# print("generating summary table")
# genCorrSummaryTable(
#   file_dir=pwd,
#   outdir=None,
#   pattern="spearman",
#   outfile_label=""
#   )






