# function maskGeneOutliers
# takes in a matrix of TPMs, rows=samples , cols=genes
# outputs files with outliers masked as Na values
# AFB 16/01/2020

def maskGeneOutliers(file_dir="", pattern=".csv", filesep=",", outdir=None, outlabel=""):
  
  # """
  # file_dir: directory containing the rpkms
  # pattern: pattern to search for, default is all .csv files
  # outdir: where should the files go, default is file_dir
  # filesep: file separating char, default is comma
  # outlabel: add label to output files
  # """
  
  import os 
  import pandas as pd
  from scipy import stats
  import numpy as np
  import re
  import multiprocessing
  from joblib import Parallel, delayed
  
  if outdir == None:
    outdir = file_dir
  
  # get all file paths
  file_paths = [file for file in os.listdir(file_dir) if pattern in file]
  print(file_paths)
  
  num_cores = len(file_paths)
  
  def dofilesInParallel(__file__):
    
    print(__file__)
    
    # import df
    df = pd.read_csv(open(file_dir+__file__), sep=filesep, encoding="utf-8", engine='python', index_col=0, header=0) # genes cols, samples rows
    
    # genes to cols 
    if len(df.index.values) > len(df.columns.values):
      df = df.transpose()
    
    if 'brain_region' in df.columns.values:
      df.drop('brain_region', axis='columns', inplace=True)
    
    def doMaskOutliers(gene):
    
      gene_clean = gene[np.logical_not(np.isnan(gene))]
      
      # get gene iqr
      iqr = stats.iqr(gene_clean) # inter-quartile range
      
      # upper and lower outlier cut-offs - quartile +/- 3*iqr
      lower_lim = np.quantile(gene_clean, 0.25) - (iqr*3) # lower outlier cut-off
      upper_lim = np.quantile(gene_clean, 0.75) + (iqr*3) # upper outlier cut-off
      
      # masking outliers with NaN
      gene.values[gene.values > upper_lim] = np.nan
      gene.values[gene.values < lower_lim] = np.nan
      
      # returning outlier-masked row
      return(gene)
      
    # applying across genes
    df = df.apply(doMaskOutliers, axis=0)
    
    print("total NaN: ", df.isna().sum().sum())
    
    # write out
    os.chdir(outdir)
    fname = __file__.replace(pattern, '')
    df.to_csv(fname + "_" + outlabel + "_masked_outliers.csv", index=True, header=True)

  # run files in parallel
  Parallel(n_jobs=num_cores)(delayed(dofilesInParallel)(__file__=i) for i in file_paths)

    
