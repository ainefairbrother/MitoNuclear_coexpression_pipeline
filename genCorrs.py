#!/usr/bin python

# function genSpearmanCorrs
# takes in a matrix of TPMs, rows=samples , cols=genes
# generates all gene pairwise spearman corrs
# outputs filtered files
# AFB 16/01/2020


def genCorrs(file_dir="", pattern=".csv", outdir=None, outlabel="", method="pearson", random_shuffle_cols=False, all_corrs=False):
  
  # """
  # file_dir: directory containing the rpkms
  # pattern: pattern to search for, default is all .csv files
  # outdir: where should the files go, default is file_dir
  # outlabel: add label to output files
  # method: correlation method - spearman or pearson 
  # random_shuffle_cols: creates random mt-nuc associations that assign a random nuc gene to the r statistic. Useful only for v specific
  # analyses, normally should remain False
  # all_corrs: performs pairwise correlations between all gene pairs as opposed to all gene pairs * 13 mt genes. Takes considerably longer. 
  # Not useful for generation of a simple nuc-mt correlation matrix. 
  # """

  # import libs
  import os
  import numpy as np
  import pandas as pd
  from scipy import stats
  import re
  import multiprocessing
  from joblib import Parallel, delayed

  if outdir == None:
    outdir = file_dir

  # get all file paths
  file_paths = [file for file in os.listdir(file_dir) if pattern in file]
  print(file_paths)
  
  # define num of cores, here it's one core per file
  num_cores = len(file_paths)
  
  # defining mito genes - allows generation of a nuc x mito corr matrix
  mito_genes = ['ENSG00000198888','ENSG00000198763','ENSG00000198840','ENSG00000198886',
  'ENSG00000212907','ENSG00000198786','ENSG00000198695','ENSG00000198899',
  'ENSG00000228253','ENSG00000198804','ENSG00000198712','ENSG00000198938',
  'ENSG00000198727']

  def dofilesInParallel(__file__):
    
    fn = __file__.replace(pattern, '')

    print(fn)

    # import df
    df = pd.read_csv(open(file_dir+__file__), sep=',', encoding="utf-8", engine='python', index_col=0, header=0) # genes cols, samples rows
  
    if len(df.columns.values) < len(df.index): # 
      df = df.T
      
    print(df.head())
      
    ## if random shuffle ###
    if random_shuffle_cols == True:
      import random
      cols = df.columns.tolist()
      random.shuffle(cols)
      df.columns = cols
    ########################
      
    df.columns = df.columns.str.replace("\\..*", "")
    df.index = df.index.str.replace("\\..*", "")
    
    # choice about whether to generate all correlations possible, or to generate mt-nuc only. The latter is considerably quicker
    if all_corrs == False:
      corr_matrix_rows = mito_genes
    elif all_corrs == True:
      corr_matrix_rows = df.columns
    
    if method == "pearson":
      # generate pairwise corrs
      df_corr = pd.DataFrame() # Correlation matrix
      df_p = pd.DataFrame()  # Matrix of p-values
      counter=0
      for x in corr_matrix_rows:
        for y in df.columns:
          print(counter)
          
          # get gene value lists
          x_ls = df[x]
          y_ls = df[y]
          
          # list of indices where NaNs are in both
          nan_in_both = sorted(list(set(list(np.where(np.isnan(x_ls))[0]) + list(np.where(np.isnan(y_ls))[0]))), reverse=True)
          
          # remove these elements of both 
          x_ls = [i for j, i in enumerate(x_ls) if j not in nan_in_both]
          y_ls = [i for j, i in enumerate(y_ls) if j not in nan_in_both] 
          
          # perform corr on cleaned value lists
          corr = stats.pearsonr(x_ls, y_ls) # calculate the spearman r of col x and col y
          df_corr.loc[x,y] = corr[0] # assign corrs to df_corr
          df_p.loc[x,y] = corr[1] # assign pvals to df_p
          counter+=1
          
    elif method == "spearman":
      # generate pairwise corrs
      df_corr = pd.DataFrame() # Correlation matrix
      df_p = pd.DataFrame()  # Matrix of p-values
      counter=0
      for x in corr_matrix_rows:
        for y in df.columns:
          print(counter)

          # get gene value lists
          x_ls = df[x]
          y_ls = df[y]
          
          # list of indices where NaNs are in both
          nan_in_both = sorted(list(set(list(np.where(np.isnan(x_ls))[0]) + list(np.where(np.isnan(y_ls))[0]))), reverse=True)
          
          # remove these elements of both 
          x_ls = [i for j, i in enumerate(x_ls) if j not in nan_in_both]
          y_ls = [i for j, i in enumerate(y_ls) if j not in nan_in_both] 
          
          # perform corr on cleaned value lists
          corr = stats.spearmanr(x_ls, y_ls) # calculate the spearman r of col x and col y
          df_corr.loc[x,y] = corr[0] # assign corrs to df_corr
          df_p.loc[x,y] = corr[1] # assign pvals to df_p
          counter+=1
    
    print(df_corr.shape)
    print(df_p.shape)
    
    # write out
    os.chdir(outdir)
    df_corr.to_csv(fn + "_" + outlabel +"_" + method + "_corrs.csv", index=True, header=True)
    df_p.to_csv(fn + "_" + outlabel +"_" + method + "_pvals.csv", index=True, header=True)
    
  Parallel(n_jobs=num_cores)(delayed(dofilesInParallel)(__file__=i) for i in file_paths)
    

