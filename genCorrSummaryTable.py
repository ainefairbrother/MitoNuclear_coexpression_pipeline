# function to read in correlation and pval matrix files of the same size
# makes them into long-form
# places them into a df with gene cols to indicate mt-nuc gene pairs

def genCorrSummaryTable(file_dir="", outdir=None, pattern=".csv", outfile_label=""):
  
  # import libs
  import pandas as pd
  import os
  import re
  
  file_names = [file for file in os.listdir(file_dir) if pattern+'_pvals' in file or pattern+'_corrs' in file]
  
  print(file_names)
  
  if outdir == None:
    outdir = file_dir

  # import corrs + pvals
  corrs = []
  pvals = []
  
  # collect corrs and pvals from separate files
  for i in range(0, len(file_names)):
        
    # get file name
    name_match = re.compile('(\w+)__')
    name = file_names[i]#name_match.match(file_names[i]).group(1)
    
    # import df
    df = pd.read_csv(open(file_dir+file_names[i]), index_col=0, header=0, encoding="utf-8", engine='python')
    
    # check that mt are rows
    if 'ENSG00000198888.2' not in df.index.values:
      print("mito not rows - transposing")
      df = df.transpose()
    
    # store in list
    if 'pvals' in file_names[i]:
      
      # from corr matrix --> stacked table
      df = df.stack().reset_index()
      df.columns = ["mt_gene", "nuc_gene", name+"pval"]
      pvals.append(df)
      
    if 'corrs' in file_names[i]:
        
      # from corr matrix --> stacked table
      df = df.stack().reset_index()
      df.columns = ["mt_gene", "nuc_gene", name+"corr"]
    
      corrs.append(df)
      
  # remove label cols for all but first df
  for i in range(0, len(corrs)):
    if i != 0:
      corrs[i] = corrs[i].iloc[:,[2]]

  for i in range(0, len(pvals)):
      pvals[i] = pvals[i].iloc[:,[2]]
  
    # put all together to form summary df - all R vals for all tissues
  corr_df = pd.concat(corrs, axis=1)
  pval_df = pd.concat(pvals, axis=1)
  summary_df = pd.concat([corr_df, pval_df], axis=1)
    
  print("saving file...")
  os.chdir(outdir)
  summary_df.to_csv(outfile_label+"summary_table.csv", index=True, header=True)
        
    
      










