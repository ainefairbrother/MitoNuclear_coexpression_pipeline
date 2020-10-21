
def rosmap_regress_covariates(tpm_dir, pattern, meta_path, out_dir=None, outlabel=""):  
  
  """
  :param tpm_dir: directory containing the TPM files - genes are cols, samples are rows
  :param pattern: unique pattern to select only TPM files
  :param meta_path: path to metadata file, downloadable from GTEx i.e. GTEx_Data_V6_Annotations_SampleAttributesDS.txt
  :param pheno_path: path to phenotype file, downloadable from GTEx, gives age, sex, COD associated with patient samples
  :param out_dir: directory to save residual files into
  :param outlabel: optional label to add to output files
  :return: residuals df
  
  corrects out the following covars: ['pmi', 'RIN', 'library_batch', 'race', 'msex', 'study', 'age_death', 'age_at_visit_max']
  """
  
  if out_dir == None:
    out_dir = tpm_dir

  # import libs
  import pandas as pd
  import re
  import numpy as np
  from sklearn import preprocessing
  from sklearn.preprocessing import StandardScaler
  from sklearn.linear_model import LinearRegression
  from scipy.stats import shapiro
  import os
  #from tqdm import tqdm, tqdm_notebook
  
  # get full file paths and file names
  file_names = [file for file in os.listdir(tpm_dir) if pattern in file]
  
  print(file_names)
  
  print(len(file_names))

  for i in range(0, len(file_names)): 
    
    # get file name
    fn = re.sub("__.*", "", file_names[i])
    print(i, file_names[i])
    
    # prepping tpm file --------------------------------------------------------------------------------
    
    TPM = pd.read_csv(open(tpm_dir+file_names[i]), encoding="utf-8", engine='python', index_col=0, header=0)
  
    if 'ENS' in TPM.index[0]:
      TPM = TPM.transpose()

    # clean preceding X from sample 
    if 'X' in TPM.index[0]:
      TPM.index = TPM.index.str.extract(r'X(.*)').iloc[:,0].tolist()
    
    # prepping metadata file --------------------------------------------------------------------------------

    # importing metadata file
    meta = pd.read_csv(open(meta_path),encoding="utf-8",engine='python',index_col=0,header=0)
    
    # subset metadata
    if meta.shape[0] != TPM.shape[0]:
      meta = meta[meta.index.isin(TPM.index)]
    
    # filter for covs to regress out
    meta = meta = meta[['pmi', 'RIN', 'library_batch', 'race', 'msex', 'study', 'age_death', 
    'age_at_visit_max', 'Unknown', 'InNeurons', 'Oligodendrocytes', 'Endothelial', 'Microglia', 
    'Astrocytes', 'OPC', 'ExNeurons']]
      
    # scaling the age cols
    scaler = StandardScaler()
    scaler.fit(meta.loc[:,['age_death', 'age_at_visit_max']])
    meta_transformed = pd.DataFrame(scaler.transform(meta.loc[:,['age_death', 'age_at_visit_max']]))
    meta_transformed.index = meta.index
    meta_transformed.columns = ['age_death', 'age_at_visit_max']
    meta_transformed.head()  
    
    # replacing unscaled cols with scaled cols
    meta.loc[:,['age_death', 'age_at_visit_max']] = meta_transformed.loc[:,['age_death', 'age_at_visit_max']]
      
    # setting correct datatypes
    meta['library_batch'] = meta['library_batch'].astype(int)
    
    print(meta.head(n=10))
  
    print(meta.shape)
    print(TPM.shape)
      
    # performing mlr across genes --------------------------------------------------------------------------------
    
    # defining empty array to hold residual tpm values
    residual_df = pd.DataFrame(np.zeros(shape=(len(TPM.index.values), len(TPM.columns.values))), 
                            index = TPM.index, 
                            columns = TPM.columns)

    shapiro_df = pd.DataFrame(np.zeros(shape=(len(TPM.columns.values), 1)), 
                            index = TPM.columns, 
                            columns = ['shapiro_pval'])
                            
    # initiating the lm
    lm = LinearRegression()

    def regress_covs_from_gene(y):

      """
      function to regress out covariates from each gene - apply over gene axis
      inserts residuals into predefined residual dataframe with correct row (sample id) and col (gene) labels
      """
      gene_name = y.name
      
      # predictors
      X = meta
      
      # response var
      y = np.array(y)
      
      # get list of locations where gene has nan values
      nan_locs = np.where(np.isnan(y))[0]
      
      # masking NAs
      finiteYmask = np.isfinite(y)
      
      # cleaning NAs from predictors and response vars
      Yclean = y[finiteYmask]
      Xclean = X[finiteYmask]
      
      # fit the model 
      lm.fit(Xclean, Yclean)
      
      # make predictions from covs - i.e. expected TPMs
      predictions = lm.predict(Xclean)
      
      # reinsert removed nan values one-by-one to make y and predictions the same length
      for i in range(0, len(nan_locs)):
        predictions = np.insert(predictions, nan_locs[i], np.nan)
      
      # residuals = observations - predictions
      residuals = y - predictions
      
      # adding residuals to residual df
      residual_df[gene_name] = residuals
      
      # # cleaning residuals of NaNs for shapiro fn. 
      # shapiro_resids = residuals[np.logical_not(np.isnan(residuals))]
      # 
      # # get shapiro pvals and add to shapiro_df
      # shapiro_pval = shapiro(shapiro_resids)[1]
      # shapiro_df.loc[gene_name, 'shapiro_pval'] = shapiro_pval
    
    # import progress bar
    from tqdm import tqdm, tqdm_notebook
    tqdm.pandas(tqdm_notebook)

    # applying regress_covs_from_gene fn. across genes
    TPM.progress_apply(regress_covs_from_gene, axis=0);
    
    # # find % genes with normally distributed residuals
    # sig_norm_gene_count = (int(shapiro_df[shapiro_df <= 0.05].count()) / len(shapiro_df.index.values)) * 100
    # print(round(sig_norm_gene_count, 2), "% of gene residuals have a significantly normal distribution")

    # export output --------------------------------------------------------------------------------

    # export filtered file
    os.chdir(out_dir)
    residual_df.to_csv(fn + "_" + outlabel + "_residuals.csv", index=True, header=True)
        
      
                            
    
    
    
    
    
    
    
    
    
    
    
