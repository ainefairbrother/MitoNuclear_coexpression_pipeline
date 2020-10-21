
def gtex_regress_covariates(tpm_dir, pattern, meta_path, pheno_path, out_dir=None, outlabel=""):

    """
    :param tpm_dir: directory containing the TPM files - genes are cols, samples are rows
    :param pattern: unique pattern to select only TPM files
    :param meta_path: path to metadata file, downloadable from GTEx i.e. GTEx_Data_V6_Annotations_SampleAttributesDS.txt
    :param pheno_path: path to phenotype file, downloadable from GTEx, gives age, sex, COD associated with patient samples
    :param out_dir: directory to save residual files into
    :param outlabel: optional label to add to output files
    :return: residuals df
    
    Corrects for the following hardcoded covs: RIN, SMNABTCHT, SMNABTCH, SMGEBTCH, SMGEBTCHD, SMCENTER, AGE, GENDER, DTHHRDY
    
    """
    
    if out_dir == None:
      out_dir = tpm_dir

    # import libs
    import pandas as pd
    import re
    import numpy as np
    from sklearn import preprocessing
    from sklearn.linear_model import LinearRegression
    from scipy.stats import shapiro
    import os
    #from tqdm import tqdm, tqdm_notebook

    # prepping metadata file --------------------------------------------------------------------------------

    # importing metadata file
    meta = pd.read_csv(open(meta_path),encoding="utf-8",sep="\t",engine='python',index_col=0,header=0)

    del meta.index.name

    # index to col
    meta.index.name = 'long_id'
    meta.reset_index(inplace=True)

    # adding short ids
    meta_ids = meta['long_id'].tolist()
    for i in range(0, len(meta_ids)):
        meta_ids[i] = re.sub(r"(GTEX-\w+)-\w+-\w+-.+", r"\g<1>", meta_ids[i])

    # adding short id
    meta.insert(0, 'short_id', meta_ids)

    # prepping pheno file --------------------------------------------------------------------------------

    # importing phenotype data file
    pheno = pd.read_csv(open(pheno_path), encoding="utf-8", sep="\t", engine='python', index_col=0, header=0)

    pheno.index.name = 'short_id'
    pheno.reset_index(inplace=True)

    # looping over TPM files --------------------------------------------------------------------------------

    # get full file paths and file names
    file_paths = [tpm_dir + file for file in os.listdir(tpm_dir) if pattern in file]
    file_names = [file for file in os.listdir(tpm_dir) if pattern in file]
    
    print(len(file_names))

    for i in range(0, len(file_paths)):

        # prepping tpm file --------------------------------------------------------------------------------

        # get file name
        fn = re.sub("_.*", "", file_names[i])
        print(i, fn)

        # importing tpm file
        TPM = pd.read_csv(open(file_paths[i]),
                          encoding="utf-8",
                          engine='python',
                          index_col=0,
                          header=0)

        del TPM.index.name
        
        TPM.index = TPM.index.str.replace(r'\.', '-', regex=True)
        
        print(TPM.head(n=5))
        
        # index to col
        TPM.index.name = 'long_id'
        TPM.reset_index(inplace=True)

        # adding short ids
        TPM_ids = TPM['long_id'].tolist()
        for j in range(0, len(TPM_ids)):
            TPM_ids[j] = re.sub(r"(GTEX-\w+)-.+", r"\g<1>", TPM_ids[j])

        # adding short id
        TPM.insert(0, 'short_id', TPM_ids)

        # cleaning covs --------------------------------------------------------------------------------

        # extracting variables from meta table
        meta = meta[["short_id", "long_id", "SMRIN", "SMNABTCHT", "SMNABTCH", "SMGEBTCH", "SMGEBTCHD", "SMCENTER"]]
        # was previously just "SMCENTER" - changed 24/04/20
        
        # joining meta and pheno by short_id
        covs = pd.merge(meta, pheno, on="short_id")

        # getting covs for ids present in tpm file
        covs = covs[covs['long_id'].isin(TPM['long_id'])]

        # set index cols for covs
        covs = covs.set_index('short_id')

        # set index cols for tpm table
        TPM = TPM.set_index('short_id')
        TPM = TPM.drop('long_id', axis='columns')

        # order tpm and cov tables the same
        TPM = TPM.reindex(covs.index)

        # encoding covariates  --------------------------------------------------------------------------------
        
        # separating the numeric and textual covs
        covs_with_labels = covs.loc[:,covs.columns.isin(covs.columns[covs.dtypes == 'object'].tolist())]
        covs_with_nums = covs.filter(covs.columns[covs.dtypes != 'object'].tolist())
        
        encoder = preprocessing.LabelEncoder()

        # fit the transformer
        for i in range(0, len(covs_with_labels.columns)):
            encoder.fit(covs_with_labels.iloc[:,i])
            covs_with_labels.iloc[:,i] = encoder.transform(covs_with_labels.iloc[:,i])
        
        encoded_covs = pd.merge(covs_with_nums, covs_with_labels, on="short_id")
        
        # fill encoded covs with mean for tissues without data for these - doesn't do any filling for brain, but does 
        # limited filling for other GTEx tissues i.e. wholeblood
        # fills with most common value for the variable in question - should make very little difference to the correction 
        encoded_covs = encoded_covs.fillna(encoded_covs.mean())
        
        # performing mlr across genes --------------------------------------------------------------------------------

        # defining empty array to hold predicted tpm values
        residual_df = pd.DataFrame(np.zeros(shape=(len(TPM.index.values), len(TPM.columns.values))),
                                   index=TPM.index,
                                   columns=TPM.columns)
        # defining empty array to hold shapiro pvals
        shapiro_df = pd.DataFrame(np.zeros(shape=(len(TPM.columns.values), 1)),
                                  index=TPM.columns,
                                  columns=['shapiro_pval'])

        # initiating the lm
        lm = LinearRegression()

        # define fn. to apply column-wise
        def regress_covs_from_gene(y):

            """
            function to regress out covariates from each gene - apply over gene axis
            inserts residuals into predefined residual dataframe with correct row (sample id) and col (gene) labels
            """
            gene_name = y.name

            # predictors
            X = encoded_covs

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

            # residuals = observations - predictions, normalising so they're normally distributed
            residuals = y - predictions

            # adding residuals to residual df
            residual_df[gene_name] = residuals
            
            # cleaning residuals of NaNs for shapiro fn. 
            shapiro_resids = residuals[np.logical_not(np.isnan(residuals))]

            # get shapiro pvals and add to shapiro_df
            shapiro_pval = shapiro(shapiro_resids)[1]
            shapiro_df.loc[gene_name, 'shapiro_pval'] = shapiro_pval

        # initiate progress bar
        #tqdm.pandas(tqdm_notebook)

        # applying regress_covs_from_gene fn. across genes
        TPM.apply(regress_covs_from_gene, axis=0);

        # find % genes with normally distributed residuals
        sig_norm_gene_count = (int(shapiro_df[shapiro_df <= 0.05].count()) / len(shapiro_df.index.values)) * 100
        print(round(sig_norm_gene_count, 2), "% of gene residuals have a significantly normal distribution")

        # export output --------------------------------------------------------------------------------

        # export filtered file
        os.chdir(out_dir)
        residual_df.to_csv(fn + "_" + outlabel + "_residuals.csv", index=True, header=True)

