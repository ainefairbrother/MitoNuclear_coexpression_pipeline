# takes in a matrix of TPMs/RPKMs, rows=samples , cols=genes
# applies gene rpkm > 0 filter to each region i.e. gene must be expressed in all samples
# outputs filtered files
# AFB 16/01/2020

def filterNullGenesAndSamps(file_dir="", pattern=".csv", filesep=",", outdir=None, outlabel=""):
  
  # """
  # file_dir: directory containing the rpkms
  # pattern: pattern to search for, default is all .csv files
  # outdir: where should the files go, default is file_dir
  # filesep: file separating char, default is comma
  # outlabel: add label to output files
  # """
  
  import os 
  import pandas as pd
  
  if outdir == None:
    outdir = file_dir
  
  # get all file paths
  file_paths = [file_dir+file for file in os.listdir(file_dir) if pattern in file]
  file_names = [file for file in os.listdir(file_dir) if pattern in file]
  print(file_names)
  
  for i in range(0,len(file_paths)):
    
    print(file_names[i])
    
    # import df
    df = pd.read_csv(open(file_paths[i]), sep=filesep, encoding="utf-8", engine='python', index_col=0, header=0) # genes cols, samples rows
    
    # genes to cols 
    if len(df.index.values) > len(df.columns.values):
      df = df.transpose()
    
    # removing samples with all vals = 0
    df = df.loc[~(df==0).all(axis=1)]

    # removing genes with all 0s
    df = df.loc[:, ~(df==0).all(axis=0)] 
    
    # retaining cols where all vals are not equal to 0
    df = df.loc[:, (df!=0).all(axis=0)]
    
    # export filtered file
    os.chdir(outdir)
    df.to_csv(file_names[i].replace('.csv', '') + "_" + outlabel + "_0filtered.csv", index=True, header=True)
  
# def filterGenesInAllFiles(file_dir="", pattern=".csv", filesep=",", outdir=None, outlabel="", med_fill_na=False):
#   
#   # """
#   # file_dir: directory containing the rpkms
#   # pattern: pattern to search for, default is all .csv files
#   # outdir: where should the files go, default is file_dir
#   # filesep: file separating char, default is comma
#   # outlabel: add label to output files, default is none
#   # med_fill_na: fill NAs in genes with the gene median val, default is False
#   # """
#   
#   import os 
#   import pandas as pd
#   
#   if outdir == None:
#     outdir = file_dir
#   
#   # get all file paths
#   file_paths = [file_dir+file for file in os.listdir(file_dir) if pattern in file]
#   file_names = [file for file in os.listdir(file_dir) if pattern in file]
#   print(file_names)
#   
#   print("getting genes in {} files".format(len(file_names)))
# 
#   master_genelist = []
#   for i in range(0,len(file_paths)):
#     
#     print(i, file_names[i])
#     
#     # import df
#     df = pd.read_csv(open(file_paths[i]), sep=filesep, encoding="utf-8", engine='python', index_col=0, header=0, nrows=1) # genes cols, samples rows
#     
#     # genes to cols 
#     if 'ENS' not in df.columns.values[1]:
#       print('warning: DF not in correct orientation, so loading whole file')
#       df = pd.read_csv(open(file_paths[i]), sep=filesep, encoding="utf-8", engine='python', index_col=0, header=0) # genes cols, samples rows
#       df = pd.DataFrame(df.transpose())
# 
#     # adding to master list
#     genelist = df.columns.values.tolist()
#     master_genelist.append(genelist)
#   
#   # get genes in common to all region files and save
#   list_intersect = set.intersection(*[set(list_) for list_ in master_genelist])
#   
#   print("filtering files for gene set of length {}".format(len(list_intersect)))
# 
#   # loop to filter and output filtered files
#   for i in range(0,len(file_paths)):
# 
#     # import df
#     df = pd.read_csv(open(file_paths[i]), sep=',', encoding="utf-8", engine='python', index_col=0, header=0) # genes cols, samples rows
#     
#     if(len(df.index.values) > len(df.columns.values)):
#       df = pd.DataFrame(df.transpose())
#     
#     if 'brain_region' in df.columns.values:
#       df.drop('brain_region', axis='columns', inplace=True)
#     
#     # median filling df
#     if med_fill_na == True:
#       print("na count: ", df.isnull().values.sum())
#       df.fillna(df.median(), inplace=True)
#       print("after median impute: ", df.isnull().values.sum())
# 
#     df = df.loc[:,list_intersect]
#     print(i, file_names[i], df.shape)
# 
#     # export filtered file
#     os.chdir(outdir)
#     df.to_csv(file_names[i].replace('.csv', '') + "_" + outlabel + "_genefiltered.csv", index=True, header=True)




