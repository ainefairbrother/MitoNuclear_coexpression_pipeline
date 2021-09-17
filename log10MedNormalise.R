#! /usr/bin/Rscript

# function to log10 median normalise counts
# log10 transformation to make sample distributions normal prior to median normalisation 
# median normalisation to make the sample expression medians the same for inter-sample comparability
# runs all files in parallel
# takes a dataframe in the form samples = cols, genes = rows

log10MedNormalise = function(file_dir="", pattern=".csv", filesep=",", outlabel="", med_norm=TRUE){
  
  # file_dir: directory where input file is stored
  # pattern: string pattern to recognise multiple files - if using one file, simply give full file name
  # filesep: separator of file (i.e. \t or ,)
  # outlabel: label to add to the output file. Allows next pipeline fn. to find it easily
  # med_norm: can disable median normalisation by setting this to FALSE
  
  library(beadarray)
  library(parallel)
  
  # import files and store in list
  file_list = list.files(path=file_dir, pattern=pattern)
  
  # run all files in parallel - use n cores = n files
  num_cores = length(file_list)
  
  # defining normalisation fn 
  doNormalisation = function(f){
    
    file_name = gsub("__.*", "", f)
    
    print(file_name)
    
    df = read.csv(file=paste0(file_dir, f), sep=filesep, row.names=1, header=T)
    
    # samples to cols
    if(length(colnames(df)) > length(rownames(df))){
      df = data.frame(t(df))
    }
    
    print("loading df")
    print(dim(df))
    print(df[1:4,1:4])
    
    print("getting all NA in df")
    print(sum(is.na(df)))
    
    # remove individuals that have zero reads
    print("removing individuals with 0 reads")
    test = apply(df, 2, function(x) { any(x > 0)})
    new = df[,test]
    df = new
    print(dim(df))
    
    # remove genes that have zero reads
    print("removing genes with 0 reads")
    test = apply(df, 1, function(x) { all(x > 0)})
    new = df[test,]
    df = new
    print(dim(df))
    
    print("Removing genes that are all NA")
    df = df[rowSums(is.na(df)) != ncol(df), ]
    print(dim(df))
    
    print("getting all NA in filtered df")
    print(sum(is.na(df)))
    
    print("log10 normalising")
    pseudoCount = log10(df + 1)
    
    print("pseudoCount table: ")
    
    print("writing outfiles")
    if(med_norm==TRUE){
      # median normalise with the "beadarray" module in R
      # medianNormalise normalises expression across columns - so samples should be cols
      MedianpseudoCount = medianNormalise(pseudoCount, log=F)
      print(dim(MedianpseudoCount))
      print("exporting file")
      write.csv(MedianpseudoCount, file=paste0(file_dir, gsub(pattern, "", file_name), "_", outlabel, "_log10_mediannorm_TPM.csv"))
    }
    
    if(med_norm==FALSE){
      write.csv(pseudoCount, file=paste0(file_dir, gsub(pattern, "", file_name), "_", outlabel, "_log10_norm_TPM.csv"))
    }
  }
  
  # parallelise normalisation over all files
  mclapply(file_list, doNormalisation, mc.cores=num_cores)
  #mclapply(file_list, doNormalisation, mc.cores=num_cores+1)
}



