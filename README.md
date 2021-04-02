
Decon-coreg is a modification of Decon-eQTL to allow interaction with gene expression instead of with SNPs



### For compiling
commons-cli v1.3.1  
commons-csv v1.2  
commons-io v2.4  
commons-lang3 v3.4  
commons-math3 v3.6  
jsci-core  
jsci-mathmlimpl  
jsci-sci

## minimal usage example
    
    java -jar Decon-eQTL-v*.*.*-jar-with-dependencies.jar -c <file containing cellcounts> \
                            -e <file containing expression data> \
                            -o <output directory> \
                            -sn <file with gene gene combinations to test>
    

## Options overview

    -c,--cellcount <file>                     Cellcount file name
    -e,--expression <file>                    Expression file name
    -gc,--genotypeConfigurationType <arg>     Which genotype configuration type to use (either all or two)
    -help                                     print this message
    -no,--no_console                          Do not output logging info to the console
    -o,--outfolder <path>                     Path to folder to write output to
    -oe,--outputPredictedExpression           Write output file with predicted expression
    -of,--outfile <file>                      Outfile name of deconvolution results (will be written in outfolder)
    -r,--round_dosage                         Round the dosage to the closest int
    -sn,--snpsToTest <file>                   Tab delimited file with first column gene name, second column SNP name. Need to match with names from genotype and expression files.
    -t,--test_run                             Only run deconvolution for 100 QTLs for quick test run
    -w,--whole_blood_qtl                      Add whole blood eQTL (pearson correlation genotypes and expression)
