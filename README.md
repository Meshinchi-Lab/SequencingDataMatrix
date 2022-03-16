# Analysis Directory 
### for SequencingDataMatrix

Objective: create a sample manifest for RNA-seq, miRNA-seq, WGS, DNA methylation datasets. These will be used associated clinical data elements with the sample IDs in the sequencing data. These can be used to generate multi-assay experiments or summarized experiment objects, dge objects, etc. 

1) Append cleaned sample manifests from BCCA to the current manifests found on the Fred Hutch Meshinchi Fast drive. The code in `TARGET_AML_Sequencing_Manifests_05.05.21.Rmd` can be used a template. 

The sample manifests from BCCA are provided typically with on the FTP site and often have filenames such as. "library.summary", and are simple tab delimeted text files mapping the provided sample ID with the BCCA sample IDs. 


2) Create a symlink to copy the sample manifest for the RNA-seq datasets to the working copies of the RNA-seq countr matrices. 

```
cp -s $PWD/TARGET_AML_Ribodepleted_Manifest_08.12.21.csv /fh/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/0000.00.03_ExpressionMatrices/TARGET_AML_Ribodepleted_Manifest_08.12.21.csv
```


Author: Jenny Leopoldina Smith<br>
ORCID: [0000-0003-0402-2779](https://orcid.org/0000-0003-0402-2779)
<br>
