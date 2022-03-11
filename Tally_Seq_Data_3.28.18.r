#Jenny Smith
#11/27/17

#combine the tallied data together. 


library(dplyr)
library(magrittr)
library(tidyr)
library(stringr)



setwd("X:/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/analysis/2017.08.18_RNAseq_TallyperPatient/")



#Function to rename some columns with a suffix. 
set_colNames <- function(df){
  suffix <- substitute(df) #use the name of the dataframe as suffix
  
  idx <- colnames(df) != "USI"
  colnames(df)[idx] <- paste(colnames(df)[idx], suffix, sep=".")
  
  return(df)
}


#CDE

merged <- read.csv("H:/reference_mapping-files/TARGET_AML_1031_0531_Merged_CDE_3.22.18.csv",
                   stringsAsFactors = FALSE)

merged.subset <- merged %>%
  filter(!is.na(TARGET.USI.1))

head(merged.subset)


NBM <- read.csv("H:/reference_mapping-files/Inventory_of_normal_bone_marrows.csv", 
                stringsAsFactors = FALSE)

NBM <- NBM %>%
  mutate(Protocol=rep("NBM", nrow(.)),
         SEX=ifelse(grepl("F", SEX), "Female", "Male")) %>%
  select(TARGET.USI.1=USI,Age.Yrs=AGE,Gender=SEX, everything())

head(NBM)
dim(NBM) #72



#Merge the 0531 and 1031 CDEs 

merged.nbm <- merged.subset %>%
  bind_rows(.,NBM[,1:3]) %>%
  select(USI=TARGET.USI.1, everything())


dim(merged.nbm) #2178 AML samples + 72 NBM == 2250 Samples
head(merged.nbm)


# write.csv(merged.nbm, "TARGET_AML_1031_0531_NBM_Merged_CDE_3.28.18.csv", row.names = FALSE)




#Tallied Data matrices

#WGS
wgs <- read.csv("byDataType/TARGET_AML_0531_1031_WGS_DataAvailability.csv", stringsAsFactors = FALSE) 
wgs <- set_colNames(wgs)

head(wgs)
dim(wgs)



#TCS
tcs <- read.csv("byDataType/TARGET_AML_0531_1031_TCS_DataAvailability.csv", stringsAsFactors = FALSE)


tcs <- tcs %>%
  select(1:7,Indels_Somatic=Somatic_Indels, SNVs_Somatic=Somatic_SNVs)

tcs <- set_colNames(tcs)


head(tcs)
dim(tcs)


#DNA methylation  
DNAme <- read.csv("byDataType/TARGET_AML_0531_1031_Methylation_Updated_DataAvailability.csv", stringsAsFactors = FALSE) 

DNAme <- DNAme %>%
  mutate(CopyNumber=ifelse(InfiniumMethylation_EPIC == 1, 1, 0)) %>%
  select(1:8,13,9:12)


DNAme <- set_colNames(DNAme) 

head(DNAme)
dim(DNAme) #2250 by 13

#NOTE: Had to manually update the Methylation data in "split_27k_methylation_matrix.r". 


#mRNA
mRNAseq.tally <- read.csv("byDataType/TARGET_AML_0531_1031_mRNAseq_DataAvailability.csv", stringsAsFactors = FALSE) 


#Update 1031 with more detailed RNAseq information
mRNAseq.1031 <- read.csv("byDataType/HD_1031_RNAseq_Samples.csv", stringsAsFactors = FALSE) #colnames from TPM matrix

HD.0531 <- dir("X:/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/gene/2014Aug_BCCA_0531_Illumina_data/",
               pattern = "*.gene.*") %>%
  str_split_fixed(.,pattern = "\\-",n=5) %>%
  subset(., select=3) %>%
  as.character() %>% 
  unique()

LD.0531 <- dir("X:/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/RNA/mRNAseq/level3/gene/2016Apr_BCCA_0531_LowDepth_Illumina_data/", 
               pattern = "*transcript*") %>%
  str_split_fixed(.,pattern = "\\-",n=5) %>%
  subset(., select=3) %>%
  as.character() %>% 
  unique()


#Need fusion data to be up to date
fus1 <- read.csv("../2017.09.15_Fusions_K.Mungall/Fusions_with_subsets_0531.csv",
                 stringsAsFactors = FALSE)

fus2 <- read.csv("../2017.09.15_Fusions_K.Mungall/Merged_Results/No_Filtering/TARGET_AML_1031_Fusions_batch123_Formated_11.3.17.csv",
                 stringsAsFactors = FALSE)

length(unique(fus1$sampleID)) + length(unique(fus2$USI)) #approximately 1,277


#Combine the more detailed 0531 and 1031 annotations to 
mRNAseq <- mRNAseq.tally %>%
  mutate(fusion=ifelse(USI %in% fus1$sampleID | USI %in% fus2$USI, 1,fusion),
         LowDepth_RNASeq=ifelse(USI %in% LD.0531, 1, 0),
         HighDepth_RNASeq=ifelse(USI %in% HD.0531, 1, 0)) %>%
  mutate(HighDepth_RNASeq=ifelse(USI %in% mRNAseq.1031$USI, 1,HighDepth_RNASeq))


mRNAseq <- set_colNames(mRNAseq) 
head(mRNAseq)
dim(mRNAseq)

sum(mRNAseq$fusion.mRNAseq) #1,308, the majority of which are from TARGET-21 substudy 
sum(mRNAseq$LowDepth_RNASeq.mRNAseq) #466, bc 446 AML samples had CDE, and 20 NBMs 
sum(mRNAseq$HighDepth_RNASeq.mRNAseq) #1,279 due to 1 sample in 0531 & 1031, and 3 samples ineligable/no CDE. 

#Find which USIs not in CDE
# "PANISJ" "PANMTU" "PATLDZ"
#PATLDZ was removed  due to down-sydrome constitutional trisomy 21 (I beleive)
#cannot find CDE in 0531 or 1031, because they are 03P1. must request CDE if available. 



#miRNAseq
miRNAseq <- read.csv("byDataType/TARGET_AML_0531_1031_miRNAseq_DataAvailability.csv", stringsAsFactors = FALSE) 
miRNAseq <- set_colNames(miRNAseq)

head(miRNAseq)
dim(miRNAseq)


#Combine the individual data frames

USI.Proctols <- merged.nbm[,c("USI", "Protocol")]
dim(USI.Proctols)


SeqData <- Reduce(function(x,y) inner_join(x,y, by="USI"), list(wgs,tcs,DNAme,mRNAseq,miRNAseq)) %>%
  inner_join(., USI.Proctols, by="USI") 
  
head(SeqData)
dim(SeqData) #2250   and 66 Columns


#Clean up the combined data
SeqData.tidy <- SeqData %>%
  select(-which(grepl("X\\.|\\.1\\.", colnames(.)))) %>%
  mutate(TARGET.Project=ifelse(USI %in% AML0531$USI, "Yes", "No")) %>%
  select(USI,Protocol,TARGET.Project,
         WGS.DataAvailable=WGS.wgs,
         TCS.DataAvailable=TCS.tcs,
         Methylation.DataAvailable=Methylation.DNAme,
         mRNAseq.DataAvailable=mRNAseq.mRNAseq,
         miRNAseq.DataAvailable=miRNAseq.miRNAseq,
         everything(),
         -matureMiRNA.miRNAseq, #matureMiRNA.miRNAseq will need to be updated by opening the cated mature miRNA matrices
         -raw_idat.DNAme) %>% #raw_idat files are not currently renamed with TARGET USIs
  set_rownames(.$USI)



head(SeqData.tidy)
dim(SeqData.tidy) #2250 by 61 cols 


table(SeqData.tidy$TARGET.Project)

sapply(SeqData.tidy[,4:8], table)

idx <- grep("Rela|Refr|TARGET21",  colnames(SeqData.tidy))
sapply(SeqData.tidy[,idx], table)

# write.csv(SeqData.tidy, "TARGET_AML_0531_1031_SeqDataAvailability_3.29.17.csv", row.names = FALSE)




#Merge the Clinical Data and the Tallied counts

withClinical <- inner_join(merged.nbm, SeqData.tidy, by = "USI")

head(withClinical)
dim(withClinical) #2250 samples by 96 cols

# write.csv(withClinical,"TARGET_AML_0531_1031_SeqDataAvailability_CDEs_3.29.18.csv")



######## Compare with the SRA table ##########################

SeqData.tidy <- read.csv("X:/fast/meshinchi_s/workingDir/TARGET/AML_TARGET/SequencingDataMatrix/TARGET_AML_0531_1031_SeqDataAvailability_3.29.17.csv",
                         stringsAsFactors = FALSE)


sratab <- read.table("TARGET_AML_SRA_Run_Table_4.2.18.txt",
                     stringsAsFactors = FALSE, sep="\t", header = TRUE)


head(sratab)
dim(sratab) #4,859 by 41

barcode <- str_split_fixed(sratab$Sample_Name, "-", n=5)

sratab.1 <- sratab %>%
  mutate(USI=barcode[,3],
         DiagnosticSample=ifelse(grepl("09A|03A", barcode[,4]), "Yes","No"),
         RemissionSample=ifelse(grepl("10A|14A", barcode[,4]), "Yes","No"),
         RelapseSample=ifelse(grepl("04A|40A", barcode[,4]), "Yes", "No"),
         RefractorySample=ifelse(grepl("41A|42A", barcode[,4]), "Yes", "No"),
         TARGET21=ifelse(grepl("21", barcode[,2]), "Yes", "No")) %>%
  select(-which(grepl("X", colnames(.)))) %>%
  arrange(USI) %>%
  select(USI,DiagnosticSample,RemissionSample, RefractorySample, RelapseSample,TARGET21, everything()) 


# head(sratab.1[,1:5], n=25)
length(unique(sratab.1$USI)) #2,071



onlyOurs <- setdiff(SeqData.tidy$USI, sratab.1$USI) #220 unique samples present only in our data
onlySRA <- setdiff(sratab.1$USI, SeqData.tidy$USI) #41 


SRA.miRNAseq <- sratab.1 %>%
  filter(Assay_Type == "miRNA-Seq") %>%
  select(1:7,"Sample_Name", "Library_Name","LoadDate", "SRA_Sample")

dim(SRA.miRNAseq) #1900 samples by  11 cols


miRNAseq.withLibs <- SeqData.tidy %>%
  filter(miRNAseq.DataAvailable == 1) %>%
  select(USI,Protocol,miRNAseq.DataAvailable) %>%
  full_join(., SRA.miRNAseq, by="USI")



dim(miRNAseq.withLibs) #1908 by  13
# write.csv(miRNAseq.withLibs, "TARGET_AML_miRNAseq_Samples.csv", row.names = FALSE)


NoProtocol <- miRNAseq.withLibs %>% 
  filter(is.na(Protocol)) %>%
  filter(! grepl("^BM|^RO|^Kas|^MV",USI)) %>%
  select(USI) %>%
  unlist()

NoProtocol

sum(NoProtocol %in% merged$TARGET.USI.1) #NO CDE 






