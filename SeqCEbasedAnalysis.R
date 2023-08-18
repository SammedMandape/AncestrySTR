library(tidyverse)
library(data.table)
library(ggplot2)
library(janitor)
#library(vtree)
library(data.table)
library(ggalluvial)

## data transformation for STRUCTURE tool

######################################################################################
## SEQUENCE based STR and SNP analysis
######################################################################################
mydatstruc <- fread("AncestrySTRs_genotypescodedasnumbers_NNedited071323.txt", 
                    sep = "\t",
                    colClasses = 'character', data.table = FALSE)

mydatstruc %>%
  pivot_longer(cols = -`Sample ID`, names_to = "Locus", values_to = "Alleles") %>%
  filter(!(Locus =="Amelogenin"), !(Locus =="HPRTB")) %>%
         #!(`Sample ID`=="10561"), !(`Sample ID`=="17467"),!(`Sample ID`=="47812"),
         #!(`Sample ID` %in% c("11960","C1","19971","78755","C62"))) %>%
  filter(!grepl("^DX|^DY|^Y", Locus)) %>% #-> mydatstruc_1 #%>%
  separate(Alleles, into = c("Allele1","Allele2"), sep = "/") %>%
  mutate(Allele1 = as.numeric(Allele1), Allele2 = as.numeric(Allele2))-> mydatstruc_1 #%>%
  #filter(!is.na(Allele2))

#mydatstruc_1 %>% filter(!is.na(Allele2))
mydatstruc_1[is.na(mydatstruc_1)] <- -9

mydatstruc_1 %>% 
  pivot_longer(cols = c(Allele1, Allele2),
               names_to = "Genotype",
               values_to = "Repeat") %>%
  mutate(Repeat = as.numeric(Repeat)) %>%
  pivot_wider(names_from = Locus,
              names_sep = ".",
              values_from = (Repeat)) %>% select(-Genotype) %>%
  rename("SampleID"="Sample ID")-> mydatstruc_2

mydatstruc_2 %>% 
  left_join(forjoinStr,
            by=c("SampleID"="Sample ID")) %>%
  select(SampleID, Popcode, everything())-> mydatstruc_2.1

mydatstruc_2.1 %>% distinct(SampleID)
fwrite(mydatstruc_2.1, "Structure_input.txt", sep = " ", col.names = F)

mydatstruc_2.1 %>% group_by(Popcode) %>% summarise(n())
#######################################################################################

#######################################################################################
## CE_STR based analysis
#######################################################################################
mydatstrucCE <- fread("AncestrySTRs_genotypescodedasnumbers_CE-STR.txt", 
                    sep = "\t",
                    colClasses = 'character', data.table = FALSE)

mydatstrucCE %>%
  pivot_longer(cols = -`Sample ID`, names_to = "Locus", values_to = "Alleles") %>%
  filter(!(Locus =="Amelogenin"), !(Locus =="HPRTB")) %>%
  #!(`Sample ID`=="10561"), !(`Sample ID`=="17467"),!(`Sample ID`=="47812"),
  #!(`Sample ID` %in% c("11960","C1","19971","78755","C62"))) %>%
  filter(!grepl("^DX|^DY|^Y", Locus)) %>% #-> mydatstruc_1 #%>%
  separate(Alleles, into = c("Allele1","Allele2"), sep = "/") %>%
  mutate(Allele1 = as.numeric(Allele1), Allele2 = as.numeric(Allele2))-> mydatstrucCE_1 #%>%
#filter(!is.na(Allele2))

#mydatstruc_1 %>% filter(!is.na(Allele2))
mydatstrucCE_1[is.na(mydatstrucCE_1)] <- -9

mydatstrucCE_1 %>% 
  pivot_longer(cols = c(Allele1, Allele2),
               names_to = "Genotype",
               values_to = "Repeat") %>%
  mutate(Repeat = as.numeric(Repeat)) %>%
  pivot_wider(names_from = Locus,
              names_sep = ".",
              values_from = (Repeat)) %>% select(-Genotype) %>%
  rename("SampleID"="Sample ID")-> mydatstrucCE_2

mydatstrucCE_2 %>% 
  left_join(forjoinStr,
            by=c("SampleID"="Sample ID")) %>%
  select(SampleID, Popcode, everything())-> mydatstrucCE_2.1

mydatstrucCE_2.1 %>% distinct(SampleID)
fwrite(mydatstrucCE_2.1, "StructureCE_input.txt", sep = " ", col.names = F)

mydatstrucCE_2.1 %>% group_by(Popcode) %>% summarise(n())

#####################################################################################


