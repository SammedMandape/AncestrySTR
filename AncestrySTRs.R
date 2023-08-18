library(tidyverse)
library(data.table)
library(ggplot2)
library(janitor)
#library(vtree)
library(data.table)
library(ggalluvial)

mydata <- fread("STRs.txt")

mydata %>%
  mutate(Pop = ifelse(`Sample ID`=="C1", "ASA", Pop)) -> mydata1

mydata1 %>% 
  mutate(`Final Haplotype 2` = ifelse(`Final Haplotype 2` ==  "", 
                                      `Final Haplotype 1`,
                                      `Final Haplotype 2`),
         `Final A2` = ifelse(`Final A2`=="",
                             `Final A1`,
                             `Final A2`),
  `Final Read Depth 1.1` = ifelse(is.na(`Final Read Depth 2`), `Final Read Depth 1`/2,
                                  `Final Read Depth 1`
                                  ),
  `Final Read Depth 2.1` = ifelse(is.na(`Final Read Depth 2`), `Final Read Depth 1`/2,
                                  `Final Read Depth 2`)
  )  -> foo

foo %>% group_by(Pop, `Final Haplotype 1`) %>% mutate(PopHap1Count = n()
                                ) -> foo1

foo %>% group_by(Pop, `Final Haplotype 2`) %>% mutate(PopHap2Count = n()
) -> foo2

full_join(foo1, foo2, by = c("Sample ID", "Locus", "Pop")) -> mydata2

fwrite(mydata2, col.names = T, file = "Transformed_data.csv")

mydata2 %>% filter(is.na(`Final A3.x`), !(Locus =="Amelogenin"),!(Locus =="HPRTB")) %>%
  filter(!grepl("^DX|^DY|^Y", Locus)) %>%
  pivot_longer(cols = c(`Final Haplotype 1.x`,`Final Haplotype 2.x`), 
               values_to = "Haplotype", names_to = "Haplotype_type") -> foo.temp

##################################################################################
## getting how many populations each allele belongs to, the more is found in just
## one population, it contributes to ancestry
foo.temp %>%
  group_by(Haplotype, Pop, Locus) %>%
  summarise(Count=n()) -> foo1

foo1 %>%
  group_by(Haplotype, Locus) %>% count() -> foo3.1 #foo2

# foo2 %>% ggplot(aes(y=n)) +
#   geom_bar()
# 
# foo3 %>% ggplot(aes(y=n)) +
#   geom_bar()

foo3.1 %>% ggplot(aes(y=n)) +
  geom_bar() + facet_wrap(.~Locus, scale="free")

foo3.1 %>% group_by(Locus,n) %>% summarise(Population = n()) %>% 
ggplot(aes(x=reorder(Locus, -Population),y=Population,fill=as.character(n))) +
  geom_bar(stat = "identity") +
  geom_text(stat="identity",aes(label= after_stat(y)),
            size =5,
            position = position_stack(0.5),
            color="white")
#mydata2 %>% 
#  filter(`Final Haplotype 1.x`=="CAGAGAGAAAGAATCAACAGGATCAATGGATGCATAGGTAGATGATAGATAGATAGATAGATAGATAGATAGATAGATAGACAGACAGACAGACAGACAGACAGACAGATGAGAGGGGATTTATTAGAGGAATTAGCTCAAGTGATATGGAGGCTGAAAAATCTCATGACAGTCCATCTGCAA") -> bar1
##################################################################################

##################################################################################
## get number of sequence based alleles
SEBased <- foo.temp %>% 
  ungroup() %>% 
  group_by(Locus) %>% 
  distinct(Haplotype) %>% 
  summarise(NumAllele=n()) %>%
  mutate(Type = "SE")

## get number of length based alleles
LEBased <- mydata2 %>% filter(is.na(`Final A3.x`), !(Locus =="Amelogenin"),!(Locus =="HPRTB")) %>%
  filter(!grepl("^DX|^DY|^Y", Locus)) %>%
  pivot_longer(cols=c(`Final A1.x`,`Final A2.x`), 
               values_to = "CEAlleles", names_to = "AlleleType") %>% 
  ungroup() %>%
  group_by(Locus) %>%
  distinct(CEAlleles) %>%
  summarise(NumAllele=n()) %>%
  mutate(Type = "LE")

## combine se and le based alleles
LESE <- bind_rows(SEBased, LEBased)
LESE %>% ggplot(aes(reorder(Locus, -NumAllele, sum),NumAllele, fill = forcats::fct_rev(Type))) +
  geom_bar(stat = "identity")



##################################################################################
# tranforming for infocalc 
# filter sex alleles, filter any multiallelic sites
mydata2 %>% mutate(Popcode = case_when(Pop=="CAU" ~ 1,
                                       Pop=="ASA" ~2,
                                       Pop== "HIS" ~3,
                                       Pop=="AFA"~4),
                   Country = "USA",
                   GeographicReg = "America") -> mydata3
mydata3 %>% filter(is.na(`Final A3.x`), !(Locus =="Amelogenin")) %>% select(`Sample ID`, Popcode, #`Final A1.x`, `Final A2.x`,
                                                   Pop,Country, GeographicReg,Locus,
                                                  `Final Read Depth 1.1.x`, `Final Read Depth 2.1.x`
                                                   )  %>%
  rename(Allele1=`Final Read Depth 1.1.x`,
         Allele2=`Final Read Depth 2.1.x`) %>%
  filter(!grepl("^DX|^DY|^Y", Locus))-> mydata4

str_remove(mydata4$`Sample ID`,"^0+") -> mydata4$`Sample ID`

###########################################################
# this is for STRUCTURE software to get the popcode
mydata4 %>% ungroup() %>% select(`Sample ID`, Popcode) %>% distinct() -> forjoinStr
###########################################################
# --------> but this is WRONG way as it just considers read depth rather 
# than unique haplotypes
# convert into long format for all alleles in a locus and
# then into a wide format
# set any missing values to -9 as in doc of infocalc
# mydata4 %>% 
#   pivot_longer(cols=c(Allele1, Allele2),
#                names_to = "Genotype", 
#                values_to = "AlleleRD") %>%
# pivot_wider(names_from = Locus, 
#             names_sep = ".",
#             values_from = (AlleleRD)) %>% select(-Genotype) -> mydata5
# 
# mydata5[is.na(mydata5)] <- -9
# #mydata5 %>% is.na() -> foo
# fwrite(mydata5,"infocalc_input.txt", sep = " ")
##################################################################################

##################################################################################
## --------> but this is the RIGHT way
## SE based infocalc analysis
# Changing names for pentaD and pentaE
mydatstruc_1 %>% mutate(Locus = ifelse(Locus=="PENTAD","PentaD",Locus),
                        Locus = ifelse(Locus== "PENTAE","PentaE",Locus)
                        ) -> mydatstruc_1.1

# 2 diff ways to left-join as there is some diff in the # of alleles
# some alleles are more in SEbased repeat file..further analysis 
# showed about 5 loci were more in SEbased file. So decided
# to go with first approach of left-join 
mydata6 <- mydata4 %>% 
  left_join(mydatstruc_1.1, by=c("Sample ID"="Sample ID","Locus"="Locus"))

# mydata6 <- mydatstruc_1.1 %>% 
#   left_join(mydata4, by=c("Sample ID"="Sample ID","Locus"="Locus")) %>%
#   #filter(is.na(Popcode))
#   filter(!(Locus=="PentaE-long")) %>% filter(is.na(Popcode))


mydata6 %>% select(-c(Allele1.x, Allele2.x)) %>%
  filter(Locus!="HPRTB") %>%
  pivot_longer(cols = c("Allele1.y","Allele2.y"),
               names_to = "Genotype",
               values_to = "SERepeat") %>%
  pivot_wider(names_from = Locus,
              values_from = SERepeat) %>% select(-Genotype) -> mydata7

mydata7[is.na(mydata7)] <- -9
fwrite(mydata7, "infocalc_input_SERepeat.txt", sep=" ")
###############################################
###############################################
## CE based infocalc analysis

mydatstrucCE_1 %>% mutate(Locus = ifelse(Locus=="PENTAD","PentaD",Locus),
                          Locus = ifelse(Locus== "PENTAE","PentaE",Locus)
) -> mydatstrucCE_1.1

mydata6CE <- mydata4 %>% 
  left_join(mydatstrucCE_1.1, by=c("Sample ID"="Sample ID","Locus"="Locus"))

mydata6CE %>% select(-c(Allele1.x, Allele2.x)) %>%
  filter(Locus!="HPRTB") %>%
  pivot_longer(cols = c("Allele1.y","Allele2.y"),
               names_to = "Genotype",
               values_to = "CERepeat") %>%
  pivot_wider(names_from = Locus,
              values_from = CERepeat) %>% select(-Genotype) -> mydata7CE

mydata7CE[is.na(mydata7)] <- -9
fwrite(mydata7CE, "infocalc_input_CERepeat.txt", sep=" ")
###############################################
###############################################
## Plotting rosenbergs informativeness

myinfo <- fread("Infocalc_Rosenbergs_informativeness_CE-SE.txt", sep = "\t",
                header = T)
mypopstr <- fread("PopSTR_loci.txt", sep = "\t", header = T) %>%
  mutate(Locus=Loci, I_n = In, Type="popSTR",) %>% select(-In, -Loci)
myinfo %>% 
  mutate(Type=ifelse(Type=="CE","Length based","Sequence based"))->myinfo1

myinfo_popstr <- bind_rows(mypopstr, myinfo1)
myinfo_popstr %>%
  ggplot(aes(reorder(Locus, -I_n, sum),I_n, group=Type, color=Type)) +
  geom_line(size=1.2, linetype="dashed") +
  geom_point( size = 3) +
  labs(#title = "Allelic Diveristy",
    y = expression(Rosenbergs~informativeness~(I[n])),
    x= "Loci",
    color= "Allele type")+
  theme_bw()+
  scale_color_viridis_d(begin = 0.25, end = 0.75, direction = -1) +
  theme(legend.position = c(0.7, 0.8), 
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(angle=90, vjust = 0.5, hjust = 1),
        text = element_text(size = 14)) #+
  #scale_colour_manual(values = c("#f98e09","#0d0887"))

### separating LE vs popstr and SE vs popstr
myinfo_popstr %>% filter(Type=="Sequence based") %>%
  mutate(Type1="Length based vs popSTR") %>% 
  ggplot(aes(reorder(Locus, -I_n, sum),I_n, group=Type, color=Type)) +
  geom_line(size=1.2, linetype="dashed") +
  geom_point( size = 3) +
  labs(#title = "Allelic Diveristy",
    y = expression(Rosenbergs~informativeness~(I[n])),
    x= "Loci",
    color= "Allele type")+
  theme_bw()+
  scale_color_viridis_d(begin = 0.25, end = 0.75, direction = -1) +
  theme(legend.position = c(0.7, 0.8), 
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(angle=90, vjust = 0.5, hjust = 1),
        text = element_text(size = 14))
