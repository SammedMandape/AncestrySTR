#Creator: Mary Anne Panoyan 
# Date of last edit: Nov 18 2022
#STR ANCESTRY PROJ

#goal of the R script to produce counts of the alleles to determine how often they appear in 
#varying populations 
#GO TO LINE 135 FOR NEW CODE 

#HRPTB is an x chromosomal linked loci make sure to remove it from the data!!

#Load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(janitor)
library(vtree)
library(data.table)
library(ggalluvial)
#remotes::install_github("davidsjoberg/ggsankey")

#check if we are in the current working directory 
getwd()

#load the dataset 
#df is the master file from excel 
df <-  read.csv(
  "/Users/maryannepanoyan/Desktop/projects/current projects/STR_Ancestry_Project/RE__[ZIP_ATTACHED][EXT]_STR_Ancestry_Project/STR_DATA_MASTER_FILE copy.csv")
#remove hprtb
df<-subset(df, Locus!="HPRTB")
#df3 is the transformed file from Sammed's code 
df3 <-  read.csv(
  "/Users/maryannepanoyan/Desktop/projects/current projects/STR_Ancestry_Project/RE__[ZIP_ATTACHED][EXT]_STR_Ancestry_Project/Transformed_data_MP.csv")
df3<-subset(df3, Locus!="HPRTB")

############## 

#using df3 which contains pophapcount from Sammed's code creating facet wrapped 
#bar plots displaying the  counts of haplotypes per locus by population
ggplot(df3, aes(fill=population, y=PopHap1Count, x=Locus)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~Locus, scales="free_x", "free_y") +
  labs( title = "Number of haplotype 1 for each population by locus",
        y = "Number of alles",
        x = "Locus ") +
  theme_bw()+
  scale_fill_hue()

ggplot(df3, aes(fill=population, y=PopHap2Count, x=Locus)) + 
  geom_bar(position="dodge", stat="identity") +
  facet_wrap(~Locus, scales="free_x", "free_y") +
  labs( title = "Number of haplotype 2 for each population by locus",
        y = "Number of  alles",
        x = "Locus ") +
  theme_bw()+
  scale_fill_hue()


#creating a duplicate of the original STR master file to manipulate 
df2 <- data.frame(df)


#creating the df as a tibble for dyplr
df2 <- as_tibble(df2)

#determining the number of unique counts per population including both haplotypes across all loci
df2 %>%
  group_by(population) %>%
  summarize(unique_counts = n_distinct(unlist(across(Final.Haplotype.1: Final.Haplotype.2))))

# when checking how many alleles in  haplotype 1 and 2 across all loci per population 
# 1 AFA                  526
# 2 ASA                  416
# 3 CAU                  431
# 4 HIS                  438

#creating a df to plot 
pop <- c("AFA", 'ASA', 'CAU', 'HIS')
unique_count <- c(526, 416, 431, 438)
count_df <- data.frame(pop, unique_count)
#creating a bar plot for total number of unique plots 
ggplot(data = count_df, aes(y = unique_count, x = pop, fill = pop)) +
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs( title = "Number of alleles for each population across all loci",
        y = "Number of alles",
        x = "Poplation ") +
  scale_fill_hue()

#number of unique haplotypes including haplotype 1 and 2 
summary_df <- df2 %>%
  count(population, Locus, Final.Haplotype.1, Final.Haplotype.2, sort = TRUE)

ggplot(data = summary_df, aes(y = n, x = Locus, fill = population)) +
  geom_bar(stat = "identity", position = "dodge") + 
  theme_bw() + 
  labs( title = "Number of  alleles (including both haplotypes) for each population by locus",
        y = "Number of  alles",
        x = "Locus ") +
  scale_fill_hue() +
  facet_wrap(~Locus, scales="free_x", "free_y") +
  theme_bw()

#determining the number of unique hap1 by how often it appears in each population by locus
#using the tabyl function from the janitor package 
table1<-tabyl(df2,Locus, population, Final.Haplotype.1)
head(table1)

#same as above using hap
table2<-tabyl(df2,Locus, population, Final.Haplotype.2)

#creating a table showing how often the unique haplotype appears in each population
#does not differentiate by locus
table3<-tabyl(df2, population, Final.Haplotype.1)
fwrite(table3, col.names = T, file = "table_of_haplotype1_by_pop.csv")

table4<-tabyl(df2, population, Final.Haplotype.2)
fwrite(table4, col.names = T, file = "table_of_haplotype2_by_pop.csv")



#potential plot 
#using the vtree package to display a variable tree
#the vtree will display each unique haplotype present for each population by locus
#for the example below we start with the asa population looking just at the d10 locus
#we can then see that there are 8 alleles present in this population for that specific locus
#there is one allele that is very present in the population contributing to 50% of the alleles 
#with 3 alleles being represented under 1%
#the main issue with this diagram is that it would produce 4x28 figures, which would be time 
#consuming 

#vtree example 
vtree(df2, c("population","Locus","Final.Haplotype.1"),
      fillcolor = c(population = "#e7d4e8",
                                  Locus = "#99d8c9",
                                  Final.Haplotype.1 = "#9ecae1"),
      keep = list(population = "ASA", Locus ="D10S1248"))

####################################
#new section, created after Nov 18 meeting
#Last updated: Nov 28

#pivot the data frame into a long format

df1 <- df[,c('Final.Haplotype.1', 'Final.Haplotype.2', "Final.A1", "Final.A2", "Locus","population")]
df1<-df1 %>% pivot_longer(cols=c('Final.Haplotype.1', 'Final.Haplotype.2'),
                    names_to='haplotype',
                    values_to='sequence_') 
df1<-df1 %>% pivot_longer(cols=c('Final.A1', 'Final.A2'),
                          names_to='allele',
                          values_to='length_') 

#group by loci and population shows us the number of unique haplotypes in each pop
# df_grp1 = df1 %>% group_by(Locus, population) %>%
#   summarise(counts = n_distinct(sequence_))

#number of unique haplotypes by locus
# ggplot(df_grp1, aes(y=counts, x=Locus, fill= population)) + 
#   facet_wrap(~Locus, scales="free_x") +
#   geom_bar(position="dodge", stat="identity") +
#   labs( title = "Number of unique haplotypes each population by locus",
#         y = "Number of unique haplotypes",
#         x = "Locus ") +
#   theme_bw()+
#   scale_fill_hue()

#USE THIS GRAPH
#number of populations for haplotypes by locus
#group by loci and sequence to get how many times each sequence occurs in pop
df_grp2 = df1 %>% group_by(Locus, sequence_, length_) %>%
  summarise(pops = n_distinct(population))

df_grp2<- df_grp21[order(-df_grp21$pops),]

ggplot(df_grp2, aes(y=pops, x=reorder(sequence_,-pops), fill=factor(length_))) + 
  facet_wrap(~Locus, scales="free_x",) +
  geom_bar(position="dodge", stat="identity") +
  labs( title = "Number of times each haplotype occurs in a population ",
        y = "Number of populations",
        x = "sequence ",
        fill = "number of repeats for sequence based haplotype") +
  theme_bw()+
  theme(axis.text.x=element_blank())


#USE THIS GRAPH
#number of populations for haplotypes by locus
df_grp3 = df1 %>% group_by(Locus, sequence_) %>%
  summarise(counts = n_distinct(population))

df_grp3 = df_grp3 %>% group_by(Locus) %>%
  count(counts)

df_grp3$counts <- as.character(df_grp3$counts)

ggplot(df_grp3,aes(x = reorder(Locus, -n,sum), y=n, fill=counts)) + 
  geom_bar(position="stack", stat="identity") +
  labs( title = "Number of times alleles occurs in a population ",
        y = "Number of alleles",
        x = "STR Loci ",
        fill = "Number of Populations") +
  theme_bw()+
  scale_fill_brewer(palette="Paired") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = c(0.7, 0.7), 
        plot.title = element_text(hjust=0.5))+
  geom_text(aes(label= after_stat(y)),
            size =2,
            position = position_stack(0.5))

#USE THIS GRAPH
#recreating Laurence's graphs 
#displaying allelic diversity by loci 

#manipulating the df so that it is in the way that we want
df_grp4 = df1 %>% group_by(Locus, length_) %>%
 count(length_)

df_grp4 = df_grp4 %>% group_by(Locus) %>%
  count(Locus)

df_grp5 = df1 %>% group_by(Locus, sequence_) %>%
  count(sequence_)

df_grp5 = df_grp5 %>% group_by(Locus) %>%
  count(Locus)

df_merge <- merge(df_grp4,df_grp5,by="Locus")

df_merge <-df_merge%>%rename(length_based = n.x,
                    sequence_based = n.y)

df_merge$total_count <-  c(df_merge$length_based + df_merge$sequence_based)
df_merge<- df_merge[order(-df_merge$total_count),]

#graph the data 
df_merge %>% select(-total_count) %>% 
  gather(type, count, length_based:sequence_based) %>%
  ggplot(.,
         aes(x= reorder(Locus, -count), y=count, fill=forcats::fct_rev(type)))+
  geom_bar(stat="identity")+
  labs(title = "Allelic Diveristy",
       y = "Number of alleles observed",
       x= "STR Loci",
       fill= "Allele type")+
  theme_bw()+
  scale_fill_brewer(palette="Paired") +
  theme(legend.position = c(0.7, 0.8), 
        plot.title = element_text(hjust=0.5),
        axis.text = element_text(angle=90, vjust = 0.5, hjust = 1))+
  geom_text(aes(label= after_stat(y)),
            size =2,
            position = position_stack(0.5))

#Alluvial plot 
#group data by population and locus, count the number of haplotypes present
df_grp1 = df1 %>% group_by(Locus, population) %>%
  count(sequence_)


ggplot(data = df_grp1,
       aes(axis1 = Locus,
           axis2 = population)) +
  scale_x_discrete(limits = c("Locus", "population"), expand = c(.2, .05)) +
  geom_alluvium(aes(fill = sequence_)) +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.y=element_blank())+
  scale_fill_viridis_d() +
  ggtitle("Population specific variation")
