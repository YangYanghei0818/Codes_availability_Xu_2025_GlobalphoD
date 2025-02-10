library(data.table)

otu_unrare <- fread("feature-table_unrarefied.tsv",sep = "\t",skip = 1)
#function to calculate relative abundance of OTU_table and selecting ubiquitous OTUS

filter_ubiquitous <- function(otu_table, ubiquitous_rate) {
  # OTU_table must be sample in col and OTUID in row
  # the first column should be "OTU_ID" or "ASV_ID"
  
  stopifnot(is.double(ubiquitous_rate), ubiquitous_rate >= 0 && ubiquitous_rate <= 1)
  
  otu_table %>%
    rename(OTU_ID = 1) %>%  # Assuming the first column is OTU_ID
    dplyr::select(-OTU_ID) %>%  # Remove the OTU_ID column for calculation (select for column)
    mutate(across(everything(), ~ . / sum(.))) %>%  # Normalize columns 
    filter(rowSums(. > 0) >= ncol(.) * ubiquitous_rate) # filter for row
}


filter_ubiquitous_absolute <- function(otu_table, ubiquitous_rate) {
  # OTU_table must be sample in col and OTUID in row
  # the first column should be "OTU_ID" or "ASV_ID"
  
  stopifnot(is.double(ubiquitous_rate), ubiquitous_rate >= 0 && ubiquitous_rate <= 1)
  
  otu_table %>%
    rename(OTU_ID = 1) %>%  # Assuming the first column is OTU_ID
    dplyr::select(-OTU_ID) %>%  # Remove the OTU_ID column for calculation
    filter(rowSums(. > 0) >= ncol(.) * ubiquitous_rate) # filter for row
}


       
# Calculating 

otu_unrare1 <- as.data.frame(otu_unrare)
rownames(otu_unrare1) <- otu_unrare1[,1]

otu_unrare2 <- otu_unrare1[,2:ncol(otu_unrare1)]
OTU_uniquity <- rowSums(otu_unrare2 > 0) %>% as.data.frame() 

# ubiquity rate

OTU_ubiquity<-(OTU_uniquity/ncol(otu_unrare2))*100  
names(OTU_ubiquity)<-"ubiquity"
options(scipen = 999)



otu_morethan100 <- otu_unrare1 %>% 
  filter_ubiquitous(100/3115)

