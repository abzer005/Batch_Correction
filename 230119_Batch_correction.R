# Packages ----------------------------------------------------------------
library(magrittr)
library(vegan)
library(tidyverse)
library(ggpubr)

setwd(normalizePath(readline("Enter the path of the folder with input files: "),"/",mustWork=FALSE))

#Load input data:-----------------------------------------------------------------
ft <- read.csv("DOM_Tissue_GNPS_quant.csv",check.names = F)
md <- read.csv("DOM_Tissue_meta.csv", check.names = F)

colnames(ft) <- gsub(' Peak area','',colnames(ft)) #Removing Peak area extensions from the column names of ft
ft <- ft[,colSums(is.na(ft))<nrow(ft)] #Removing if any NA columns present in the file

new_ft <- ft
new_md <- md
#Changing the row names of the files
rownames(new_ft) <- paste(new_ft$'row ID', 
                          round(new_ft$'row m/z',digits = 3),
                          round(new_ft$'row retention time',digits = 3), sep = '_') 
new_ft <- new_ft[,grep('mzML',colnames(new_ft))] #Picking only the files with column names containing 'mzML'

#Transposing the ft
ft_t <- as.data.frame(t(new_ft))
ft_t <- ft_t %>% mutate_all(as.numeric)  #converting all values to numeric
colnames(ft_t) <- paste0('X',colnames(ft_t)) # since our colnames are numbers, it is better to have a character suffix

ft_t <- ft_t[order(rownames(ft_t)),,drop=F] #ordering the ft_t by its row names
md <- md[order(md$filename),,drop=F] #ordering the md by its row names
identical(md$filename,rownames(ft_t)) #should return TRUE now

table(md$ATTRIBUTE_Tray) # This attribute contains the batch info

#Picking only Tray 1-5:
md_batch <- md %>% filter(ATTRIBUTE_Tray == 'Tray1'|
                            ATTRIBUTE_Tray == 'Tray2'|
                            ATTRIBUTE_Tray == 'Tray3'|
                            ATTRIBUTE_Tray == 'Tray4'|
                            ATTRIBUTE_Tray == 'Tray5')

ft_batch <- ft_t  %>% filter(row.names(ft_t) %in% md_batch$filename)

identical(md_batch$filename,rownames(ft_batch))

#get PCoA on raw data:
ft_batch2 <- ft_batch + 0.1
ft_s <- scale(ft_batch2, scale = T, center = F)
identical(rownames(ft_s),md_batch$filename)

distm <- vegdist(ft_s,  method = "bray") # calculating a distance matrix
PcoA <- cmdscale(distm, k = 10, eig = T, add = T)
PcoA_points <- as.data.frame(PcoA$points)
variance <- round(PcoA$eig*100/sum(PcoA$eig),1)
names(PcoA_points)[1:10] <- paste0('PCoA',1:10)
identical(rownames(PcoA_points), md_batch$filename)


p1 <- ggplot(PcoA_points, aes(x = PCoA1, 
                        y = PCoA2, 
                        colour = as.factor(md_batch$ATTRIBUTE_Tissue_Study), 
                        label = row.names(PcoA))) +
  geom_point(size=2, alpha=0.6) +
  scale_colour_manual(values = c('orange','darkgreen','red','blue','black')) +
  labs(title = "Attribute: Tissue Study") + 
  xlab(paste('PCoA1: ',variance[1],'%', sep = ' ')) +
  ylab(paste('PCoA2: ',variance[2],'%', sep = ' ')) +
  theme(legend.position = 'bottom',legend.title=element_blank())
p1

p2 <- ggplot(PcoA_points, aes(x = PCoA1, 
                              y = PCoA2, 
                              colour = as.factor(md_batch$ATTRIBUTE_Tray), 
                              label = row.names(PcoA))) +
  geom_point(size=2, alpha=0.6) +
  scale_colour_manual(values = c('orange','darkgreen','red','blue','black')) +
  labs(title = "Attribute: Tray") + # labs(title ='Given title')
  xlab(paste('PCoA1: ',variance[1],'%', sep = ' ')) +
  ylab(paste('PCoA2: ',variance[2],'%', sep = ' ')) +
  theme(legend.position = 'bottom',legend.title=element_blank())

p2 

p <- ggarrange(p1, p2, labels = c("A", "B"),ncol = 2, align = 'v')

annotate_figure(p,top = text_grob("PCoA Scores Plots on scaled raw data before Batch correction \n using Bray-Curtis Dissimilarity Index", 
                                  face = "bold", size = 14))

ggsave('19_01_23_scores_plot_before_Inter_batch.svg',width=13, height=10)

# For inter-batch correction:-------------------------------------------------------------------

#merging metadata (new_md) and transposed feature table based on the sample names
ft_merged <- merge(md_batch,ft_batch, by.x= "filename", by.y=0,all.x=TRUE) #by.y =0 indicates the rownames of ft_t
ft_merged2 <- ft_merged %>% select(`filename`,`ATTRIBUTE_Tray`,starts_with("X")) 

#STEP 1: Calculate the overall mean of the feature 
fm <- as.data.frame(rbind(colMeans(ft_merged2[,-(1:2)]))) #getting the columnwise mean for ft_merged2 except its 1st 2 columns

#STEP 2: Batch-specific feature mean
bm <- ft_merged2[,-1] %>%  #excluding filename column as we are geting only batchwise mean value
  group_by(`ATTRIBUTE_Tray`) %>%  # grouping them by Batch
  summarise_all(mean) %>% # getting column-wise mean
  column_to_rownames('ATTRIBUTE_Tray') %>%
  as.data.frame() # storing it as dataframe

#STEP 3: inter-batch corrected data
##The intensities in each batch are then divided by the batch mean and multiplied by the overall mean.

batch_df <- ft_merged2 %>%
  group_split(`ATTRIBUTE_Tray`) %>% #group_split splits & stores the batchwise info as individual dataframes inside a list
  lapply(., function(x) {  
    x <- column_to_rownames(x,'filename') # then, we make "filename" as the rownames of each dataframe within the list
  }) 

sapply(batch_df, dim) # gives the dimension of each list element columnwise.

ib <- list()
for (i in 1:length(batch_df)){
  ib[[i]] <- sweep((batch_df[[i]][,-1]), 2, as.numeric(bm[i,]+1), "/") #dividing each batch dataframe by batchwise feature-mean
  ib[[i]] <- sweep(ib[[i]], 2, as.numeric(fm+1), "*") # multiplying by overall mean
}

ib <- bind_rows(ib) #binding all the list elements together

#Testing with PCoA

Final <- ib + 1 #To avoid 0 values in scaled data, we replace 0s with 1s
Final_S <-scale(Final, scale = T, center = F) #scaling the table
Final_S <- Final_S[order(rownames(Final_S)),,drop=F] #ordering by rownames
identical(rownames(Final_S),md_batch$filename)

distm <- vegdist(Final_S, method = "bray")# compute distance
PcoA <- cmdscale(distm, k = 2, eig = T, add = T)
PcoA_points <- as.data.frame(PcoA$points)
variance <- round(PcoA$eig*100/sum(PcoA$eig),1)
names(PcoA_points)[1:2] <- c('PCoA1', 'PCoA2')
identical(rownames(PcoA_points), md_batch$filename)

g1 <- ggplot(PcoA_points, aes(x = PCoA1, y = PCoA2, colour = as.factor(md_batch$ATTRIBUTE_Tissue_Study), label = row.names(PcoA))) +
  geom_point(size=2, alpha=0.6) +
  scale_colour_manual(values = c('orange','darkgreen','red','blue','black')) +
  xlab(paste('PCoA1: ',variance[1],'%', sep = ' ')) +
  ylab(paste('PCoA2: ',variance[2],'%', sep = ' ')) +
  labs(title = "Attribute: Tissue Study") + # labs(title ='Given title')
  theme(legend.position = 'bottom',legend.title=element_blank())
  

g2 <- ggplot(PcoA_points, aes(x = PCoA1, y = PCoA2, colour = as.factor(md_batch$ATTRIBUTE_Tray), label = row.names(PcoA))) +
  geom_point(size=2, alpha=0.6) +
  scale_colour_manual(values = c('orange','darkgreen','red','blue','black')) +
  xlab(paste('PCoA1: ',variance[1],'%', sep = ' ')) +
  ylab(paste('PCoA2: ',variance[2],'%', sep = ' ')) +
  labs(title = "Attribute: Tray") + # labs(title ='Given title')
  theme(legend.position = 'bottom',legend.title=element_blank())


g <- ggarrange(g1, g2, labels = c("A", "B"),ncol = 2, align = 'v')
annotate_figure(g, top = text_grob("PCoA Scores Plots on scaled raw data after Inter-Batch correction \n using Bray-Curtis Dissimilarity Index", 
                                   face = "bold", size = 14))

ggsave('19_01_23_scores_plot_after_Inter_batch.svg',width=13, height=10)



#Intra-Batch correction:----------------------------------------------------------------------------------

ft_split <- ft_merged %>%  group_split(ATTRIBUTE_Tray)

#to plot tryptophan
trp <- ft %>% filter(`row m/z` >=204, `row m/z` < 205)

#after looking at the different mz values, lets pick mz 204.1380 --> ID 54258
grep('X54258',colnames(ft_merged),value = TRUE) # lets pick either 'X930_204.957_0.562'or "X54258_204.138_9.808"

for(i in 1:5){
  x <- ggplot(ft_split[[i]], 
              aes(x=`ATTRIBUTE_Injection_Order`, 
                  y=`X930_204.957_0.562`, #paste the y axis name from the previous cell output
                  color=`ATTRIBUTE_Tissue_Study`, 
                  group=`ATTRIBUTE_Tissue_Study`)) + 
    geom_point(size=2.5, alpha=0.9)+
    geom_line() + labs(title = paste('Tray ',i,' Before intra batch correction'))
  
  ggsave(paste('Tray_',i,'_before_intra_batch_corr.svg'),width=10, height=10)
  show(x)
}

for(i in 1:5){
  x <- ggplot(ft_split[[i]], 
         aes(x=`ATTRIBUTE_Injection_Order`, 
             y=`X54258_204.138_9.808`, 
             color=`ATTRIBUTE_Tissue_Study`, 
             group=`ATTRIBUTE_Tissue_Study`)) + #paste the y axis name from the previous cell output
    geom_point(size=2.5, alpha=0.9)+
    geom_line() + labs(title = paste('Tray ',i,' Before intra batch correction'))
  
  show(x)
}


# Trying to perform intra-batch only on the 1st batch:
batch_corrected <- ft_split[[1]] %>% 
  mutate(is_QC = (ATTRIBUTE_Tissue_Study=="Pooled_QC")) %>% #here, I provided the QC information
  mutate(data = map(data,
                    ~ suppressWarnings(mutate(..1, # I do not understand '..1'
                                              into_driftcor = drift_correct_with_LOESS_on_QCs_vec_p(`X54258_204.138_9.808`,
                                                                                                          is_QC,
                                                                                                          ATTRIBUTE_Injection_Order)$x.corrected))))
#---------------------------------------------------------