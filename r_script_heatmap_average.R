heatmap_template
library(tidyverse)
library(ggplot2)
library(dplyr)
library(tidyr)

#Clean all environment 
rm(list = ls())


setwd("~/Analises R/r_array_huvec_120121/IL6_72hrs_vs_PBS_batch")

save.image(file = "heatmap_high_low_beta_121021")
load("heatmap_high_low_beta_121021")


#Upload the two csv files
sample_ids <- read.csv(file="./DMP_analysis_result.csv")
my_first_df <- read.csv(file="./raw_result_bmiq.csv")

#Merge two datasets
#select only the variables that are in the sample_ids file list in the SampleID column
my_second_df <- filter(my_first_df, SampleID %in% sample_ids$SampleID)
my_third_df <- merge(x= sample_ids, y= my_first_df, by=c("SampleID"))
head(my_third_df)

df <- my_second_df

####### calculations############
df_2 <- my_third_df
df_2 <- mutate(df_2, mean_col_pbs = rowMeans(select(df_2,  starts_with("PBS")), na.rm = TRUE)) #( IL6_each colunm - average controls)
head(my_third_df)
write.csv(df_2,file="./heatmap_template.csv")


###Load 
heatmap_test <- read.csv(file="./heatmap_template.csv")

#### Getting top 75 values
head(heatmap_test)
top_75beta <- heatmap_test%>% slice_max(deltaBeta, n = 75)

write.csv(top_50beta,file="./increase_75dmp.csv",quote=F,row.names = F)

top_75beta_heat <- top_75beta[,c( "gene",	"sIL6_exp2_1",	"sIL6_exp3_1",	"sIL6_exp2_2",	"sIL6_exp3_2",	"sIL6_exp1_1",	
                      "sIL6_exp1_2",	"sIL6_exp1_3", "sPBS_exp2_1", "sPBS_exp3_1",  "sPBS_exp2_2",  "sPBS_exp3_2",  "sPBS_exp1_1",  "sPBS_exp1_2")]

# Pivot make the matrix to heatmap
top_75beta_heat <- top_75beta_heat %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)

top_50beta_heat_pivot <- pivot_longer(data = top_50beta_heat,
                               cols = -c(1), 
                               names_to = "Sample", 
                               values_to = "Beta")

#New heatmap 
sapply(top_75beta_heat_pivot, class)
top_50beta_heat_pivot$Beta <- as.numeric(top_75beta_heat_pivot$Beta)

#Heat map
p_top50_heat <- ggplot(data = top_50beta_heat_pivot, mapping = aes(x = Sample, y = gene, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("beta value high on IL-6") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                        mid = "#f6f805",
                        high = "red",
                        guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.text.x = element_text(angle = 90)) 

p_top50_heat


################## lowest 75 heatmap###########

head(heatmap_test)
low_50beta <- heatmap_test%>% slice_min(deltaBeta, n = 75)
write.csv(low_50beta,file="./decrease_75dmp.csv",quote=F,row.names = F)



low_50beta_heat <- low_50beta[,c( "gene",	"sIL6_exp2_1",	"sIL6_exp3_1",	"sIL6_exp2_2",	"sIL6_exp3_2",	"sIL6_exp1_1",	
                                  "sIL6_exp1_2",	"sIL6_exp1_3", "sPBS_exp2_1", "sPBS_exp3_1",  "sPBS_exp2_2",  "sPBS_exp3_2",  "sPBS_exp1_1",  "sPBS_exp1_2")]

# Pivot make the matrix to heatmap
low_50beta_heat <- low_50beta_heat %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)

low_50beta_heat_pivot <- pivot_longer(data = low_50beta_heat,
                                      cols = -c(1), 
                                      names_to = "Sample", 
                                      values_to = "Beta")

#New heatmap 
sapply(low_50beta_heat_pivot, class)
low_50beta_heat_pivot$Beta <- as.numeric(low_50beta_heat_pivot$Beta)



p_low50_heat <- ggplot(data = low_50beta_heat_pivot, mapping = aes(x = Sample, y = gene, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.text.x = element_text(angle = 90)) 

p_low50_heat

#Lowest b-value - top 10
top_10_decr<- filter(heatmap_test, gene %in% c("TNFAIP8", "ARL15", "DGKE", "C6orf126", "RANBP9", "KCNQ1", "PHLDB3", "SLIT2", "FBXW7", "SFRS4"))

top_10_decr <- top_10_decr[,c( "gene",	"sIL6_exp2_1",	"sIL6_exp3_1",	"sIL6_exp2_2",	"sIL6_exp3_2",	"sIL6_exp1_1",	
                             "sIL6_exp1_2",	"sIL6_exp1_3", "sPBS_exp2_1", "sPBS_exp3_1",  "sPBS_exp2_2",  "sPBS_exp3_2",  "sPBS_exp1_1",  "sPBS_exp1_2")]

# Pivot make the matrix to heatmap
top_10_decr <- top_10_decr %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)

top_10_decr <- pivot_longer(data = top_10_decr,
                           cols = -c(1), 
                           names_to = "Sample", 
                           values_to = "Beta")

top_10_decr$Beta <- as.numeric(top_10_decr$Beta)
sapply(top_10_decr, class)


p_top_10_decr <- ggplot(data = top_10_decr, mapping = aes(x = Sample, y = gene, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  theme(text = element_text(size = 14)) +
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("Top 10 CpGs decrease B value") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.text.x = element_text(angle = 90)) 

p_top_10_decr

#Top_10_increase methylation

#Increase b-value - top 10
top_10_incre<- filter(heatmap_test, gene %in% c("UACA", "MFAP2", "GALNT7", "EPOR", "C10orf79", "RABGAP1L", "TMED10", "ZC3H12C", "ERC2", "MINA", "MYO3B", "MEF2A", "SSH1", "MAP3K10"))

top_10_incre <- top_10_incre[,c( "gene",	"sIL6_exp2_1",	"sIL6_exp3_1",	"sIL6_exp2_2",	"sIL6_exp3_2",	"sIL6_exp1_1",	
                               "sIL6_exp1_2",	"sIL6_exp1_3", "sPBS_exp2_1", "sPBS_exp3_1",  "sPBS_exp2_2",  "sPBS_exp3_2",  "sPBS_exp1_1",  "sPBS_exp1_2")]

# Pivot make the matrix to heatmap
top_10_incre <- top_10_incre %>% 
  mutate_if(is.numeric,as.character, is.factor, as.character)

top_10_incre <- pivot_longer(data = top_10_incre,
                            cols = -c(1), 
                            names_to = "Sample", 
                            values_to = "Beta")

top_10_incre$Beta <- as.numeric(top_10_incre$Beta)
sapply(top_10_incre, class)


p_top_10_incre <- ggplot(data = top_10_incre, mapping = aes(x = Sample, y = gene, fill = Beta)) +
  geom_tile() +
  xlab(label = "Sample") + 
  #remove x and y axis labels
  labs(x="",y="")+
  theme_grey(base_size=8)+
  theme(text = element_text(size = 14)) +
  #remove extra space
  scale_y_discrete(expand=c(0,0))+ 
  ggtitle("Top 10 CpGs increase B value") +
  scale_alpha( trans = "log" ) +
  scale_fill_gradient2(low = "#353436",
                       mid = "#f6f805",
                       high = "red",
                       guide = "colorbar") +
  coord_equal() + #square 
  theme(axis.text.x = element_text(angle = 90)) 

p_top_10_incre

################################################
#Selectin 25 rows (1 - 25)
newdata <- beta_heat %>% slice(1:800)
newdata2 <- beta_heat %>% slice(801:1601)
newdata3 <- beta_heat %>% slice(1602:2402)
newdata4 <- beta_heat %>% slice(2403:3203)
newdata5 <- beta_heat %>% slice(3203:3848)

### Frequency of DMP positions
library(ggplot2)
library(ggpubr)
theme_set(theme_pubr())

histo <- sample_ids

ggplot(sample_ids, aes(feature)) +
  geom_bar(fill = "#0073C2FF") +
  theme_pubclean()

## Compute the frequency 

df <- histo %>%
  group_by(feature) %>%
  summarise(counts = n())
df
head(df)

ggplot(df, aes(x = feature, y = counts)) +
  geom_bar(fill = "#0073C2FF", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  ylim(0, 200) + # xlim(5, 40) - Altera o eixo X 
  theme_pubclean()

ggsave("myplot.pdf")


