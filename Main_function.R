library(tidyverse)
library(ggbeeswarm)


df <- readxl::read_xlsx("C:\\Users\\86156\\Desktop\\Qpcr_data.xlsx",col_names = F)

## Extract first line
header_group  <- df[1, -1]  ## remove replication
header_target <- df[2, -1]

## Fill first line
for(i in seq_along(header_group)) {
  if(is.na(header_group[i]) || header_group[i] == "NA") {
    header_group[i] <- header_group[i-1]
  }
}

## Generate new colname
new_names <- paste0(header_group, "_", header_target)
names(df) <- c("rep", new_names)

## Cut redundant column
df <- df[-c(1,2), ]

## Fill column
df <- df %>%
  fill(rep, .direction = "down")

## Change to long dataframe
df_long <- df %>%
  pivot_longer(
    cols = -rep,
    names_to = c("group", "target"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  mutate(value = as.numeric(value))

## Filter NoCq data
df_long<-df_long[which(!is.na(df_long$value)),]

## Summarise data
Data_summary<-df_long %>% group_by(rep,group,target) %>% summarise(Cq=mean(value))


NC<-"GAPDH"
Target<-"TTN"
Ctrl_group<-"Non-target control"
#Rank<-c("Non-target control","TTN-1","TTN-2","TTN-3")

## Colname
colnames(Data_summary)<-c("rep","group","target","cq")

## Calculate delta Cq
Neican<-Data_summary[which(Data_summary$target==NC),] %>% group_by(group) %>% summarise(mean=mean(cq))
Cq_frame<-merge(Data_summary[which(Data_summary$target==Target),],Neican,by="group") %>% mutate(dCq=cq-mean)

## Calculate baseline
Duizhao<-mean(Cq_frame$dCq[which(Cq_frame$group==Ctrl_group)])

## Calculate fold change
Res<-Cq_frame[,c(1,6)] %>% mutate(fold=log10(2^(Duizhao-dCq))+1) 

## Summarise results
Res_summary<-Res %>% group_by(group) %>% summarise(mean=mean(fold),sd=sd(fold))
Res_summary$group<-factor(Res_summary$group,levels=Rank)
ggplot(Res_summary)+
  geom_col(aes(x=group,y=mean),color=NA,fill="yellow",width=0.6)+
  geom_errorbar(aes(x=group,ymax=mean+sd,ymin=mean-sd),width=0.4)+
  geom_quasirandom(data=Res,
                   aes(x=group,y=fold),shape=21,fill="white",cex=2,dodge.width=0.5)+
  #scale_y_continuous(limits=c(0,7),expan=c(0,0))+
  #coord_fixed(ratio = 0.7)+ 
  #scale_y_continuous(limits=c(0,5),expan=c(0,0))+
  #coord_fixed(ratio = 1)+ 
  scale_y_continuous(limits=c(0,6),expan=c(0,0))+
  coord_fixed(ratio = 1.2)+ 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45,vjust=1,hjust=1))
ggsave("TH-1_res.pdf",height=4.3,width=4.3)
