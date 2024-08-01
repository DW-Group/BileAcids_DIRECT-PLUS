################################################################################################################################################
#
# The improvement in body adiposity and lipid profiles induced by the Mediterranean diet intervention 
#   varied based on baseline levels of fecal bile acids.
#
##############################################################################################################################################

library(pacman)
pacman::p_load(tibble,ggpubr,grid,dplyr,psych,tidyverse)
pacman::p_load(ggplot2,vegan,ape,ggeffects,gee,geepack)
select <- dplyr::select

setwd(paste0(path,"/harvard/DIRECT_PLUS/Analysis/"))
load("/Results/DIRECT/dat_analysis.RData")

#data pre
gdata::keep(metadata, fmets_direct_BA, fmets_direct_BA_nm,BA_info_fmets,res_inter_Diet_BAs,sure=T)

p_23dg_anno <- function(x){
  x <- as.numeric(x)
  ifelse(x<0.001,'<0.001',
         ifelse(x>=0.001 & x<0.01,'<0.01',
                ifelse(x>=0.01,paste0('=',format(round(x,3),nsmall=3)),x)))};

res_inter  <- copy(res_inter_Diet_BAs) %>% 
  left_join(BA_info_fmets[,c('NAME_ABB','BA_FID')],by='NAME_ABB')

# Average in groups------
## BMI------
(sig_inter <- res_inter %>% filter(lipid_var=='BMI' & (as.numeric(q_fdr)<0.05|'q_fdr' %in% c('<0.05','<0.001'))) )
sig_FID_met <- unique(sig_inter$BA_FID); 

show_ba <- sig_inter %>% 
  filter(NAME_ABB %in% c('12-DHCA','TCA','TLCA-S','GCOA-S')) %>% 
  mutate(NAME=NAME_ABB, 
         NAME_ABB=factor(NAME_ABB,
                         levels=c('12-DHCA','TCA','TLCA-S','GCOA-S')),
         q_fdr=p_3dg(q_fdr)) %>% 
  arrange(NAME_ABB)
group_color <- c('#649BE7','#FF8911','#00A985')

plots <- list()
for (o in lipid_var[1]){
  p_inter <- show_ba[show_ba$lipid_var==o,]$'q_fdr';p_inter
  #1
  n=1
  m=show_ba$BA_FID[n] 
  df_input=fmets_direct_BA; 

  m_2g <- paste0(m,'_2g')
  df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                   'High BA levels','Low BA levels'), levels=c('Low BA levels','High BA levels'))
  box1 <- ggplot(data=df_input,
                 aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
    geom_boxplot(alpha=0.5)+
    stat_summary(fun=mean, geom='point', alpha=1,size=3,aes_string(group='Group2', color='Group2')) +
    stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group2',color='Group2')) + 
    scale_fill_manual(values=c('white','white'))+
    labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('q')[interaction]~.(p_inter[n])~')')),
         x='',y=bquote(paste(.(o),' (kg/m'^2~')')))+
    scale_color_manual(values=c("#728B92","#FF7F0E"))+
    scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                     breaks = c('0mon', '6mon', '12mon','18mon'),
                     labels=c('0', '6', '','18'))+
    facet_grid(as.formula(paste0('.~', m_2g)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          strip.background = element_rect(fill='grey90'),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(angle = 90,size=15),
          axis.text.x = element_text(size=15),text=element_text(size=15),
          legend.position="none",
          plot.title = element_text(hjust = 0,size=15),
          plot.subtitle = element_text(hjust = 0,size=15)); box1
  plots[[n]] <- box1
  
  #2)
  n=2
  m=show_ba$BA_FID[n]
  df_input=fmets_direct_BA;
  
  m_2g <- paste0(m,'_2g')
  df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                   'High BA levels','Low BA levels'), 
                            levels=c('Low BA levels','High BA levels'))
  box4 <- ggplot(data=df_input,
                 aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
    geom_boxplot(alpha=0.5)+
    stat_summary(fun=mean, geom='point', alpha=1,size=3,
                 aes_string(group='Group2', color='Group2')) +
    stat_summary(fun=mean, geom='line', size=1,
                 aes_string(group='Group2',color='Group2')) + 
    scale_fill_manual(values=c('white','white'))+
    labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('q')[interaction]~.(p_inter[n])~')')),
         x='',y=bquote(paste(.(o),' (kg/m'^2~')')))+
    scale_color_manual(values=c("#728B92","#FF7F0E"))+
    scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                     breaks = c('0mon', '6mon', '12mon','18mon'),
                     labels=c('0', '6', '','18'))+
    facet_grid(as.formula(paste0('.~', m_2g)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          strip.background = element_rect(fill='grey90'),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(angle = 90,size=15),
          axis.text.x = element_text(size=15),
          text=element_text(size=15),
          legend.position="none",
          plot.title = element_text(hjust = 0,size=15),
          plot.subtitle = element_text(hjust = 0,size=15)); box4
  plots[[n]] <- box4
  #3)
  n=3;m=show_ba$BA_FID[n]
  df_input=fmets_direct_BA;
  
  m_2g <- paste0(m,'_2g')
  df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                   'High BA levels','Low BA levels'), 
                            levels=c('Low BA levels','High BA levels'))
  box1 <- ggplot(data=df_input,
                 aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
    geom_boxplot(alpha=0.5)+
    stat_summary(fun=mean, geom='point', alpha=1,size=3,
                 aes_string(group='Group2', color='Group2')) +
    stat_summary(fun=mean, geom='line', size=1,
                 aes_string(group='Group2',color='Group2')) + 
    scale_fill_manual(values=c('white','white'))+
    labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('q')[interaction]~.(p_inter[n])~')')),
         x='',y=bquote(paste(.(o),' (kg/m'^2~')')))+
    scale_color_manual(values=c("#728B92","#FF7F0E"))+
    scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                     breaks = c('0mon', '6mon', '12mon','18mon'),
                     labels=c('0', '6', '','18'))+
    facet_grid(as.formula(paste0('.~', m_2g)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          strip.background = element_rect(fill='grey90'),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(angle = 90,size=15),
          axis.text.x = element_text(size=15),text=element_text(size=15),
          legend.position="none",
          plot.title = element_text(hjust = 0,size=15),
          plot.subtitle = element_text(hjust = 0,size=15)); box1
  plots[[n]] <- box1
  
  #4)
  n=4;m=show_ba$BA_FID[n]
  df_input=fmets_direct_BA;
  m_2g <- paste0(m,'_2g')
  df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_inputdf_input[,m]),na.rm=T),
                                   'High BA levels','Low BA levels'), levels=c('Low BA levels','High BA levels'))
  box2 <- ggplot(data=df_input,
                 aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
    geom_boxplot(alpha=0.5)+
    stat_summary(fun=mean, geom='point', alpha=1,size=3,aes_string(group='Group2', color='Group2')) +
    stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group2',color='Group2')) + 
    scale_fill_manual(values=c('white','white'))+
    labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('q')[interaction]~.(p_inter[n])~')')),
         x='',y=bquote(paste(.(o),' (kg/m'^2~')')))+
    scale_color_manual(values=c("#728B92","#FF7F0E"))+
    scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                     breaks = c('0mon', '6mon', '12mon','18mon'),
                     labels=c('0', '6', '','18'))+
    facet_grid(as.formula(paste0('.~', m_2g)))+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          strip.background = element_rect(fill='grey90'),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(angle = 90,size=15),
          axis.text.x = element_text(size=15),text=element_text(size=15),
          legend.position="none",
          plot.title = element_text(hjust = 0,size=15),
          plot.subtitle = element_text(hjust = 0,size=15)); box2
  plots[[n]] <- box2
}


combined_plot <- plot_grid(plotlist = plots, nrow = 2); # combined_plot
pdf('results/GEE/Inter_lipids_diet_str_BA_show.pdf',height=10,width=10)
combined_plot
dev.off()


# lipids
##TG---------
show_ba <- res_inter %>% 
  filter(lipid_var==lipid_var[2] & (as.numeric(P)<0.05)) %>% 
  mutate(P=p_3dg_anno(as.numeric(P)),
         NAME=NAME_ABB) %>% 
  filter(NAME_ABB %in% c('7keto-DCA','12-DHCA')); 
#1)average in 3 times

plots <- list()
for (o in lipid_var[2]){
  p_inter <- show_ba[show_ba$lipid_var==o,]$'P';p_inter
  for(n in 1: length(show_ba[show_ba$lipid_var==o,]$BA_FID)){
    m=show_ba$BA_FID[n] #12DHCA
    df_input=fmets_direct_BA; 
    m_2g <- paste0(m,'_2g')
    df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                     'High BA levels','Low BA levels'), 
                              levels=c('Low BA levels','High BA levels'))
    box1 <- ggplot(data=df_input,
                   aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
      geom_boxplot(alpha=0.5)+
      stat_summary(fun=mean, geom='point', alpha=1,size=3,
                   aes_string(group='Group2', color='Group2')) +
      stat_summary(fun=mean, geom='line', size=1,
                   aes_string(group='Group2',color='Group2')) + 
      scale_fill_manual(values=c("#728B92","#FF7F0E"))+
      labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('p')[interaction]~.(p_inter[n])~')')),
           x='',y='Triglycerides (mg/dL)')+
      scale_color_manual(values=c("#728B92","#FF7F0E"))+
      scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                       breaks = c('0mon', '6mon', '12mon','18mon'),
                       labels=c('0', '6', '','18'))+
      facet_grid(as.formula(paste0('.~', m_2g)))+theme_bw()+
      theme(panel.grid.major = element_blank(),
            strip.background = element_rect(fill='grey90'),
            axis.title.y = element_text(size=15), 
            axis.title.x = element_blank(),
            axis.text.y = element_text(angle = 90,size=15),
            axis.text.x = element_text(size=15),
            legend.position="none",
            plot.title = element_text(hjust = 0,size=15),
            plot.subtitle = element_text(hjust = 0,size=15)); box1
    plots[[n]] <- box1
  }
}

combined_plot <- plot_grid(plotlist = plots,ncol = 2); 
pdf('results/GEE/Inter_averageTG_diet2g_str_medianBA.pdf',height=5,width=10)
combined_plot
dev.off()


##TC/HDLc ratio--------------
show_ba <- res_inter %>% 
  filter(lipid_var==lipid_var[2] & (as.numeric(P)<0.05)) %>% 
  mutate(P=p_3dg_anno(as.numeric(P)),
         NAME=NAME_ABB); 

#1)average in 3 times

plots <- list()
for (o in lipid_var[3]){
  p_inter <- show_ba[show_ba$lipid_var==o,]$'P'; p_inter
  for(n in 1: length(show_ba[show_ba$lipid_var==o,]$BA_FID)){
    m=show_ba$BA_FID[n] 
    df_input=fmets_direct_BA; 
    
    m_2g <- paste0(m,'_2g')
    df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                     'High BA levels','Low BA levels'), 
                              levels=c('Low BA levels','High BA levels'))
    box1 <- ggplot(data=df_input,
                   aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
      geom_boxplot(alpha=0.5)+
      stat_summary(fun=mean, geom='point', alpha=1,size=3,
                   aes_string(group='Group2', color='Group2')) +
      stat_summary(fun=mean, geom='line', size=1,
                   aes_string(group='Group2',color='Group2')) + 
      
      scale_fill_manual(values=c("#728B92","#FF7F0E"))+
      labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('p')[interaction]~.(p_inter[n])~')')),
           x='',y='TC/HDLc ratio')+
      scale_color_manual(values=c("#728B92","#FF7F0E"))+
      scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                       breaks = c('0mon', '6mon', '12mon','18mon'),
                       labels=c('0', '6', '','18'))+
      facet_grid(as.formula(paste0('.~', m_2g)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            strip.background = element_rect(fill='grey90'),
            axis.title.y = element_text(size=15), axis.title.x = element_blank(),
            axis.text.y = element_text(angle = 90,size=15),
            axis.text.x = element_text(size=15),
            legend.position="none",
            plot.title = element_text(hjust = 0,size=15),
            plot.subtitle = element_text(hjust = 0,size=15)); box1
    plots[[n]] <- box1
  }
}

combined_plot <- plot_grid(plotlist = plots, ncol = 2); #combined_plot
pdf('results/GEE/Inter_averageTC_HDL_diet2g_str_medianBA.pdf',height=5,width=10)
combined_plot
dev.off()



##TC--------
lipid_var2 <- c("Cholesterol","LDLc")
nn=1
(show_ba <- res_inter %>% 
    filter(lipid_var==lipid_var[nn] & (as.numeric(P)<0.05)) %>%  
    mutate(P=p_3dg_anno(as.numeric(P)),NAME=NAME_ABB) %>% 
    filter(NAME_ABB %in% c('7keto-DCA','GCDCA'))); 
#1)average in 3 times

plots <- list()
for (o in lipid_var2[nn]){
  p_inter <- show_ba[show_ba$lipid_var==o,]$'P';p_inter
  for(n in 1: length(show_ba[show_ba$lipid_var==o,]$BA_FID)){
    m=show_ba$BA_FID[n] #
    df_input=fmets_direct_BA; 
    
    m_2g <- paste0(m,'_2g')
    df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                     'High BA levels','Low BA levels'), 
                              levels=c('Low BA levels','High BA levels'))
    box1 <- ggplot(data=df_input ,
                   aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
      geom_boxplot(alpha=0.5)+
      stat_summary(fun=mean, geom='point', alpha=1,size=3,
                   aes_string(group='Group2', color='Group2')) +
      stat_summary(fun=mean, geom='line', size=1,
                   aes_string(group='Group2',color='Group2')) + 
      
      scale_fill_manual(values=c("#728B92","#FF7F0E"))+
      labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('p')[interaction]~.(p_inter[n])~')')),
           x='',y='TC (mg/dL)')+
      scale_color_manual(values=c("#728B92","#FF7F0E"))+
      scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                       breaks = c('0mon', '6mon', '12mon','18mon'),
                       labels=c('0', '6', '','18'))+
      facet_grid(as.formula(paste0('.~', m_2g)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            strip.background = element_rect(fill='grey90'),
            axis.title.y = element_text(size=15), 
            axis.title.x = element_blank(),
            axis.text.y = element_text(angle = 90,size=15),
            axis.text.x = element_text(size=15),
            legend.position="none",
            plot.title = element_text(hjust = 0,size=15),
            plot.subtitle = element_text(hjust = 0,size=15)); box1
    plots[[n]] <- box1
  }
}

combined_plot <- plot_grid(plotlist = plots, nrow = 1); 
pdf(paste0('results/GEE/Inter_average_',lipid_var2[nn],'diet2g_str_medianBA.pdf'),
    height=5,width=10)
combined_plot
dev.off()

##LDLc--------
nn=2
(show_ba <- res_inter%>% 
    filter(lipid_var==lipid_var[nn] & (as.numeric(P)<0.05)) %>%  
    mutate(P=p_3dg_anno(as.numeric(P)),NAME=NAME_ABB) %>% 
    filter(NAME_ABB %in% c('GUDCA','GCDCA'))); 

#1)average in 3 times
plots <- list()
for (o in lipid_var2[nn]){
  p_inter <- show_ba[show_ba$lipid_var==o,]$'P';p_inter
  for(n in 1: length(show_ba[show_ba$lipid_var==o,]$BA_FID)){
    m=show_ba$BA_FID[n] 
    df_input=fmets_direct_BA; 
    
    m_2g <- paste0(m,'_2g')
    df_input[,m_2g] <- factor(ifelse(df_input[,m]>median(as.numeric(df_input[,m]),na.rm=T),
                                     'High BA levels','Low BA levels'), 
                              levels=c('Low BA levels','High BA levels'))
    box1 <- ggplot(data=df_input ,
                   aes_string(x='Time',y=o,fill='Group2',color='Group2'))+
      geom_boxplot(alpha=0.5)+
      stat_summary(fun=mean, geom='point', alpha=1,size=3,
                   aes_string(group='Group2', color='Group2')) +
      stat_summary(fun=mean, geom='line', size=1,
                   aes_string(group='Group2',color='Group2')) + 
      
      scale_fill_manual(values=c("#728B92","#FF7F0E"))+
      labs(title= bquote(paste(.(show_ba$NAME[n]),' ('~italic('p')[interaction]~.(p_inter[n])~')')),
           x='Time (month)',y='LDLc (mg/dL)')+
      scale_color_manual(values=c("#728B92","#FF7F0E"))+
      scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                       breaks = c('0mon', '6mon', '12mon','18mon'),
                       labels=c('0', '6', '','18'))+
      facet_grid(as.formula(paste0('.~', m_2g)))+
      theme_bw()+
      theme(panel.grid.major = element_blank(),
            strip.background = element_rect(fill='grey90'),
            axis.title = element_text(size=15),
            axis.text.y = element_text(angle = 90,size=15),
            axis.text.x = element_text(size=15),
            legend.position="bottom",
            plot.title = element_text(hjust = 0,size=15),
            plot.subtitle = element_text(hjust = 0,size=15)); box1
    plots[[n]] <- box1
  }
}

combined_plot <- plot_grid(plotlist = plots, nrow = 1); 
pdf(paste0('results/GEE/Inter_average_',lipid_var2[nn],'diet2g_str_medianBA.pdf'),height=6,width=10)
combined_plot
dev.off()



