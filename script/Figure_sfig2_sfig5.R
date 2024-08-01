# Distribution of top 10 bile acids and microbial species
pacman::p_load(tidyverse,gmodels,ggplot2,ggsci,corrplot,ggdensity,venn,vegan,ape,gee,geepack)
setwd(paste0(path,"/harvard/DIRECT_PLUS/Analysis/"))
load("/Results/Overview/dat_overview.RData")
gdata::keep(metadata, fmets_direct_BA, BA_info_fmets,ba_rank,metadata_spe_fmets0,sure=T)


#*Overview plot-------
#*#ba 
ba_df0 <- metadata %>% filter(time==0) %>% 
  left_join(fmets_direct_BA,by=c('sno'="CLIENT_SAMPLE_ID2",'time') ) %>% 
  left_join(species_dat[,c('Sample_ID','s__Prevotella_copri')]) %>% as.data.frame() %>% 
  select(c('Sample_ID','sno','group','BMI',BA_info_fmets$CHEM_ID2)) %>% 
  arrange(group) %>% 
  mutate(Sample_ID <- factor(Sample_ID,levels=Sample_ID),
         Group=case_when(group==1~'HDG', group==2~'MedDiet',
                                            group==3~'Green-MedDiet',
                                            TRUE~as.character(group))) %>% 
  mutate(Group=factor(Group,levels=c('HDG','MedDiet','Green-MedDiet')))

#_1) Bile acids discribution--------
#ba data prepartion--
ba_df <- ba_df0 %>% select(c('sno',BA_info_fmets$CHEM_ID2[1:39])) %>% 
  column_to_rownames(var='sno') %>% rotate_df();  
ba_lg_mean <- map(ba_df[,BA_info_fmets$CHEM_ID2[1:39]],mean) %>% 
  data.frame() %>% setnames(BA_info_fmets$NAME_ABB) %>% 
  reshape2::melt(narm=F) %>% 
  left_join(BA_info_fmets[,c('NAME_ABB','SUB_PATHWAY2','Category')],
            by=c('variable'='NAME_ABB')); head(ba_lg_mean)
ba_rank=ba_lg_mean %>% arrange(desc(value)); ba_rank=head(ba_rank,10)

BA_top <- data.frame(apply(ba_df,2,log10)) %>% rownames_to_column(var='CHEM_ID2') %>% 
  left_join(BA_info_fmets[,c('CHEM_ID2','NAME_ABB')],by='CHEM_ID2') %>% 
  filter(NAME_ABB %in% ba_rank$variable) %>% 
  mutate(NAME_ABB=factor(NAME_ABB, levels=ba_rank$variable)) %>% arrange(NAME_ABB) %>% 
  select(-c('CHEM_ID2')) 

df_qt=BA_top %>% column_to_rownames(var='NAME_ABB') %>% rotate_df()
df_qt_out = apply(df_qt,2,function(x){quantile(x,probs=c(0.25,0.75))}) %>% 
  as.data.frame() %>% rotate_df() %>% setnames(c('Q1','Q3'))

df_qt_out = apply(df_qt_out,2,function(x){format(round(x,4),nsmall=4)})
df_qt_out <- data.frame(df_qt_out) %>% rownames_to_column(var='Bile_acids') %>% 
  mutate(IQR=paste0(Q1,' - ',Q3))
write.xlsx(df_qt_out,'results/Correlation/Quartiles_BAs.xlsx')


BA_top10_long <- reshape2::melt(BA_top,id.vars = 'NAME_ABB') %>% 
  mutate(NAME_ABB=factor(NAME_ABB,levels=rev(ba_rank$variable)))

## BA plot ------
##0) group---
group_color <- c("#0085CFCC","#E69206CC",'#2A8B00CC')
lca <- ba_df0[,c("F405",'Group','BMI','sno')] %>% 
  arrange(desc(F405)) %>% 
  arrange(Group) #%>% 
lca$sno <- factor(lca$sno,levels = lca$sno)
sample_order <- levels(lca$sno)
BA_top10_long$sno <- factor(BA_top10_long$sno,levels = levels(lca$sno))

p_group <- ggplot(lca,aes(x=sno,y=1,fill=factor(Group))) +
  geom_col(width = 1)+
  scale_fill_manual(values = group_color)+
  labs(y="Group")+
  theme_void()+
  theme(axis.title.y = element_text(size=15),
        legend.position = "none");p_group

#legend
p_group_lg <- ggplot(lca,aes(x=sno,y=1,fill=Group)) +
  geom_col(width = 1)+
  scale_fill_manual(values = group_color)+
  labs(y="Group")+
  theme_void()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right")
pdf("results/Correlation/Figs2_group_legend0.pdf", height = 1,width = 2)
grid.draw(ggpubr::get_legend(p_group_lg))
dev.off()

##2) bile acids---
p_ba <- ggplot(BA_top10_long, aes(x=sno,y=NAME_ABB)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Log 10 levels",
                       direction = -1, option = "D",end=1,begin = 0)+
  scale_y_discrete(expand=c(0,0)) +
  theme(legend.position = "none",
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.text.x = element_blank());p_ba
#legend
p_ba_lg <- ggplot(BA_top10_long, aes(x=sno,y=NAME_ABB)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Log 10 levels", na.value = "#6DCA55EE",
                       direction = -1, option = "D",end=1,begin = 0)+
  scale_y_discrete(expand=c(0,0)) +
  theme_void()+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),
        legend.position = "right",
        plot.title = element_blank(),
        legend.key.height =unit(0.5,'cm'),
        axis.title=element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()); p_ba_lg
pdf("results/Correlation/Figs2_ba_legend1.pdf", height = 2,width = 1)
grid.draw(ggpubr::get_legend(p_ba_lg))
dev.off()


#_2)spe_df Distributions --------------
spe_df <- species_dat %>% filter(Sample_ID %in% ba_df$Sample_ID) %>% 
  column_to_rownames(var='Sample_ID') %>% t %>% as.data.frame() 

spe_ave_relab <- 
  data.frame(species=rownames(spe_df),
             ave_relab=rowSums(spe_df,na.rm = T)/ncol(spe_df))
spe_ave_relab <- arrange(spe_ave_relab,desc(ave_relab))
head(spe_ave_relab)
(top10_spe <- head(spe_ave_relab$species,10))
spe_df <- spe_df %>%
  mutate(across(everything(),log10)) 
species_top10 <- spe_df[top10_spe,] %>% 
  rownames_to_column(var = "species") %>% 
  mutate(species=gsub('_',' ',gsub('s__','',species))) 
(spe_order <- factor(species_top10$species,levels = rev(species_top10$species)))
species_top10_long <- gather(species_top10,key = "Sample_ID",value = "reabun",colnames(species_top10)[-1]) #%>% 
tem=ba_df0[,c('sno','Sample_ID')] %>% 
  mutate(sno=factor(sno,levels=sno)) %>% 
  arrange(sno) %>% 
  mutate(Sample_ID,factor(Sample_ID,levels=Sample_ID))
sample_id_order=levels(tem$Sample_ID)
species_top10_long$Sample_ID <- factor(species_top10_long$Sample_ID,levels = sample_id_order)
species_top10_long$species <- 
  factor(species_top10_long$species,levels = spe_order)
species_top10_long <- species_top10_long %>% arrange(species,Sample_ID)
head(species_top10_long)

## species plot------
p_spe <- ggplot(species_top10_long, aes(Sample_ID,species)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10 scale", na.value = "gray80",
                       direction = -1, option = "G",begin=0,end=1)+
  labs(x='Samples')+
  scale_y_discrete(limits = rev(spe_order)) +
  theme(legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.position = "none",
        plot.title = element_blank(),
        axis.title.y=element_blank(),axis.title.x=element_text(size=15),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=15, face = "italic"),
        axis.text.x = element_blank()); p_spe

#legend
p_spe_lg <- ggplot(species_top10_long, aes(Sample_ID,species)) +
  geom_tile(aes(fill = reabun)) +
  scale_fill_viridis_c(name ="Relative abundance\nlog10", na.value = "gray50",
                       direction = -1, option = "G",end=1,begin = 0)+
  scale_y_discrete(limits = rev(spe_order)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 15),legend.box.just='center',
        legend.position = "right",
        legend.key.height =unit(0.5,'cm'),
        plot.title = element_blank(),
        axis.title=element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y=element_text(size=15, face = "italic"),
        axis.text.x = element_blank());p_spe_lg
pdf("results/Correlation/Fig0_spe_legend2.pdf", height = 2,width = 1)
grid.draw(ggpubr::get_legend(p_spe_lg))
dev.off()

##combine----
blank <- ggplot()+theme_void()
pdf("results/Correlation/Fig0_BA_sample.pdf",width = 10, height = 8, onefile = F)
egg::ggarrange(p_group,p_ba,blank,p_spe,ncol = 1,heights = c(0.1,2,0.05,2))
dev.off()

##Fig S5_BAs change-Time---------------------
df=metadata %>%   inner_join(fmets_direct_BA,by=c('sno'="CLIENT_SAMPLE_ID2",'time') )
data=metadata %>%   inner_join(fmets_direct_BA,by=c('sno'="CLIENT_SAMPLE_ID2",'time') )
data0=copy(data) %>% filter(time==0); setnames(data0,BA_info_fmets$CHEM_ID2,paste0(BA_info_fmets$CHEM_ID2,'_0'))
data6=copy(data) %>% filter(time==6); setnames(data6,BA_info_fmets$CHEM_ID2,paste0(BA_info_fmets$CHEM_ID2,'_6'))
data18=copy(data) %>% filter(time==18); setnames(data18,BA_info_fmets$CHEM_ID2,paste0(BA_info_fmets$CHEM_ID2,'_18'))

dfch=data0[,c('sno','Group',paste0(BA_info_fmets$CHEM_ID2,'_0'))] %>% 
  left_join(data6[,c('sno',paste0(BA_info_fmets$CHEM_ID2,'_6'))],by='sno') %>% 
  left_join(data18[,c('sno',paste0(BA_info_fmets$CHEM_ID2,'_18'))],by='sno')

for (i in BA_info_fmets$CHEM_ID2){
  ba0=paste0(i,'_0'); ba6=paste0(i,'_6'); ba18=paste0(i,'_18')
  ba60=paste0(i,'_60');ba180=paste0(i,'_180')
  dfch[,ba60]=log(dfch[,ba6],10)-log(dfch[,ba0],10)
  dfch[,ba180]=log(dfch[,ba18],10)-log(dfch[,ba0],10)
}

dfch2=bind_rows(dfch[,c('Group',paste0(BA_info_fmets$CHEM_ID2,'_60'))] %>% mutate(time='60') %>% 
                  setnames(c('Group',BA_info_fmets$CHEM_ID2,'Time')),
                dfch[,c('Group',paste0(BA_info_fmets$CHEM_ID2,'_180'))] %>% mutate(time='180') 
                %>% setnames(c('Group',BA_info_fmets$CHEM_ID2,'Time')) )

dbase=data.frame(dfch[,c('Group',paste0(BA_info_fmets$CHEM_ID2,'_60'))] %>% mutate(time='00') %>% 
                   setnames(c('Group',BA_info_fmets$CHEM_ID2,'Time')))
dbase[,BA_info_fmets$CHEM_ID2] <- apply(dbase[,BA_info_fmets$CHEM_ID2],2,function(x){x=0})
data=bind_rows(dfch2,dbase) %>% mutate(Time=factor(Time,levels=c('00','60','180')))


plots <- list()
time_col=c('#3FA446','#006E4A')
for (n in 1:(length(BA_info_fmets$CHEM_ID2))){
  for (i in BA_info_fmets$CHEM_ID2[n]){
    #print(i)
      box <- ggplot(data=data,
                  aes_string(x='Time',y=BA_info_fmets$CHEM_ID2[n],color='Group',fill='Group'))+
      stat_summary(fun=mean, geom='point', size=4.5,alpha=1,aes_string(group='Group')) +
      stat_summary(fun=mean, geom='line', size=1.5,aes_string(group='Group',color='Group')) +
      geom_hline(yintercept=0,color='gray60',linetype="dashed")+
      scale_fill_manual(values=group_color)+
      labs(title=BA_info_fmets$NAME_ABB[n],x='',y="")+
      scale_color_manual(values=group_color)+
      scale_x_discrete(limits = c('00','60', '180'), 
                       breaks = c('00','60', '180'),
                       labels=c('0mon','6-0mon', '18-0mon'))+
      theme_classic()+
      theme(panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(size=15),
            axis.text.x = element_blank(),
            legend.position="none",
            plot.title = element_text(hjust = 0,size=20,face='bold')); box
    plots[[n]] <- box
  }
}

for (n in 37:(length(BA_info_fmets$CHEM_ID2))){
  for (i in BA_info_fmets$CHEM_ID2[n]){
    #print(i)
    box <- ggplot(data=data,
                  aes_string(x='Time',y=BA_info_fmets$CHEM_ID2[n],color='Group',fill='Group'))+
      stat_summary(fun=mean, geom='point', size=4.5,alpha=1,aes_string(group='Group')) +
      stat_summary(fun=mean, geom='line', size=1.5,aes_string(group='Group',color='Group')) +
      geom_hline(yintercept=0,color='gray60',linetype="dashed")+
      scale_fill_manual(values=group_color)+
      labs(title=BA_info_fmets$NAME_ABB[n],x='',y="")+
      scale_color_manual(values=group_color)+ 
      scale_x_discrete(limits = c('00','60', '180'), 
                       breaks = c('00','60', '180'),
                       labels=c('0mon','6mon-0mon', '18mon-0mon'))+
      theme_classic()+
      theme(panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.text.y = element_text(size=15),
            axis.text.x = element_text(size=20,angle=30,hjust=0.8,vjust=0.9,color='black'),
            legend.position="none",
            plot.title = element_text(hjust = 0,size=20,face='bold')); box
    plots[[n]] <- box
  }
}

lay_ba <- lay_new(
  mat = matrix(1:42, ncol = 6),widths = rep(1,6),heights = c(1,1,1,1,1,1,1.5))  
pdf('results/Missing/Fig_S5_BAs_Time_LOG.pdf',height=20,width=20)
plot_grid(plotlist = plots,ncol=6,
          rel_widths = rep(1,6),rel_heights = c(rep(1,6),1.4))
dev.off()



box_lg <- ggplot(data=data,
                 aes_string(x='Time',y=BA_info_fmets$CHEM_ID2[n],color='Group',fill='Group'))+
  stat_summary(fun=mean, geom='point', size=3,alpha=1,aes_string(group='Group')) +
  stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group',color='Group')) +
  scale_fill_manual(values=group_color)+
  labs(title=BA_info_fmets$NAME_ABB[n],x='',y="")+  
  scale_color_manual(values=group_color)+
  theme_classic()+
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        legend.text= element_text(size=15),
        legend.position="bottom",
        legend.title = element_text(size=15,face='bold')); box_lg
pdf('results/Missing/Fig_S1_BAs_Time_LOG_legend.pdf',height=0.5,width=4)
grid::grid.draw(ggpubr::get_legend(box_lg))
dev.off()





