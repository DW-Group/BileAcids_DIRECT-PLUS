########################################################
#
# Dietary impact on fecal levels of bile acids
#
########################################################

pacman::p_load(tidyverse,gmodels,ggplot2,ggsci,corrplot,ggdensity,venn,vegan,ape,gee,geepack)
setwd(paste0(path,"/harvard/DIRECT_PLUS/Analysis/"))
load("/Results/Overview/dat_overview.RData")
gdata::keep(metadata, fmets_direct_BA, ba_diet_06,ba_diet_018, diet_inter_time, fmets_direct_BA_nm,sure=T)


##Fig S5_BAs change-Time---------------------
df=metadata %>% inner_join(fmets_direct_BA,by=c('sno'="CLIENT_SAMPLE_ID2",'time') )
data=metadata %>% inner_join(fmets_direct_BA,by=c('sno'="CLIENT_SAMPLE_ID2",'time') )
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
pdf('results/Overview/Fig_S1_BAs_Time_LOG_legend.pdf',height=0.5,width=4)
grid::grid.draw(ggpubr::get_legend(box_lg))
dev.off()


#fig S6----
ba_diet_06 <- ba_diet_06 %>%
  select(c('NAME_ABB','SUB_PATHWAY2','Category','Diet','Estimate','P','q_fdr')) %>%
  mutate(time='6m')

ba_diet_018 <- ba_diet_018 %>%
  select(c('NAME_ABB','SUB_PATHWAY2','Category','Diet','Estimate','P','q_fdr')) %>%
  mutate(time='18m')

diet_inter_time <- diet_inter_time %>% 
  mutate(NAME_ABB=ifelse(NAME_ABB %in% 'P/SBAs ratio','PBA/SBA ratio',NAME_ABB))
ba_diet_06[ba_diet_06$q_fdr<0.05,]
ba_diet_018[ba_diet_018$q_fdr<0.05,]
diet_inter_time[diet_inter_time$q<0.05,]

ba_diet <- rbind(ba_diet_06,ba_diet_018) %>% 
  unite('Diet_time',c('Diet','time'),sep='_')
ba_diet$NAME_ABB <- ifelse(ba_diet$NAME_ABB %in% 'P/SBAs ratio','PBA/SBA ratio',ba_diet$NAME_ABB)
ba_diet_est <- ba_diet %>% select(-c('P','q_fdr')) %>% 
  reshape2::dcast(NAME_ABB+SUB_PATHWAY2+Category~Diet_time) ; head(ba_diet_est)
ba2_lvl <- c(BA_info_fmets$NAME_ABB[-42],'PBA/SBA ratio')
ba_diet_est$NAME_ABB <- factor(ba_diet_est$NAME_ABB,levels=ba2_lvl)

ba_diet_est <- ba_diet_est[order(ba_diet_est$NAME_ABB),]

#Panel A: dietary impact on fecal bile acid levels
dheat <- ba_diet_est %>% column_to_rownames(var='NAME_ABB') %>% 
  select(c("PA+MED_6m","PA+GreenMED_6m",  "PA+MED_18m","PA+GreenMED_18m" )) %>% 
  setnames(c("MedDiet vs. HDG","Green-MedDiet vs. HDG","MedDiet vs. HDG","Green-MedDiet vs. HDG"))
ba_diet_p <- ba_diet %>% select(-c('Estimate','P')) %>% 
  reshape2::dcast(NAME_ABB+SUB_PATHWAY2+Category~Diet_time); head(ba_diet_p)
ba_diet_p$NAME_ABB <- factor(ba_diet_p$NAME_ABB,levels=ba2_lvl)
ba_diet_p <- ba_diet_p[order(ba_diet_p$NAME_ABB),] 
rownames(ba_diet_p) <- NULL

p_anno <-  ba_diet_p %>% column_to_rownames(var='NAME_ABB') %>% 
  select(c("PA+MED_6m","PA+GreenMED_6m","PA+MED_18m","PA+GreenMED_18m" )) %>% 
  setnames(c("MedDiet vs. HDG","Green-MedDiet vs. HDG","MedDiet vs. HDG","Green-MedDiet vs. HDG"))
#star marks
if (!is.null(p_anno)){
  ssmt <- p_anno <0.05
  p_anno[ssmt]<- '*'
  p_anno[!ssmt] <- ''
} else {
  p_anno <- F
};p_anno

ann_col <- list('Time'=c(rep('6mon',2),rep('18mon',2))) %>% 
  as.data.frame() %>% 
  mutate(Time=factor(Time,levels=c('6mon','18mon')))
time_col=c('#3FA446','#006E4A')
(row_color <- time_col[c(1,3)]);
names(row_color) <- unlist(unique(ann_col));row_color 

ann_row <- list(ba_diet_est$SUB_PATHWAY2,ba_diet_est$Category) %>% 
  as.data.frame() %>% 
  setnames(c('Sub pathway','Category'))
ann_row$Category <- factor(ann_row$Category,
                           levels=c('Unconjugated BA','Conjugated BA','Sum'))
col_color1 <- col_BA_cat; 
col_color2 <- col_BA_cat2
names(col_color1) <- unlist(unique(ba_diet_est$SUB_PATHWAY2));
names(col_color2) <- unlist(unique(ba_diet_est$Category))
names(time_col) <- c('0mon','6mon','18mon')
(ann_color <- list('Time'=time_col,'Sub pathway'=col_color1,'Category'=col_color2))

summary(dheat)

col=colorRamp2(c(-0.2,0,0.2),c("#0035A2",'white','#CB2C3E'))
ht <- ComplexHeatmap::pheatmap(as.matrix(dheat),
                               name='Beta',
                               legend_breaks=c(-0.2,-0.1,0,0.1,0.2),
                               cellheight=15,cellwidth=30,angle_col = '315',
                               fontsize_row = 12,fontsize_col=12,
                               show_rownames = T,show_colnames = T,border_color=NA,
                               cluster_rows=F,
                               cluster_cols=F,gaps_col = c(2),
                               annotation_row = ann_row,
                               annotation_col=ann_col,
                               annotation_colors=ann_color,
                               color= col, 
                               row_names_side = "left", column_names_side="bottom",
                               display_numbers = as.matrix(p_anno), fontsize_number=15,number_color='black',
                               legend=T,
                               heatmap_legend_param = list(legend_height=unit(2,"cm"), 
                                                           legend_direction='vertical',
                                                           title_position="topleft",
                                                           labels_gp = gpar(fontsize = 15), 
                                                           title_gp = gpar(fontsize = 15,
                                                                           fontface = "bold"))
); ht



pdf('results/GEE/Heatmap_beta_diet_time.pdf',width=15,height=20)
ht
dev.off()

## inter -p annotation----
ba_diet_est_anno <- ba_diet %>% select(-c('P','q_fdr')) %>% 
  reshape2::dcast(NAME_ABB+SUB_PATHWAY2+Category~Diet_time) %>% 
  left_join(diet_inter_time[,c('NAME_ABB','q')]); head(ba_diet_est_anno)

ba_diet_est_anno$NAME_ABB <- factor(ba_diet_est_anno$NAME_ABB,levels=rev(ba2_lvl))

ba_diet_est_anno <- ba_diet_est_anno[order(ba_diet_est_anno$NAME_ABB),]
dheat_anno <- ba_diet_est_anno %>%  
  select(c('NAME_ABB',"q" )) %>%  
  mutate(q_anno=case_when(q<0.05~'*',.default='')) %>% 
  mutate(q=log10(q)) 

ba_diet_p_ano <- dheat_anno
ba_diet_p_ano$NAME_ABB <- factor(ba_diet_p_ano$NAME_ABB,levels=ba2_lvl)
ba_diet_p_ano <- ba_diet_p_ano[order(ba_diet_p_ano$NAME_ABB),] 

q_anno <-  ba_diet_p_ano %>% 
  select(c("q" )) %>% 
  setnames(c("q for interaction"))

#star marks
if (!is.null(q_anno)){
  ssmt <- q_anno <0.05
  q_anno[ssmt]<- '*'
  q_anno[!ssmt] <- ''
} else {
  q_anno <- F
};q_anno

summary(ba_diet_est_anno$q)

ht_q <- ggplot(dheat_anno, aes(x =1 , y=NAME_ABB,fill=q)) +
  geom_tile(width = 1, height = 1, color='gray40',size=0.45)+
  geom_text(aes(label=q_anno), color="black", size=6,nudge_y = -0.2)+
  scale_fill_gradient2(low="#FFA859",mid='white',high='#005985', 
                       midpoint= -0.5)+
  scale_y_discrete(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  theme_few() +
  labs(title = "", y = "", x = "q",fill='q for interaction\n(log10)') +
  theme(axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=15,face='bold'),
        legend.title=element_text(size=15,face='bold'),
        legend.text=element_text(size=15),
        legend.position='none'); ht_q 

pdf('results/GEE/Heatmap_q_diet_time.pdf',width=1,height=9)
ht_q
dev.off()

#legend--
ht_q_lg <- ggplot(dheat_anno, aes(x =1 , y=NAME_ABB,fill=q)) +
  geom_tile(width = 1, height = 1, color='gray40',size=0.45)+
  geom_text(aes(label=q_anno), color="black", size=6,nudge_y = -0.2)+
  scale_fill_gradient2(low="#FFA859",mid='white',high='#005985', 
                       midpoint= -0.5)+
  scale_y_discrete(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  theme_few() +
  labs(title = "", y = "", x = "q",fill='q for interaction\n(log10)') +
  theme(axis.text = element_blank(),
        axis.title.x=element_text(size=15,face='bold'),
        legend.title=element_text(size=15,face='bold'),
        legend.text=element_text(size=15),
        legend.position='right'); ht_q_lg
pdf('results/GEE/Heatmap_q_diet_time_lg.pdf',width=2,height=2)
grid.draw(ggpubr::get_legend(ht_q_lg))
dev.off()


#panel B-boxplot--
plots <- list()
for (n in c(2,17,18)){
  for (i in BA_info_fmets$CHEM_ID2[n]){
    #print(i)
    data=fmets_direct_BA
    box <- ggplot(data,
                  aes_string(x='Time',y=i,fill='Group',color='Group'))+
      geom_boxplot(width=0.7,alpha=0.5,outlier.size = 0.5)+
      stat_summary(fun=mean, geom='point', size=5,alpha=1,aes_string(group='Group', color='Group')) +
      stat_summary(fun=mean, geom='line', size=1.5,aes_string(group='Group',color='Group')) +
      scale_fill_manual(values=c("#728B92",'#F47508',"#2A8B00"))+
      labs(title=BA_info_fmets$NAME_ABB[n],x='Time (Month)',y="Fecal levels (log10)")+
      scale_color_manual(values=c("#728B92",'#F47508',"#2A8B00"))+
      scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                       breaks = c('0mon', '6mon', '12mon','18mon'),
                       labels=c('0', '6', '','18'))+
      scale_y_continuous(expand=c(0,0.1))+
      theme_classic()+
      theme(panel.grid = element_blank(),
            axis.title = element_text(size=15),
            axis.text.y = element_text(angle = 90,size=15),
            axis.text.x = element_text(size=15),
            legend.position="none",
            plot.title = element_text(hjust =0,size=20,face='bold')); box
    plots[[i]] <- box
  }
}

pdf("results/GEE/Diet_BAs_time_box.pdf",width = 5, height = 15, onefile = F)
egg::ggarrange(plots,ncol=1)
dev.off()



#legend
box_lg <- ggplot(data,
                 aes_string(x='Time',y=i,fill='Group',color='Group'))+
  geom_boxplot(width=0.7,alpha=0.5,outlier.size = 0.5)+
  stat_summary(fun=mean, geom='point', size=5,alpha=1,
               aes_string(group='Group', color='Group')) +
  stat_summary(fun=mean, geom='line', size=1.5,aes_string(group='Group',color='Group')) +
  scale_fill_manual(values=c("#728B92",'#F47508',"#2A8B00"))+
  labs(title=BA_info_fmets$NAME_ABB[n],x='Time (Month)',y="Fecal levels (log10)")+
  scale_color_manual(values=c("#728B92",'#F47508',"#2A8B00"))+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=15),
        axis.text.y = element_text(angle = 90,size=15),
        axis.text.x = element_text(size=15),
        legend.position="bottom",
        plot.title = element_text(hjust =0,size=20,face='bold')); box_lg

pdf("results/GEE/Diet_BAs_time_box_legend.pdf",width = 5, height = 1, onefile = F)
grid::grid.draw(ggpubr::get_legend(box_lg))
dev.off()

