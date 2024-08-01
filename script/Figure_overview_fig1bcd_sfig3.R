############################################################################################################
# Overview of the multi-omics data in the 18 month randomized controlled trial：
#     fecal bilea acid, body adopisty, serum lipid biomarkers, and gut microbiome
############################################################################################################

pacman::p_load(tidyverse,gmodels,ggplot2,ggsci,corrplot,ggdensity,venn,vegan,ape,gee,geepack)
setwd(paste0(path,"/harvard/DIRECT_PLUS/Analysis/"))
load("/Results/Overview/dat_overview.RData")
gdata::keep(metadata, fmets_direct_BA, fmets_direct_BA_nm,fpmets_name,direct_tax,metadata_spe_fmets0,sure=T)

#Panel B: Time trend of outcome--------
#data pre
outcome=c("BMI","Triglycerides","TC_HDL")
summary(metadata[,outcome])

metadata <- within(metadata,{
  Time <- factor(case_when(time==0~'0m',
                           time==6~'6m',
                           time==18~'18m',
                           .default = time),
                 levels=c('0mon','6mon','18mon'))
  Group <- factor(case_when(group==1~'HDG',
                            group==2~'MedDiet',
                            group==3~'Green-MedDiet',.default = group),
                  levels=c('HDG','MedDiet','Green-MedDiet'));
  sex <- case_when(sex==0~'Male',sex==1~'Female',.default =sex)
  sno <- paste0('S',sno)
  TC_HDL <- Cholesterol/HDLc
  }
)


#1)--
group_color <- c('#649BE7','#FF8911','#00A985')
o='BMI'
box1 <- ggplot(data=df_gee,
               aes_string(x='Time',y=o,fill='Group',color='Group'))+
  geom_boxplot(aes(x=Time),width=0.77,size=0.43,alpha=0.8)+
  stat_summary(fun=mean, geom='point', alpha=1,size=2.5,
               aes_string(group='Group', color='Group')) +
  stat_summary(fun=mean, geom='line', size=1,
               aes_string(group='Group',color='Group')) + 
  labs(y='',x='Time (Month)',title=o)+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  scale_y_continuous(expand = c(0.01,0.0))+
  annotate('text',x=2,y=max(df_gee[,o]),size=4,
           label=bquote(paste(italic('p')~.(p_23dg(pt1$'P')))))+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90,size=14),
        axis.title.x = element_text(size=15),
        axis.text = element_text(size=15),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size=17),
        plot.subtitle = element_text(hjust = 0,size=13)); box1


#2)--
o='Triglycerides'
box2 <- ggplot(data=df_gee,
               aes_string(x='Time',y=o,fill='Group',color='Group'))+
  geom_boxplot(aes(x=Time),width=0.77,size=0.43,alpha=0.8)+
  stat_summary(fun=mean, geom='point', alpha=1,size=2.5,
               aes_string(group='Group', color='Group')) +
  stat_summary(fun=mean, geom='line', size=1,
               aes_string(group='Group',color='Group')) + 
  labs(y='',x='Time (Month)',title=o)+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  scale_y_continuous(expand = c(0.01,0.0))+
  annotate('text',x=2,y=max(df_gee[,o]),size=4,
           label=bquote(paste(italic('p')~.(p_23dg(pt2$'P')))))+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90,size=14),
        axis.title.x = element_text(size=15),
        axis.text = element_text(size=15),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size=17),
        plot.subtitle = element_text(hjust = 0,size=13)); box2


#3)--
o='TC_HDL'
box3 <- ggplot(data=df_gee,
               aes_string(x='Time',y=o,fill='Group',color='Group'))+
  geom_boxplot(aes(x=Time),width=0.77,size=0.43,alpha=0.8)+
  stat_summary(fun=mean, geom='point', alpha=1,size=2.5,
               aes_string(group='Group', color='Group')) +
  stat_summary(fun=mean, geom='line', size=1,
               aes_string(group='Group',color='Group')) + 
  labs(y='',x='Time (Month)',title=o)+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  scale_y_continuous(expand = c(0.01,0.0))+
  annotate('text',x=2,y=max(df_gee[,o]),size=4,
           label=bquote(paste(italic('p')~.(p_23dg(pt3$'P')))))+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90,size=14),
        axis.title.x = element_text(size=15),
        axis.text = element_text(size=15),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size=17),
        plot.subtitle = element_text(hjust = 0,size=13)); box3

pdf("results/Overview/Fig1_time_trend_Outcome.pdf",width = 12, height = 4, onefile = F)
egg::ggarrange(box1,box2,box3,ncol = 3, widths = c(1,1,1))
dev.off()


#legend--
box_lg <- ggplot(data=df_gee,
                 aes_string(x='Time',y=o,fill='Group',color='Group'))+
  geom_boxplot(aes(x=Time),width=0.77,size=0.43,alpha=0.8)+
  stat_summary(fun=mean, geom='point', alpha=1,size=2.5,
               aes_string(group='Group', color='Group')) +
  stat_summary(fun=mean, geom='line', size=1,
               aes_string(group='Group',color='Group')) + 
  labs(y='',x='Time (Month)',title=o)+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  theme_void()+
  theme(legend.position="right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 14)); box_lg


pdf("results/Overview/Fig1_time_trend_Outcome_legend.pdf", height = 1.5,width = 2)
grid.draw(ggpubr::get_legend(box_lg))
dev.off()


##-P-outcome-time--------------
#BMI
fit2 <- geeglm(BMI~Group*Time+age+sex+antibio+metformin+ Lipid_lowering,
               id=ID,data=metadata,
               corstr = constr,family = fam); summary(fit2)
fit1 <- geeglm(BMI~Group+Time+age+sex+antibio+metformin+ Lipid_lowering,
               id=ID,data=metadata,
               corstr = constr,family = fam); summary(fit1)
pt1 <- anova(fit1,fit2,test='lrtest');p_23dg(pt1$'P') #p<0.001

#TG--
fit2 <- geeglm(Triglycerides~Group*Time+age+sex+antibio+metformin+ Lipid_lowering,
               id=ID,data=metadata,
               corstr = constr,family = fam); summary(fit2)
fit1 <- geeglm(Triglycerides~Group+Time+age+sex+antibio+metformin+ Lipid_lowering,
               id=ID,data=metadata,
               corstr = constr,family = fam); summary(fit1)
pt2 <- anova(fit1,fit2,test='lrtest');p_23dg(pt2$'P')#p:0.046

#TC_HDL
fit2 <- geeglm(TC_HDL~Group*Time+age+sex+antibio+metformin+ Lipid_lowering,
               id=ID,data=metadata,
               corstr = constr,family = fam); summary(fit2)
fit1 <- geeglm(TC_HDL~Group+Time+age+sex+antibio+metformin+ Lipid_lowering,
               id=ID,data=metadata,
               corstr = constr,family = fam); summary(fit1)
pt3 <- anova(fit1,fit2,test='lrtest');p_23dg(pt3$'P')#p:<0.01


#Panel C: correlations among fecal bile acids-----
##cor: fmets-fmets----
fmets_direct_BA_qc_nm <- fmets_direct_BA_qc_nm %>% setnames(fpmets_name$NAME_ABB) 

df_cor <- fmets_direct_BA_qc_nm %>% filter(time==0) %>% 
  select(c('CLIENT_SAMPLE_ID2',BA_info_fmets$NAME_ABB)) %>% 
  column_to_rownames(var='CLIENT_SAMPLE_ID2')
df_cor <- df_cor[,levels(ba_lg_mean$variable)]
cor_ff <- corr.test(df_cor[,levels(ba_lg_mean$variable)], 
                    method='spearman',adjust="BH",alpha=0.05)
cor_ff.r <- cor_ff$r %>% data.frame() %>% setnames(colnames(cor_ff$r));
#quick look
tem <- cor_ff.r %>% as.data.frame(); 
tem[,1:ncol(tem)] <- map(tem[,1:ncol(tem)], function(x){ifelse(as.numeric(x)>=0.9999,NA,x)}); 
summary(tem[1:nrow(tem),1:ncol(tem)])
max(tem[1:nrow(tem),1:ncol(tem)],na.rm=T);# 0.94
min(tem[1:nrow(tem),1:ncol(tem)], na.rm=T);#-0.52
#CA-7keto-DCA
tem <- cor_ff$p %>% as.data.frame();
p.adjust(tem$CA,method='BH')[26];cor_ff.r[26,'CA'] #q<0.001
#7keto-DCA-UCA
tem <- cor_ff$p %>% as.data.frame();
p.adjust(tem$UCA,method='BH')[26];cor_ff.r[26,'UCA']#q<0.001


#1)triangle corrplot
get_upper_tri <- function(cormat){
  cormat[upper.tri(cormat,diag = F)]<- NA
  return(cormat)
}
(upper_tri <- as.matrix(rev(get_upper_tri(cor_ff.r)))); 
banm=rev(rownames(upper_tri))
res=data.frame(upper_tri[banm,banm]) %>% 
  rownames_to_column(var='Bile acids')
melted_cormat <- reshape2::melt(upper_tri,narm = T)


p3 <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+ 
  geom_tile(color = "white")+ 
  scale_fill_gradient2(low = "#2166AC", high = "#B10C1F", mid = "white", 
                       midpoint = 0, limit = c(-1,1), 
                       space = "Lab", na.value = "white",
                       name="Spearman  \nCorrelation  ") +
  theme_minimal()+ 
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position='none')+ 
  coord_fixed(ratio=1) ;print(p3)

#2.1) col annotation 1
cor_colnames <- BA_info_fmets %>% filter(NAME_ABB %in% colnames(df_cor))

#check names
identical(levels(ba_lg_mean$variable),colnames(df_cor))
col_anno <- cor_colnames[,c("NAME_ABB","SUB_PATHWAY2")] %>% 
  mutate(y=1,NAME_ABB=factor(NAME_ABB,level=rev(levels(ba_lg_mean$variable))),
         SUB_PATHWAY2=factor(SUB_PATHWAY2,levels=c("Primary BA","Secondary BA","Sum"))) %>% 
  arrange(col_anno$NAME_ABB)


p1 <- ggplot(col_anno, aes(x=NAME_ABB, y=y)) + 
  geom_col(aes(fill=SUB_PATHWAY2),width=1,color="white",linewidth=0.2)+#) +
  scale_fill_manual(values=col_BA_cat)+
  labs(fill="Sub Pathway",y='Sub Pathway') +
  scale_y_continuous(breaks=c(0,1), limits=c(0,1),expand=c(0,0)) +
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(panel.border = element_rect(size=0.5), 
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y.left = element_text(angle=0,vjust=0.5, size = 16,face="bold"),
        legend.position = "none");p1



#2.2) col annotation 2
col_anno2 <- cor_colnames[,c("NAME_ABB","Category")] %>% 
  mutate(y=1,NAME_ABB=factor(NAME_ABB,level=rev(levels(ba_lg_mean$variable))),
         Category=factor(Category,levels=c("Unconjugated BA","Conjugated BA","Sum")))
col_anno2 <- col_anno2[order(col_anno2$NAME_ABB),]


p2 <- ggplot(col_anno2, aes(x=NAME_ABB, y=y)) + 
  geom_col(aes(fill=Category),width=1,color="white",linewidth=0.2)+#,color="white",size=0) +
  scale_fill_manual(values=col_BA_cat2)+
  labs(fill="Conjugation",y='Conjugation') +
  scale_y_continuous(breaks=c(0,1), limits=c(0,1),expand=c(0,0)) +
  coord_cartesian(expand = FALSE)+
  theme_bw()+
  theme(panel.border = element_rect(size=0.5), 
        axis.ticks = element_line(colour = "black", size = 0.8),
        axis.text = element_blank(),axis.ticks.y = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y.left = element_text(angle=0,vjust=0.5, size = 16,face="bold"),
        legend.position = "none");p2


#4) barplot mean abun
data=fmets_direct_BA_qc %>% filter(time==0) %>% select(BA_info_fmets$CHEM_ID2)
data[,BA_info_fmets$CHEM_ID2] <- map(data[,BA_info_fmets$CHEM_ID2],log10)
ba_lg_mean <- map(data,mean) %>% 
  data.frame() %>% setnames(BA_info_fmets$NAME_ABB) %>% 
  reshape2::melt(narm=F) %>% 
  left_join(BA_info_fmets[,c('NAME_ABB','SUB_PATHWAY2','Category')],
            by=c('variable'='NAME_ABB')); head(ba_lg_mean)

ba_rank=ba_lg_mean %>% arrange(desc(value)); 
ba_lg_mean$SUB_PATHWAY2 <- factor(ba_lg_mean$SUB_PATHWAY2,levels=rev(c('Primary BA','Secondary BA','Sum')))
ba_lg_mean$Category <- factor(ba_lg_mean$Category,levels=rev(c("Unconjugated BA","Conjugated BA","Sum")))
ba_lg_mean <- ba_lg_mean[order(ba_lg_mean$SUB_PATHWAY2,ba_lg_mean$Category,ba_lg_mean$value),]
ba_lg_mean$variable <- factor(ba_lg_mean$variable,levels=ba_lg_mean$variable)
ba_lg_mean <- ba_lg_mean[order(ba_lg_mean$SUB_PATHWAY2,ba_lg_mean$Category,ba_lg_mean$value),]
summary(ba_lg_mean$value)



col_BA_cat <- c('#6681D8D9',"#75D1C9","#ffb4a8")
col_BA_cat2 <- c( "#F7E987","#B5CDA3","#ffb4a8") 
p4 <- ggplot(ba_lg_mean, aes(x=variable, y=value)) + 
  geom_hline(yintercept = 0, linetype="dashed", size=1.2) +
  geom_col(aes(fill=SUB_PATHWAY2)) +
  scale_fill_manual(values=rev(col_BA_cat)) +
  scale_y_continuous(limits=c(0,10),breaks=c(0,2,4,6,8,10),
                     expand=c(0,0),position="right") +
  coord_flip(expand = FALSE) +
  labs(y=expression("Mean levels (log10)")) +
  theme_half_open()+
  theme(panel.border = element_rect(size=0.5), 
        panel.grid.major = element_blank(),
        axis.ticks  = element_line(colour = "black", size = 0.8),
        axis.text.x = element_text(size=13, colour = "black"),
        axis.title.x = element_text(hjust=0, vjust=-5,size=15, colour = "black"),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none");p4


pdf('results/Overview/Correlation/Cor_fmets_fmets_baseline.pdf',
    width=12,height=10,onefile=FALSE)
egg::ggarrange(p1,blank,
               blank,blank,
               p2,blank,
               blank,blank,
               p3, p4,
               nrow=5, ncol=2, 
               widths = c(1,0.2), heights = c(0.03,0.025,0.03,0.06,1))

dev.off()


#legend 1)-sub pathway
legend_cat  <- ggplot(col_anno, aes(x=NAME_ABB, y=y)) + 
  geom_col(aes(fill=SUB_PATHWAY2),width=0.99) +
  scale_fill_manual(values=col_BA_cat)+
  labs(fill="Sub Pathway") +
  theme(axis.text = element_blank(),axis.title=element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 17),
        legend.title = element_text(colour = "black", size = 17,face="bold"),
        legend.spacing = unit(0.1, 'cm'));legend_cat

pdf('results/Overview/Correlation/Cor_fmets_fmets_baseline_legend1.pdf',
    width=2.2,height=1.5,onefile=FALSE)
grid::grid.draw(ggpubr::get_legend(legend_cat))
dev.off() 


#legend 2.2)-category
legend_cat2  <- ggplot(col_anno2, aes(x=NAME_ABB, y=y)) + 
  geom_col(aes(fill=Category),width=0.99) +
  scale_fill_manual(values=col_BA_cat2)+
  labs(fill="Conjugation",y='Conjugation') +
  theme(axis.text = element_blank(),axis.title=element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 17),
        legend.title = element_text(colour = "black", size = 17,face="bold"),
        legend.spacing.x = unit(0.1, 'cm'),legend.spacing.y = unit(0.1, 'cm'));legend_cat2

pdf('results/Overview/Correlation/Cor_fmets_fmets_baseline_legend2.pdf',
    width=2.2,height=1.5,onefile=FALSE)
grid::grid.draw(ggpubr::get_legend(legend_cat2))
dev.off() 


#legend 2.3)--coefficient
legend_heat <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+ 
  geom_tile(color = "white")+ 
  scale_fill_gradient2(low = "#2166AC", high = "#B10C1F", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", na.value = "white",
                       name="Spearman \ncorrelation  ") +
  theme_minimal()+ # minimal theme 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.title = element_blank(), panel.grid = element_blank(),
        legend.position='right',legend.key.height=unit(0.7,'cm'),legend.direction = "vertical",
        legend.text = element_text(colour = "black", size = 17),
        legend.title = element_text(colour = "black", size = 17,face="bold"))+ 
  coord_fixed(ratio=1) ;legend_heat 
pdf('results/Overview/Correlation/Cor_fmets_fmets_baseline_legend_heatmap.pdf',
    width=1.5,height=3,onefile=FALSE)
grid::grid.draw(ggpubr::get_legend(legend_heat))
dev.off() 



#Panel D: pcoa-phylum & species-----
# top phylum
phylum <- direct_tax %>% 
  filter(!is.na(phylum)) %>% 
  select(-c("kingdom","class","order","family","genus","species")) %>% 
  column_to_rownames(var = "phylum") %>% t(); 
phylum <- phylum[,order(colSums(phylum)/nrow(phylum),decreasing=T)] %>% 
  as.data.frame() %>% 
  rownames_to_column(var='Sample_ID') %>% 
  mutate(Sample_ID=gsub('_taxonomic_profile','',Sample_ID)); 
head(colnames(phylum));#"p__Firmicutes","p__Bacteroidetes"


#data
df_spe <- metadata_spe_fmets0 %>% 
  select(all_of(c('sno',grep('s__.*_0$',colnames(metadata_spe_fmets0),value=T))) ) %>% 
  column_to_rownames(var='sno')
bc_dist_spe = vegdist(df_spe,"bray",na.rm=T) 
pcoa_spe <- cmdscale(bc_dist_spe,eig = T,k=2) 
dat_pcoa0 <- data.frame(pcoa_spe$points) %>% rownames_to_column(var = "sno")
(explainedvar1 <- round(pcoa_spe$eig[1]/sum(pcoa_spe$eig),4)*100)
(explainedvar2 <- round(pcoa_spe$eig[2]/sum(pcoa_spe$eig),4)*100)
(sum_eig <- sum(explainedvar1,explainedvar2)) 
dat_pcoa_phylum <- dat_pcoa0 %>% left_join(metadata[metadata$time==0,
                                                    c('sno','Sample_ID')],by='sno') %>% 
  inner_join(df_phylum,by='Sample_ID')


#top 2 phylum--
#Firmicutes
div_perm1 <- adonis2(bc_dist_spe~p__Firmicutes,dat_pcoa_phylum, 
                     permutations = 999,method = "bray"); div_perm1
(ado_txt <- paste0('PERMANOVA, ','P-value=',round(div_perm1$`Pr(>F)`[1],3)))#R2=0.042, P-value=4e-04
summary(dat_pcoa_phylum$p__Firmicutes);summary(dat_pcoa_phylum$p__Bacteroidetes)

pcoa_plot1 <-ggplot(dat_pcoa_phylum,aes(x=X1,y=X2))+
  geom_point(aes(fill=p__Firmicutes),size=4.7,shape=21,alpha=0.85)+#scale_fill_manual(name='',values=lvl2_color)+
  scale_fill_gradient2(name="",low ="#FFFFFF", mid = "#0DC7F3FF", high = "#1E276D",
                       midpoint=50,na.value="#0035A2FF")+
  labs(title='Firmicutes',
       x=paste0("PCo1 (", explainedvar1, "%)"),
       y=paste0("PCo2 (", explainedvar2, "%)"))+
  theme_few()+
  theme(panel.grid  = element_blank(),
        axis.text.y = element_text(size=20,angle=90),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0,size=22),
        legend.position='none',
        axis.title=element_text(size=21)); pcoa_plot1

pdf("results/Correlation/Fig0_spe_Firmicutes_PCOA.pdf", height = 5.1,width = 5)
pcoa_plot1
dev.off()


#2--
div_perm2 <- adonis2(bc_dist_spe~p__Bacteroidetes,dat_pcoa_phylum, 
                     permutations = 999,method = "bray"); div_perm2
(ado_txt <- paste0('PERMANOVA, ','P-value=',round(div_perm2$`Pr(>F)`[1],3)))#R2=0.042, P-value=4e-04

pcoa_plot2 <-ggplot(dat_pcoa_phylum,aes(x=X1,y=X2))+
  geom_point(aes(fill=p__Bacteroidetes),size=4.7,shape=21,alpha=0.9)+#scale_fill_manual(name='',values=lvl2_color)+
  scale_x_continuous(breaks=c(-0.25,0,0.25))+
  scale_fill_gradient2(name="",low ="#FFFFFF", mid = "#0DC7F3FF", high = "#1E276D",
                       midpoint=50,na.value="#0035A2FF")+
  labs(title='Bacteroidetes',
       x=paste0("PCo1 (", explainedvar1, "%)"),
       y=paste0("PCo2 (", explainedvar2, "%)"))+
  theme_few()+
  theme(panel.grid  = element_blank(),
        axis.text.y = element_text(size=20,angle=90),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0,size=22),
        legend.position='none',
        legend.key.height =unit(0.6,'cm'),
        legend.spacing.y = unit(0.7, "cm"),
        axis.title=element_text(size=21)); #pcoa_plot2

pdf("results/Overview/Fig0_spe_Bacteroidetes_PCOA.pdf", height = 5.1,width = 5)
pcoa_plot2
dev.off()


#legend--
pcoa_lg <-ggplot(dat_pcoa_phylum,aes(x=X1,y=X2))+
  geom_point(aes(fill=p__Firmicutes),size=4.5,shape=21,alpha=0.9)+#
  scale_fill_gradient2(name="",low ="#FFFFFF", mid = "#0DC7F3FF", high = "#1E276D",
                       midpoint=50,na.value="#0035A2FF")+
  theme_few()+
  theme(panel.grid  = element_blank(),
        axis.text.y = element_text(size=20,angle=90),
        axis.text.x = element_text(size=20),
        plot.title = element_text(hjust = 0,size=20),
        legend.position='bottom',legend.text=element_text(size=15),
        legend.key.width =unit(0.7,'cm'),
        axis.title=element_text(size=21)); pcoa_lg

pdf("results/Overview/Fig0_Firm_Bacte_PCOA_legend.pdf", height = 1,width = 2)
grid.draw(ggpubr::get_legend(pcoa_lg))
dev.off()


#S figure time trend of outcome ---------
df=metadata;
data=metadata %>% mutate(TC_HDL = Cholesterol/HDLc)
data0=copy(data) %>% filter(time==0); setnames(data0,lipid_var,paste0(lipid_var,'_0'))
data6=copy(data) %>% filter(time==6); setnames(data6,lipid_var,paste0(lipid_var,'_6'))
data18=copy(data) %>% filter(time==18); setnames(data18,lipid_var,paste0(lipid_var,'_18'))

dfch=data0[,c('sno','Group',paste0(lipid_var,'_0'))] %>% 
  left_join(data6[,c('sno',paste0(lipid_var,'_6'))],by='sno') %>% 
  left_join(data18[,c('sno',paste0(lipid_var,'_18'))],by='sno')

for (i in lipid_var){
  out0=paste0(i,'_0'); out6=paste0(i,'_6'); out18=paste0(i,'_18')
  out60=paste0(i,'_60');out180=paste0(i,'_180')
  dfch[,out60]=log(dfch[,out6],10)-log(dfch[,out0],10)
  dfch[,out180]=log(dfch[,out18],10)-log(dfch[,out0],10)
}

dfch2=bind_rows(dfch[,c('sno','Group',paste0(lipid_var,'_60'))] %>% mutate(time='60') %>% 
                  setnames(c('sno','Group',lipid_var,'Time')),
                dfch[,c('sno','Group',paste0(lipid_var,'_180'))] %>% mutate(time='180') 
                %>% setnames(c('sno','Group',lipid_var,'Time')) )

dbase=data.frame(dfch[,c('sno','Group',paste0(lipid_var,'_60'))] %>% mutate(time='00') %>% 
                   setnames(c('sno','Group',lipid_var,'Time')))
dbase[,lipid_var] <- apply(dbase[,lipid_var],2,function(x){x=0})
data=bind_rows(dfch2,dbase) %>% 
  mutate(Time=factor(Time,levels=c('00','60','180')),Group0='diet') %>% arrange(Time)
box <- ggplot(data=data,
              aes_string(x='Time',y='BMI',group='sno',color='Group',fill='Group'))+
  geom_point(size=1.5)+geom_line(size=1)+
  stat_summary(fun=mean, geom='point', size=3,aes_string(group='Group0'),color='#c03546') +
  stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group0'),color='#c03546') +
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
        plot.title = element_text(size=20)); box

#BMI
out='BMI'
data=metadata %>% mutate(TC_HDL = Cholesterol/HDLc,
                         Time=factor(case_when(time==0~'0mon',time==6~'6mon',time==18~'18mon'),
                                     levels=c('0mon','6mon','18mon')),
                         Group0='diet')
group_color <- c('#649BE7','#FF8911','#00A985')
p1=ggplot(data,aes(x=Time,y=BMI,group=sno,color=Group,fill=Group))+
  geom_point(size=2)+geom_line(alpha=0.9,size=1)+
  stat_summary(fun=mean, geom='point', alpha=1,size=3,aes_string(group='Group0'),color='#c03546') +
  stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group0'),color='#c03546') + 
  labs(y='',x='Time (Month)',title=bquote(paste(.(out),' (kg/m'^2~')')))+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size=15),
        plot.subtitle = element_text(hjust = 0,size=15));p1


#Triglycerides
out='Triglycerides'

p2=ggplot(data,aes(x=Time,y=Triglycerides,group=sno,color=Group,fill=Group))+
  geom_point(size=2)+geom_line(alpha=0.9,size=1)+
  stat_summary(fun=mean, geom='point', alpha=1,size=3,aes_string(group='Group0'),color='#c03546') +
  stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group0'),color='#c03546') + 
  labs(y='',x='Time (Month)',title=bquote(paste(.(out),' (mg/dL'~')')))+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size=15),
        plot.subtitle = element_text(hjust = 0,size=15));p2

#TC_HDL
out='TC_HDL'

p3=ggplot(data,aes(x=Time,y=TC_HDL,group=sno,color=Group,fill=Group))+
  geom_point(size=2)+geom_line(alpha=0.9,size=1)+
  stat_summary(fun=mean, geom='point', alpha=1,size=3,aes_string(group='Group0'),color='#c03546') +
  stat_summary(fun=mean, geom='line', size=1,aes_string(group='Group0'),color='#c03546') + 
  labs(y='',x='Time (Month)',title='TC/HDLc ratio')+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15),
        legend.position="none",
        plot.title = element_text(hjust = 0.5,size=15),
        plot.subtitle = element_text(hjust = 0,size=15));p3


pdf("results/Overview/lipid_var_trend_sample.pdf",width = 15, height = 5, onefile = F)
egg::ggarrange(p1,p2,p3,nrow=1)
dev.off()

#legend
p3_lg=ggplot(data,aes(x=Time,y=TC_HDL,group=sno,color=Group,fill=Group))+
  geom_point(size=2)+geom_line(alpha=0.9,linewidth=1)+
  scale_color_manual(values=group_color)+
  scale_fill_manual(values=group_color)+
  scale_x_discrete(limits = c('0mon', '6mon', '','18mon'), 
                   breaks = c('0mon', '6mon', '12mon','18mon'),
                   labels=c('0', '6', '','18'))+
  theme_void()+
  theme(legend.position="bottom",legend.title = element_text(size = 15),
        legend.text = element_text(size = 15,hjust=0),
        legend.key.spacing.x = unit(0.2,'cm'),
        plot.subtitle = element_text(hjust = 0,size=15));p3_lg

pdf("results/Overview/lipid_var_trend_sample_lg.pdf",width = 5, height =0.5, onefile = F)
grid.draw(ggpubr::get_legend(p3_lg))
dev.off()


