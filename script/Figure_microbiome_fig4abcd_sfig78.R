################################################################################################################################################
#
# The improvement in body adiposity and lipid profiles induced by the Mediterranean diet intervention 
#   varied based on baseline levels of fecal bile acids.
#
##############################################################################################################################################


pacman::p_load(tibble,ggpubr,grid,dplyr,psych,tidyverse,ggplot2)
select <- dplyr::select

setwd(paste0(path,"/harvard/DIRECT_PLUS/Analysis/"))
load("/Results/DIRECT/dat_analysis.RData")

#data pre
gdata::keep(metadata, dat_baseline_BA, BA_info_fmets,taxa,metadata_spe_fmets0,ec_con_spe_l,baec_show,inter_spe_SBAs,sure=T)

#species--dissimilarity matrix-bray------
metadata_spe_fmets_fmets0 <- metadata_spe_fmets %>% 
  filter(time==0) %>% 
  left_join(dat_baseline_BA[,c('sno',"Time",mets_id)],by=c('sno','Time'))

df_spe <- metadata_spe_fmets_fmets0 %>% 
  select(all_of(c('sno',grep('^s__',colnames(metadata_spe_fmets_fmets0),value=T))) ) %>% 
  column_to_rownames(var='sno')
bc_dist_spe = vegdist(df_spe,"bray",na.rm=T) 
pcoa_spe <- cmdscale(bc_dist_spe,eig = T,k=2) 

(explainedvar1 <- round(pcoa_spe$eig[1]/sum(pcoa_spe$eig),4)*100) 
(explainedvar2 <- round(pcoa_spe$eig[2]/sum(pcoa_spe$eig),4)*100)
(sum_eig <- sum(explainedvar1,explainedvar2))

dat_pcoa_spe <- data.frame(pcoa_spe$points) %>% 
  rownames_to_column(var = "sno") %>% 
  left_join(metadata_spe_fmets_fmets0,by="sno")


#metabolites dissimilarity matrix-euclidean-------
df_ba <- metadata_spe_fmets_fmets0 %>% 
  select(all_of(c('sno',BA_info_fmets$BA_FID)) ) %>% 
  column_to_rownames(var='sno')

(explainedvar1_ba <- round(pcoa_ba$eig[1]/sum(pcoa_ba$eig),4)*100) 
(explainedvar2_ba <- round(pcoa_ba$eig[2]/sum(pcoa_ba$eig),4)*100)
(sum_eig_ba <- sum(explainedvar1_ba,explainedvar2_ba)) 


df_sba <- metadata_spe_fmets_fmets0 %>% 
  select(all_of(c('sno',BA_info_fmets[BA_info_fmets$SUB_PATHWAY2 %in% 'Secondary BA',]$BA_FID)) ) %>% 
  column_to_rownames(var='sno')

bc_dist_sba = vegdist(df_sba,"euclidean",na.rm=T) 
pcoa_sba <- cmdscale(bc_dist_sba,eig = T,k=2) 
dat_pcoa_sba <- data.frame(pcoa_sba$points) %>% 
  setnames(c('D1_sba','D2_sba')) %>% 
  rownames_to_column(var = "sno")# 

#PCoA according to the PCo1 of SBAs----
(explainedvar1_sba <- round(pcoa_sba$eig[1]/sum(pcoa_sba$eig),4)*100) # 29.44%
(explainedvar2_sba <- round(pcoa_sba$eig[2]/sum(pcoa_sba$eig),4)*100) 
(sum_eig_sba <- sum(explainedvar1_sba,explainedvar2_sba))  
div_perm1 <- adonis2(bc_dist_spe~D1_sba,dat_pcoa_sba, permutations = 999,method = "bray"); div_perm1

dat_pcoa_both <- dat_pcoa_sba %>% 
  left_join(dat_pcoa_spe[,c('sno','X1','X2')],by='sno') %>% 
  left_join(dat_pcoa_both[,c('sno','FPBA','FSBA')],by='sno') 
pcoa_plot1 <- ggplot(dat_pcoa_both,aes(x=X1,y=X2))+
  geom_point(aes(color=D1_sba),size=5,alpha=0.9)+
  scale_color_gradient2(name="PCo1 ofSBAs",midpoint = 0,
                        low = "#A2FFEB", mid = "#adbcf1", high = "#1134af")+ 
  labs(title='', 
       x=paste0("PCo1 (", explainedvar1, "%)"),
       y=paste0("PCo2 (", explainedvar2, "%)"))+
  theme_classic()+
  annotate('text',x=0.25,y=0.3,size=5,
           label=bquote(paste(R^2~'=',
                              .(format(round(div_perm1$'R2'[1],3)*100,nsmall=1)),'%')))+
  annotate('text',x=0.25,y=0.2,size=5,
           label=bquote(paste(italic(p)~'=',.(round(div_perm1$`Pr(>F)`[1],3)))))+
  theme(plot.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20,angle=0)); pcoa_plot1

pdf('results/Function/PCOA_PCOA_species_PCOAsecBA.pdf',height=5,width=5)
pcoa_plot1
dev.off()

#legend--
pcoa_plot1_lg <- ggplot(dat_pcoa_both,aes(x=X1,y=X2))+
  geom_point(aes(color=D1_sba),size=5,alpha=0.9)+
  scale_color_gradient2(name="PCo1 ofSBAs",midpoint = 0,
                        low = "#A2FFEB", mid = "#adbcf1", high = "#1134af",)+ #004533
  theme_few()+
  theme(panel.grid  = element_blank(),
        plot.title = element_text(hjust = 0,size=20),
        legend.text=element_text(size=20),
        legend.position = 'right',
        axis.title=element_text(size=20,family='Arial'),
        axis.text.x = element_text(size=20),axis.text.y = element_text(size=20,angle=90)); 
pdf('results/Function/PCOA_PCOA_species_PCOAsecBA_legend.pdf',width=1,height=2)
grid::grid.draw(ggpubr::get_legend(pcoa_plot1_lg))
dev.off()


#according to groups of BAs
color_4g=rev(c("#4870b7","#A0CBE8","#FABFD2","#C6566E"))
dat_pcoa_both$PBA_2g <- ifelse(dat_pcoa_both$FPBA>median(dat_pcoa_both$FPBA,na.rm=T),
                               'High PBAs','Low PBAs'); 
table(dat_pcoa_both$PBA_2g)
dat_pcoa_both$SBA_2g <- ifelse(dat_pcoa_both$FSBA>median(dat_pcoa_both$FSBA,na.rm=T),
                               'High SBAs','Low SBAs');
table(dat_pcoa_both$SBA_2g)
dat_pcoa_both$PSBA_4g <- paste0(dat_pcoa_both$PBA_2g,' + ',dat_pcoa_both$SBA_2g); 
table(dat_pcoa_both$PSBA_4g)
dat_pcoa_both$PSBA_4g <- factor(dat_pcoa_both$PSBA_4g,
                                levels=c('Low PBAs + Low SBAs','High PBAs + Low SBAs',
                                         'Low PBAs + High SBAs','High PBAs + High SBAs'))

dat_pcoa_both$group <- factor(case_when(dat_pcoa_both$PSBA_4g=='Low PBAs + Low SBAs'~'a',
                                        dat_pcoa_both$PSBA_4g=='High PBAs + Low SBAs'~'b',
                                        dat_pcoa_both$PSBA_4g=='Low PBAs + High SBAs'~'c',
                                        dat_pcoa_both$PSBA_4g=='High PBAs + High SBAs'~'d',
                                        TRUE~as.character(dat_pcoa_both$group)),
                              levels=c('a','b','c','d'))


p1 <- ggplot(dat_pcoa_both,aes(x=X1,y=X2))+
  geom_point(aes(color=PSBA_4g),size=5,alpha=0.9)+
  scale_color_manual(values=color_4g)+
  labs(y=paste0("PCo2 (", explainedvar2, "%)"),x='')+
  theme_classic()+
  annotate('text',x=0.25,y=0.3,size=5,
           label=bquote(paste(R^2~'=',.(format(round(div_perm2$'R2'[1],3)*100,nsmall=1)),'%')))+
  annotate('text',x=0.25,y=0.2,size=5,
           label=bquote(paste(italic(p)~'=',.(round(div_perm2$`Pr(>F)`[1],3)))))+
  theme(plot.title = element_text(hjust = 0,size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        legend.position = 'none', 
        axis.title.y=element_text(size=20),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),axis.text.y=element_text(size=20,angle=0)); p1

p3 <- ggplot(dat_pcoa_both,aes(x=PSBA_4g,y=X2,fill=PSBA_4g))+
  geom_boxplot(alpha=0.9)+
  scale_fill_manual(values=color_4g)+
  scale_color_manual(values=color_4g)+
  labs(y='')+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        axis.ticks=element_blank(),
        legend.position = "none",
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_blank()); p3

p2 <- ggplot(dat_pcoa_both,aes(x=PSBA_4g,y=X1,fill=PSBA_4g))+
  geom_boxplot(alpha=0.9)+
  scale_fill_manual(values=color_4g)+
  scale_color_manual(values=color_4g)+
  scale_x_discrete(expand=c(0,0))+
  labs(y=paste0("PCo1 (", explainedvar1, "%)"))+
  theme_void()+
  coord_flip()+
  theme(plot.title = element_text(hjust = 0.5,size=20),
        legend.position = "none", 
        axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x = element_text(size=20),axis.text.y=element_text(size=20,angle=0)); p2

lay_pcoa <- lay_new(
  mat = matrix(1:4, ncol = 2),widths = c(5,1),heights = c(5,1))  
lay_show(lay_pcoa)
pdf('results/Function/Fig_PCOA_PSBA_4g.pdf',width=8,height=8)
plots2 = lapply(c(1:4), function(x) get(paste0("p", x)))
lay_grid(plots2, lay_pcoa)
dev.off()


#legend 
p_leg <- ggplot(dat_pcoa_both,aes(x=X1,y=X2))+
  geom_point(aes(color=Group),size=5,alpha=0.9)+
  scale_color_manual(values=rev(color_4g))+
  theme_void()+
  theme(legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.position = 'right', 
        axis.title=element_blank(),
        axis.text = element_text(size=20)); p_leg

pdf('results/Function/Fig_PCOA_PSBA_4g_legend.pdf',width=3,height=2)
grid::grid.draw(ggpubr::get_legend(p_leg))
dev.off()


# heatmap plot showing association of SBAs with microbiome
#Heatmap_lm-------
lm_spe_sba <- read.xlsx('results/Function/Mass_lm_spe_SBAs.xlsx')
lm_spe_sba %>% left_join(taxa[,c('species','genus','phylum')],by='species')
(df_coef <- reshape2::dcast(lm_spe_sba, species~NAME_ABB,value.var = 'Estimate') %>% 
    mutate(species=factor(species,levels=spe_tar)) %>% arrange(species) %>% 
    column_to_rownames(var='species') %>% as.data.frame())
(df_p <- reshape2::dcast(lm_spe_sba, species~NAME_ABB,value.var = 'q_fdr') %>% 
    mutate(species=factor(species,levels=spe_tar)) %>% arrange(species) %>% 
    column_to_rownames(var='species') %>% as.data.frame())
df_coef[,1:ncol(df_coef)] <- map_df(df_coef[,1:ncol(df_coef)],as.numeric)
df_p[,1:ncol(df_p)] <- map_df(df_p[,1:ncol(df_p)],as.numeric)

p_anno <- df_p

if (!is.null(p_anno)){
  ssmt <- p_anno <0.05
  p_anno[ssmt]<- '*'
  p_anno[!ssmt] <- ''
} else {
  p_anno <- F
};p_anno

tem1 <- BA_info_fmets %>% filter(NAME_ABB %in% colnames(df_coef)) %>% 
  mutate(NAME_ABB=factor(NAME_ABB,levels=colnames(df_coef)))
tem1 <- tem1[order(tem1$NAME_ABB),]
(ann_col <- list('Sub pathway'=tem1$SUB_PATHWAY2,'Category'=tem1$Category) %>% 
    as.data.frame() %>% 
    setnames(c('Sub pathway','Category')))
ann_col$Category <- factor(ann_col$Category,levels=BA_category_lvl)

names(col_BA_cat) <- c("Primary BA","Secondary BA","Sum")
names(col_BA_cat2) <- c('Unconjugated BA','Conjugated BA','Sum')

dheat=df_coef;summary(dheat)
ht <- pheatmap(as.matrix(dheat),
               name='Coefficient',
               cellheight=20,cellwidth=20,angle_col = '45',
               fontsize_row = 15,fontsize_col=15,
               show_rownames = T,show_colnames = T,border_color='gray70',
               cluster_rows=F,
               cluster_cols=T,
               annotation_colors=ann_colors,
               color= col,
               row_names_side = "left", 
               display_numbers = as.matrix(p_anno), 
               fontface_row = 'italic', 
               legend=F,
               heatmap_legend_param = list(legend_height=unit(3,"cm"),
                                           legend_direction='vertical',
                                           title_position="topleft",
                                           labels_gp = gpar(fontsize = 15),
                                           title_gp = gpar(fontsize = 15,
                                                           fontface = "bold"))
); ht

pdf('results/Function/Heatmap_spe_list_SBA.pdf',width=10,height=5)
ht
dev.off() 

#legend--
ht_legend <- pheatmap(as.matrix(dheat),
                      name='Coefficient',
                      cellheight=20,cellwidth=20,angle_col = '45',
                      fontsize_row = 15,fontsize_col=15,
                      show_rownames = T,show_colnames = T,border_color='gray70',
                      cluster_rows=F,
                      cluster_cols=T,
                      annotation_colors=ann_colors,
                      color= col,
                      row_names_side = "left", 
                      display_numbers = as.matrix(p_anno), 
                      fontface_row = 'italic', 
                      legend=T,
                      heatmap_legend_param = list(legend_height=unit(3,"cm"),
                                                  legend_direction='vertical',
                                                  title_position="topleft",
                                                  labels_gp = gpar(fontsize = 15),
                                                  title_gp = gpar(fontsize = 15,
                                                                  fontface = "bold"))
); ht_legend
pdf('results/Function/Heatmap_spe_list_SBA_legend.pdf',width=5,height=5)
ht_legend
dev.off()


#side bar
ec_spe_review <- read.xlsx('results/Function/EC_BA/EC_species_review.xlsx') %>% 
  setnames(c('species','3α-HSDH','3β-HSDH','7α-HSDH','7β-HSDH','7α-dehydroxylase','BSH'))

df_ecspe_bar <- reshape2::melt(ec_spe_review,id='species') %>% 
  setnames(c('species','enzyme','value')) %>% mutate(value=as.numeric(value)) %>% 
  left_join(ec_con_spe_l,by=c('species','enzyme')) %>% 
  filter(value==1) %>% 
  mutate(species=factor(species,levels=rev(spe_lvl))) %>% 
  arrange(species)

col_ec <- c("#8DD3C7", "#FFFFB3","#83A2FF","#FF8888","#78D4F8","#FDB462","#B3DE69")
p_ec <- ggplot(df_ecspe_bar, aes(x = species, y = value)) +
  geom_col(aes(fill=enzyme), position = position_stack(reverse = T)) +
  scale_fill_manual(values=col_ec)+
  scale_x_discrete(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  coord_flip() +
  theme_few() +
  labs(title = "", y = "Enzyme", x = "",fill='Enzyme') +
  theme(axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_text(size=20,face='bold'),
        legend.position='none'); p_ec 
pdf('results/Function/Bar_ec_spe.pdf',width=2,height=5)
p_ec
dev.off()


#legend
p_ec_lg <- ggplot(df_ecspe_bar, aes(x = species, y = value)) +
  geom_col(aes(fill=enzyme), position = position_stack(reverse = TRUE)) +
  scale_fill_manual(values=col_ec)+
  scale_x_discrete(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
  coord_flip() +
  theme_few() +
  theme(axis.text = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position='right',
        legend.title = element_text(size=20,face='bold'),
        legend.text = element_text(size=20)); p_ec_lg
pdf('results/Function/Bar_ec_spe_legend.pdf',width=5,height=5)
grid::grid.draw(ggpubr::get_legend(p_ec_lg))
dev.off()


#Panel C--enzyme--
baec_p <- read.xlsx('results/Function/EC_BAs_FDR.xlsx')
baec_show <- baec_show %>% filter(ec_nm %in% c('7β-HSDH','3β-HSDH','7α-dehydroxylase')) %>% 
  left_join(baec_p,by=c('NAME_ABB'='Bile.acids'))

col_ec_2g <- c('#4f74c2','#E56C68')
plots <- c()
for (i in 1){
  ba <- baec_show$ba_id[i]; ba_nm <- baec_show$NAME_ABB[i];
  ec_g <- baec_show$ec_g[i]; ec_nm <- baec_show$ec_nm[i]
  q <- p_23dg_snno(as.numeric(baec_show$'q.(FDR)'[i])); 
  max_v <- max(df_ec_ba[,ba],na.rm = T)
  p1 <- ggplot(df_ec_ba,
               aes_string(y=ba, x=ec_g,fill=ec_g,color=ec_g)) +
    geom_boxplot(alpha=0.85,size=0.6)+
    geom_jitter(alpha=0.85,size=2)+ #
    scale_color_manual(values=col_ec_2g)+
    scale_fill_manual(values=col_ec_2g)+
    labs(y=paste0(ba_nm," (log10)"), x=ec_nm) +
    annotate('text',x=1.5,y=max_v,
             label=bquote(paste(italic('*q')~.(q))),size=6)+
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 15,angle=90,hjust=0.5),
          axis.title.x = element_text(size = 15, color = "#2267C1"),
          axis.title.y = element_text(size = 15),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text( hjust = 0.5)); p1 
  
  plots[[i]] <- p1
  pdf(paste0('results/Function/EC_BA/Plots_baec_',ba_nm,'_',ec_g,'.pdf'),width=5,height = 5)
  print(p1)
  dev.off()
}

for (i in 3){
  ba <- baec_show$ba_id[i]; ba_nm <- baec_show$NAME_ABB[i];
  ec_g <- baec_show$ec_g[i]; ec_nm <- baec_show$ec_nm[i]
  q <- p_23dg_snno(as.numeric(baec_show$'q.(FDR)'[i])); 
  max_v <- max(df_ec_ba[,ba],na.rm = T)
  p1 <- ggplot(data=df_ec_ba,
               aes_string(y=ba, x=ec_g,fill=ec_g,color=ec_g)) +
    geom_boxplot(alpha=0.85,size=0.6)+
    geom_jitter(alpha=0.85,size=2)+ 
    scale_color_manual(values=col_ec_2g)+
    scale_fill_manual(values=col_ec_2g)+
    labs(y=paste0(ba_nm," (log10)"), x=ec_nm) +
    annotate('text',x=1.5,y=max_v,
             label=bquote(paste(italic('*q')~.(q))),size=6)+
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 15,angle=90,hjust=0.5),
          axis.title.x = element_text(size = 15, color = "#2267C1"),
          axis.title.y = element_text(size = 15),
          
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text( hjust = 0.5)); p1 
  
  plots[[i]] <- p1
  pdf(paste0('results/Function/EC_BA/Plots_baec_',ba_nm,'_',ec_g,'.pdf'),width=5,height = 5)
  print(p1)
  dev.off()
}

for (i in 2){
  ba <- baec_show$ba_id[i]; ba_nm <- baec_show$NAME_ABB[i];
  ec_g <- baec_show$ec_g[i]; ec_nm <- baec_show$ec_nm[i]
  q <- p_23dg_snno(as.numeric(baec_show$'q.(FDR)'[i])); 
  max_v <- max(df_ec_ba[,ba],na.rm = T)
  
  p1 <- ggplot(data=df_ec_ba,
               aes_string(y=ba, x=ec_g,fill=ec_g,color=ec_g)) +
    geom_boxplot(alpha=0.85,size=0.6)+
    geom_jitter(alpha=0.85,size=2)+ #
    scale_color_manual(values=col_ec_2g)+
    scale_fill_manual(values=col_ec_2g)+
    labs(y=paste0(ba_nm," (log10)"), x=ec_nm) +
    annotate('text',x=1.5,y=max_v,
             label=bquote(paste(italic('*q')~.(q))),size=6)+
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 15,angle=90,hjust=0.5),
          axis.title.x = element_text(size = 15, color = "#2267C1"),
          axis.title.y = element_text(size = 15),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text( hjust = 0.5)); p1 
  
  plots[[i]] <- p1
  pdf(paste0('results/Function/EC_BA/Plots_baec_',ba_nm,'_',ec_g,'.pdf'),width=5,height = 5)
  print(p1)
  dev.off()
}

for (i in 5){
  ba <- baec_show$ba_id[i]; ba_nm <- baec_show$NAME_ABB[i];
  ec_g <- baec_show$ec_g[i]; ec_nm <- baec_show$ec_nm[i]
  q <- p_23dg_snno(as.numeric(baec_show$'q.(FDR)'[i])); 
  max_v <- max(df_ec_ba[,ba],na.rm = T)
  p1 <- ggplot(data=df_ec_ba,
               aes_string(y=ba, x=ec_g,fill=ec_g,color=ec_g)) +
    geom_boxplot(alpha=0.85,size=0.5)+
    geom_jitter(alpha=0.85,size=2)+ #
    scale_color_manual(values=col_ec_2g)+
    scale_fill_manual(values=col_ec_2g)+
    labs(y=paste0(ba_nm," (log10)"), x='7α-dehydroxylase') +
    annotate('text',x=1.5,y=min_v+0.5,
             label=bquote(paste(italic('*q')~.(q))),size=6)+
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 15,angle=90,hjust=0.5),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text( hjust = 0.5)); p1 
  plots[[i]] <- p1
  pdf(paste0('results/Function/EC_BA/Plots_baec_',ba_nm,'_',ec_g,'.pdf'),width=5,height = 5)
  print(p1)
  dev.off()
}
for (i in 4){
  ba <- baec_show$ba_id[i]; ba_nm <- baec_show$NAME_ABB[i];
  ec_g <- baec_show$ec_g[i]; ec_nm <- baec_show$ec_nm[i]
  q <- p_23dg_snno(as.numeric(baec_show$'q.(FDR)'[i])); 
  max_v <- max(df_ec_ba[,ba],na.rm = T)
  p1 <- ggplot(data=df_ec_ba,
               aes_string(y=ba, x=ec_g,fill=ec_g,color=ec_g)) +
    geom_boxplot(alpha=0.85,size=0.5)+
    geom_jitter(alpha=0.85,size=2)+ 
    scale_color_manual(values=col_ec_2g)+
    scale_fill_manual(values=col_ec_2g)+
    labs(y=paste0(ba_nm," (log10)"), x='7α-dehydroxylase') +
    annotate('text',x=1.5,y=max_v,
             label=bquote(paste(italic('*q')~'',.(q))),size=6)+
    theme_classic() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),axis.ticks.x = element_blank(), 
          axis.text.y = element_text(size = 15,angle=90,hjust=0.5),
          axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          legend.position = "none",
          legend.title = element_blank(),
          legend.text = element_text(size = 15),
          plot.title = element_text( hjust = 0.5)); p1 
  plots[[i]] <- p1
  pdf(paste0('results/Function/EC_BA/Plots_baec_',ba_nm,'_',ec_g,'.pdf'),width=5,height = 5)
  print(p1)
  dev.off()
}


#legend
i=1
ba <- baec_show$ba_id[i]; ba_nm <- baec_show$NAME_ABB[i];
ec_g <- baec_show$ec_g[i]; ec_nm <- baec_show$ec_nm[i]
q <- p_23dg_snno(as.numeric(baec_show$'q.(FDR)'[i])); 
max_v <- max(df_ec_ba[,ba],na.rm = T)
p_lg <- ggplot(data=df_ec_ba,
               aes_string(y=ba, x=ec_g,fill=ec_g,color=ec_g)) +
  geom_boxplot(width=0.6,alpha=0.85)+
  scale_color_manual(values=col_ec_2g, name='Enzyme')+
  scale_fill_manual(values=col_ec_2g, name='Enzyme')+
  theme_void() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        plot.title = element_text( hjust = 0.5)); p_lg

pdf('results/Function/EC_BA/Plots_baec_legend.pdf',width=4,height=2)
grid::grid.draw(ggpubr::get_legend(p_lg))
dev.off()



#Panel D--errorbar plot-----
metadata_spe_fmets_fmets0 <- metadata_spe_fmets %>% 
  filter(time==0) %>% 
  left_join(fmets_direct_BA,by=c('sno'="CLIENT_SAMPLE_ID2",'time'))

spe_show <- rep(c("s__Ruminococcus_torques",'s__Bifidobacterium_longum'),3) 
spe_show_nm <- rep(c("R. torques",'B. longum'),3) 
ba_show <- c('3-HCOA','LCA','TDCA','7,12diketo-LCA','HCA','isoUDCA')
lipid_var_show=c('BMI','BMI','Triglycerides','Triglycerides','TC_HDL','TC_HDL')

show_spe <- cbind(spe_show,ba_show) %>% as.data.frame() %>% setnames(c('species','NAME_ABB')) %>% 
  mutate(lipid_var=lipid_var_show)
res_fmets_str <- inter_spe_SBAs %>% 
  inner_join(show_spe,by=c('species','NAME_ABB','lipid_var')) %>% 
  mutate(q_fdr_2g=format(round(q_fdr_2g,3),nsmall=3)) %>% 
  select(all_of(c('species','NAME_ABB','lipid_var','value','Estimate','CI_L','CI_H','p_inter'))) %>% 
  as.data.frame() 
#plot- str by spe [carriers/noncarriers]-BA--------
plots2 <- list()
n1=(1:length(spe_show))[1:length(spe_show) %% 2 == 1]
for (i in n1){
  spe <- spe_show[i]; m <- BA_info_fmets[BA_info_fmets$NAME_ABB %in% ba_show[i],]$BA_FID
  o=lipid_var_show[i];m_nm <- ba_show[i]
  data=res_fmets_str %>% filter(species==spe & BA_FID==m &lipid_var==o)
  p_inter <- unique(p_23dg_anno(data$p_inter))
  p1 <- ggplot(data=data,aes(x=Estimate,y=species,color=value))+
    geom_vline(xintercept = 0,linetype='dashed',color='gray60')+
    geom_point(size=5)+
    geom_errorbarh(aes(xmax=CI_H,xmin=CI_L),height=0,size=1)+
    labs(title=bquote(paste(.(m_nm))),
         x='Beta (95%CI)',y=spe_show_nm[i])+
    scale_color_manual(values=c('#EB6A4E', "#419e9b"))+#
    annotate('text',x=sum(data$Estimate)/2,y=1,
             label=bquote(italic('p')[inter]~.(p_inter)),size=8)+
    theme_classic()+
    theme(axis.text.y = element_blank(), 
          axis.text.x = element_text(size=20),
          axis.title.x= element_text(size=20),
          axis.title.y= element_text(size=20,face='italic'),
          legend.position = 'none',
          legend.title = element_blank(),legend.text = element_text(size=20),
          plot.title = element_text(hjust = 0,size=20),
          plot.subtitle = element_text(hjust = 0,size=0));p1 #
  
  plots2[[i]] <- p1
  
}

n2=(1:length(spe_show))[1:length(spe_show) %% 2 == 0]
for (i in n2){
  spe <- spe_show[i]; 
  m <- BA_info_fmets[BA_info_fmets$NAME_ABB %in% ba_show[i],]$BA_FID
  o=lipid_var_show[i];
  m_nm <- ba_show[i]
  data=res_fmets_str %>% filter(species==spe & BA_FID==m &lipid_var==o)
  p_inter <- unique(p_23dg_anno(data$p_inter))
  
  p1 <- ggplot(data=data,aes(x=Estimate,y=species,color=value))+
    geom_vline(xintercept = 0,linetype='dashed',color='gray60')+
    geom_point(size=5)+
    geom_errorbarh(aes(xmax=CI_H,xmin=CI_L),height=0,size=1)+
    labs(title=bquote(paste(.(m_nm))),
         x='Beta (95%CI)',y=spe_show_nm[i])+
    scale_color_manual(values=c('#EB6A4E', "#419e9b"))+
    annotate('text',x=sum(data$Estimate)/2,y=1,
             label=bquote(italic('p')[inter]~.(p_inter)),size=8)+
    theme_classic()+
    theme(axis.text.y = element_blank(), axis.text.x = element_text(size=20),
          axis.title.x= element_text(size=20),
          axis.title.y= element_text(size=20,face='italic'),
          legend.position = 'none',
          legend.title = element_blank(),
          legend.text = element_text(size=20),
          plot.title = element_text(hjust = 0,size=20),
          plot.subtitle = element_text(hjust = 0,size=0));p1 #
  
  plots2[[i]] <- p1
  
}

blank <- ggplot()+theme_void(); plots2[[11]] <- blank; 

lay_pcoa <- lay_new(
  mat = matrix(1:6, ncol = 3),widths = c(1,1,1),heights = c(1,1) )  
pdf('results/Function/EC_BA/Plots_ba_str_species.pdf',width=15,height = 8)
lay_grid(plots2,lay_pcoa) 
dev.off()

#legend-1--
for (i in 1){
  spe <- spe_show[i]; m <- BA_info_fmets[BA_info_fmets$NAME_ABB %in% ba_show[i],]$BA_FID
  o=lipid_var_show[i];m_nm <- ba_show[i]
  data=res_fmets_str %>% filter(species==spe & BA_FID==m &lipid_var==o)
  p1_lg <- ggplot(data=data,aes(x=Estimate,y=species,color=value))+
    geom_vline(xintercept = 0,linetype='dashed',color='gray60')+
    geom_point(size=5)+
    geom_errorbarh(aes(xmax=CI_H,xmin=CI_L),height=0,size=1)+
    scale_color_manual(name=spe_show_nm[i],values=c('#EB6A4E', "#419e9b"))+
    labs(fill=spe_show_nm[i])+
    theme_void()+
    theme(legend.position = 'right',legend.key.size = unit(20,'pt'),
          legend.title = element_text(size=20,face='bold.italic'),
          legend.text = element_text(size=20));p1_lg #
  pdf('results/Function/EC_BA/Plots_ba_str_species_legend1.pdf',width=4,height=2)#5.5
  grid::grid.draw(ggpubr::get_legend(p1_lg))
  dev.off()
}

#legend-2--
for (i in 2){
  spe <- spe_show[i]; m <- BA_info_fmets[BA_info_fmets$NAME_ABB %in% ba_show[i],]$BA_FID
  o=lipid_var_show[i];m_nm <- ba_show[i]
  data=res_fmets_str %>% filter(species==spe & BA_FID==m &lipid_var==o)
  p1_lg <- ggplot(data=data,aes(x=Estimate,y=species,color=value))+
    geom_vline(xintercept = 0,linetype='dashed',color='gray60')+
    geom_point(size=5)+
    geom_errorbarh(aes(xmax=CI_H,xmin=CI_L),height=0,size=1)+
    scale_color_manual(name=spe_show_nm[i],values=c('#EB6A4E', "#419e9b"))+
    labs(fill=spe_show_nm[i])+
    theme_void()+
    theme(legend.position = 'right',legend.key.size = unit(20,'pt'),
          legend.title = element_text(size=20,face='bold.italic'),
          legend.text = element_text(size=20));p1_lg 
  pdf('results/Function/EC_BA/Plots_ba_str_species_legend2.pdf',width=4,height=2)
  grid::grid.draw(ggpubr::get_legend(p1_lg))
  dev.off()
}


#supplemental fig 8------
lipid_var_show=c(rep('BMI',4),rep('Triglycerides',4),rep('TC_HDL',4))
spe_show=c('HCA: s__Eggerthella_lenta','UCA: s__Bifidobacterium_bifidum',
          'GUDCA: s__Ruminococcus_gnavus','UDCA: s__Bacteroides_plebeius',
          'LCA-S1: s__Eggerthella_lenta','TDCA: s__Eggerthella_lenta',
          'LCA-S1: s__Ruminococcus_gnavus','DCA: s__Clostridium_sp_CAG_167',
          '3-DHDCA: s__Bifidobacterium_longum','7,12diketo-LCA: s__Bifidobacterium_longum',
          'DCA: s__Clostridium_sp_CAG_167','GLCA: s__Clostridium_spiroforme')

show_sup <- data.frame('lipid_var'=lipid_var_show,'spe_x'=spe_show)

df_plt <- inter_spe_SBAs %>% inner_join(show_sup,by=c('lipid_var','spe_x'))

plots <- c()
for (i in 1:length(spe_show)){
  data=df_plt%>% filter(spe_x %in% spe_show[i] & lipid_var %in% lipid_var_show[i])
  spe <- unique(data$species); spe_show_nm <- gsub('_',' ',gsub('s__','',spe))
  o=lipid_var_show[i]; 
  m <- unique(data$BA_FID); m_nm <- unique(data$NAME_ABB)
  p_inter <- unique(p_23dg_anno(data$p_inter))
  p1 <- ggplot(data=data,aes(x=Estimate,y=species,color=value))+
    geom_vline(xintercept = 0,linetype='dashed',color='gray60')+
    geom_point(size=5)+
    geom_errorbarh(aes(xmax=CI_H,xmin=CI_L),height=0,size=1)+
    labs(title=bquote(paste(.(m_nm))),
         subtitle=bquote(paste('(Group: ',italic(.(spe_show_nm))~')')),
         x='Beta (95% CI)',y='')+
    scale_color_manual(values=c("#EB6A4E", "#419e9b"))+#
    annotate('text',x=sum(data$Estimate)/2,y=1,
             label=bquote(italic('p')[inter]~''~.(p_inter)),size=8)+
    theme_classic()+
    theme(axis.text.y = element_text(size=20,angle=90),
          axis.text.x = element_text(size=15),
          axis.title.x= element_text(size=20),axis.title.y= element_text(size=15),
          legend.position = 'none',
          legend.title = element_blank(),legend.text = element_text(size=15),
          plot.title = element_text(hjust = 0,size=20),
          plot.subtitle = element_text(hjust = 0,size=20));p1 #
  
  plots[[i]] <- p1
  
}


lay_pcoa <- lay_new(
  mat = matrix(1:12, ncol = 3),widths = c(1,1,1),heights = c(1,1,1,1) )  
pdf('results/Function/Fig_supple_ba_str_species.pdf',width=20,height = 15)
lay_grid(plots,lay_pcoa) #
dev.off()

#legend--
for (i in length(spe_show)){
  data=tem %>% filter(spe_x %in% spe_show[i] & lipid_var %in% lipid_var_show[i])
  spe <- unique(data$species); spe_show_nm <- gsub('_',' ',gsub('s__','',spe))
  o=lipid_var_show[i]; m <- unique(data$BA_FID); m_nm <- unique(data$NAME_ABB)
  p_sup_lg <- ggplot(data=data,aes(x=Estimate,y=species,color=value))+
    geom_vline(xintercept = 0,linetype='dashed',color='gray60')+
    geom_point(size=5)+
    geom_errorbarh(aes(xmax=CI_H,xmin=CI_L),height=0,size=1)+
    scale_color_manual(name='Microbial species',values=c('#EB6A4E', "#419e9b"))+
    theme_void()+
    theme(legend.position = 'bottom',
          legend.title = element_text(size=15,face='bold'),
          legend.direction = 'horizontal',
          legend.text = element_text(size=15));p_sup_lg 
  
  pdf('results/Function/Plots_ba_str_species_legend_supple.pdf',width=5,height=2)#5.5
  grid::grid.draw(ggpubr::get_legend(p_sup_lg))
  dev.off()
  
}


#supplemental fig 7------
#other species-ba res---
dat_meta_spe_fmets0 <- fmets_direct_BA %>% 
  filter(time==0) %>% as.data.frame() %>% 
  select(-all_of(spe_tar)); dim(dat_meta_spe_fmets0)
head(colnames(dat_meta_spe_fmets0),10)

spe_al=grep('^s__',colnames(dat_meta_spe_fmets0),value=T)
ba_id=c(BA_info_fmets[BA_info_fmets$SUB_PATHWAY2=='Secondary BA',]$CHEM_ID2)
lm_spe_sba_sup <- read.xlsx('results/Function/Mass_lm_spe_SBAs_sup.xlsx')
res_final <- lm_spe_sba_sup %>% left_join(BA_info_fmets[,c('CHEM_ID2','NAME_ABB')],by='CHEM_ID2') 
res_filter=unique(res_final[res_final$q<0.05,]$species); 
data <- res_final %>% filter(species %in% res_filter) %>% mutate(species=gsub('_',' ',gsub('s__','',species)))
res_spe_ba_other=res_final %>% filter(species %in% res_filter) 
(df_coef <- reshape2::dcast(data, species~NAME_ABB,value.var = 'Beta/SD') %>% 
    column_to_rownames(var='species') %>% as.data.frame()); dim(df_coef)
(df_p <- reshape2::dcast(data, species~NAME_ABB,value.var = 'q') %>% 
    column_to_rownames(var='species') %>% as.data.frame())
p_anno <- df_p

if (!is.null(p_anno)){
  sssmt <- p_anno <0.001 
  p_anno[sssmt]<- "***"
  ssmt <- p_anno >=0.001 & p_anno <0.01
  p_anno[ssmt]<- '**'
  smt <- p_anno >=0.01 & p_anno <0.05 
  p_anno[smt]<- '*'
  p_anno[!ssmt &!sssmt & !smt] <- ''
} else {
  p_anno <- F
};p_anno

tem1 <- BA_info_fmets %>% filter(NAME_ABB %in% colnames(df_coef)) %>% mutate(NAME_ABB=factor(NAME_ABB,levels=colnames(df_coef)))
tem1 <- tem1[order(tem1$NAME_ABB),]
(ann_col <- list('Sub_pathway'=tem1$SUB_PATHWAY2,'Category'=tem1$Category) %>% as.data.frame())
ann_col$Category <- factor(ann_col$Category,levels=BA_category_lvl)

#taxa
taxa=taxa %>% mutate(spe=gsub('_',' ',gsub('^s__','',species)))
tem_row=data.frame(spe=rownames(df_coef)) %>% left_join(taxa,by='spe') %>% 
  mutate(Genus=gsub('_',' ',gsub('^g__','',genus)),
         Phylum=gsub('_',' ',gsub('^p__','',phylum))) %>% 
  mutate(spe=factor(spe,levels=rownames(df_coef))) %>% arrange(spe)
ann_row=data.frame('Phylum'=tem_row$Phylum)
col_phy=rev(c('#A0ACBD','#70DAF4','#BF83A5',"#F7D42A","#6a60a9",'#5AAE61'))
names(col_phy)=unique(anno_row$Phylum)

names(col_BA_cat) <- c("Primary BA","Secondary BA","Sum")
names(col_BA_cat2) <- c('Unconjugated BA','Conjugated BA','Sum')
(ann_colors=list('Category'=col_BA_cat2[1:3], 
                 'Phylum'=col_phy))

dheat=df_coef;summary(dheat)
col=colorRamp2(c(-4,0,4),c("#1C63A1","white","#CB2C3E")) 
ht <- ComplexHeatmap::pheatmap(as.matrix(dheat), name='Beta/SD', 
                               cellheight=15,cellwidth=30,angle_col = '315',
                               fontsize_row = 15,fontsize_col=15,
                               show_rownames = T,show_colnames = T,border_color='white',#scale='none',
                               cluster_rows=T,
                               cluster_cols=T,
                               annotation_row = ann_row,
                               annotation_col=ann_col, annotation_colors=ann_colors,
                               color= col,
                               row_names_side = "right", 
                               display_numbers = as.matrix(p_anno), fontsize_number=15,number_color='white',
                               legend=T,fontface_row = 'italic',
                               heatmap_legend_param = list(legend_height=unit(3,"cm"),#legend高度
                                                           legend_direction='vertical',
                                                           title_position="topleft",
                                                           labels_gp = gpar(fontsize = 12),#legend label字体
                                                           title_gp = gpar(fontsize = 12,
                                                                           fontface = "bold"))
); ht


pdf('results/Function/Heatmap_spe_allBAs_linear.pdf',width=30,height=20)
ht
dev.off()








