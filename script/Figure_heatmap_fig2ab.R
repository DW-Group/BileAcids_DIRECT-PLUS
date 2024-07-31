################################################################################################################################################
#
# Prospective Association of Baseline Levels of Fecal Bile Acids with Longitudinally Measured Body Adiposity and Lipid Biomarkers
#
##############################################################################################################################################

library(pacman)
pacman::p_load(tibble,ggpubr,grid,dplyr,psych,tidyverse)
pacman::p_load(ggplot2,vegan,ape,ggeffects,gee,geepack)
select <- dplyr::select

setwd(paste0(path,"/harvard/DIRECT_PLUS/Analysis/"))
load("/Results/DIRECT/dat_analysis.RData")

#data pre
gdata::keep(metadata, fmets_direct_BA, fmets_direct_BA_nm,fpmets_name,direct_tax,metadata_spe_fmets0,sure=T)

int <- function(x){
  qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))
}

geeglm_fun <- function(lipid_var, mets_id,df_input,form0,constr,family){
  res_fmets_out <- c()
  for (o in lipid_var){
    res_fmets <- c()
    for (m in mets_id){ 
      print(m)
      data <- df_input
      form <- formula(paste0(o,'~',m,form0)) 
      
      fit <- geeglm(form, id=ID, data=data,
                    corstr = corstr0,
                    family = family0);
      res <- summary(fit);
      res2 <- res$coefficients[2,] %>% 
        unlist()
      res_fmets <- rbind(res_fmets,c(o,gsub('\\_[0-9]$','',m),res2))
    }
    res_fmets <- res_fmets  %>% 
      as.data.frame() 
    res_fmets_out <- rbind(res_fmets_out,res_fmets) 
  }; 
  res_fmets_out <- res_fmets_out %>% 
    as.data.frame() %>% 
    setnames(c('lipid_var','F_mets','Estimate','S','Wald','P')) %>% 
    mutate(q_fdr=p.adjust(P,method='BH')) 
    left_join(BA_info_fmets[,c('CHEM_ID2','NAME','NAME_ABB', 'SUB_PATHWAY2','Category')],
              by=c('F_mets'='CHEM_ID2'))
  return(res_fmets_out)
}

geeglm_fun_T <- function(lipid_var, mets_id,df_input,form0,constr,family){
  res_fmets_out <- c()
  for (o in lipid_var){
    res_fmets <- c()
    for (m in mets_id){ #fecal_direct
      print(m)
      data <- df_input
      
      form <- formula(paste0(o,'~',m,form0)) 
      data[,m] <- cut(data[,m], 
                      breaks=c(tert[1]-10,tert[2],tert[3],tert[4]+10),
                      labels=c('T1','T2','T3'))
      
      fit <- geeglm(form, id=ID, data=data,
                    corstr = corstr0,
                    family = family0);
      (res <- summary(fit));
      res2 <- res$coefficients[2:3,] %>% 
        as.data.frame() %>% 
        mutate(F_mets=gsub('\\_[0-9]$','',m),
               lipid_var=o,
               var=c('T2 vs. T1','T3 vs. T1'))
      res_fmets <- rbind(res_fmets,res2)
    }
    res_fmets <- res_fmets  %>% as.data.frame() %>% 
      setnames(c('Estimate','S','Wald','P','F_mets','lipid_var','var'))
    res_fmets_out <- rbind(res_fmets_out,res_fmets) 
  }; 
  res_fmets_out <- res_fmets_out %>% 
    as.data.frame() %>% 
    left_join(BA_info_fmets[,c('CHEM_ID2','NAME','NAME_ABB', 'SUB_PATHWAY2','Category')],
              by=c('F_mets'='CHEM_ID2')) %>% 
    mutate(P=as.numeric(P),
           q_fdr=as.numeric(q_fdr)); 
  
  res_fmets_out <- res_fmets_out %>% mutate(CI_L=format(round(Estimate-1.96*SE,1),nsmall=1),
                                            CI_H=format(round(Estimate+1.96*SE,1),nsmall=1)) %>% 
    mutate(Beta_CI=paste0(format(round(Estimate,1),nsmall=1),' (',CI_L,', ',CI_H,')'))
  
  return(res_fmets_out)
}

dat_inp <- copy(fmets_direct_BA)
dat_inp[,mets_id] <- map(dat_inp[,mets_id], int)

lipid_var <- c("BMI","Triglycerides","TC_HDL","weight","WC","Cholesterol","LDLc","HDLc")

#MV1
form0 <- '+Group2+MED_idx+age+sex+antibio+metformin+Lipid_lowering+Time';#
corstr0 = constr; family0 = fam
gee_fmets_baseline <- geeglm_fun(lipid_var, mets_id, dat_inp,form0,constr0,family0); 
gee_fmets_baseline$NAME <- factor(gee_fmets_baseline$NAME,levels=BA_name_fmets_lvl)
#write.xlsx(gee_fmets_baseline, file = "results/GEE/GEE_lipid_var_baselineBA_MV1.xlsx") 

gee_fmets_baseline_T <- geeglm_fun_T(lipid_var, c(tb2_out_BA_id), dat_inp,form0,constr0,family0); 
gee_fmets_baseline_T$NAME <- factors(gee_fmets_baseline_T$NAME,levels=BA_name_fmets_lvl)
#write.xlsx(out,'results/GEE/Table_GEE_baselineBA_lipid_var_Tertile.xlsx')

#plot
df=gee_fmets_baseline %>%
  filter(lipid_var %in% c("BMI","Weight","Triglycerides","TC/HDLc","TC","LDLc","HDLc")) %>% 
  mutate(T_value=Estimate/S,
         lipid_var=factor(lipid_var,
                          levels=c("BMI","Weight","Triglycerides","TC/HDLc","HDLc","TC","LDLc")))

ht_est <- df %>% select(c('lipid_var','NAME_ABB','T_value')) %>% 
  reshape2::dcast(NAME_ABB~lipid_var,value.var='T_value') %>% 
  column_to_rownames(var='NAME_ABB');ht_est
ht_est[,1:ncol(ht_est)] <- map(ht_est[,1:ncol(ht_est)],as.numeric)
ht_q <- df %>% select(c('lipid_var','NAME_ABB','q_fdr')) %>% 
  reshape2::dcast(NAME_ABB~lipid_var,value.var='q_fdr') %>% 
  column_to_rownames(var='NAME_ABB')

q_anno <- ht_q
if (!is.null(q_anno)){
  sssmt <-  q_anno < 0.001
  q_anno[sssmt] <-'****'
  ssmt <-  
    q_anno >= 0.001 & q_anno < 0.05
  q_anno[ssmt] <-'***'
  smt <-  
    q_anno >= 0.05 & q_anno <0.15
  q_anno[smt] <- '**'
  mt <-  
    q_anno >= 0.15 & q_anno <0.25
  q_anno[mt] <- '*'
  q_anno[!sssmt &!ssmt & !smt & !mt]<- '' #
} else {
  q_anno <- F
};q_anno

#row annotation
(row_name <- data.frame('NAME_ABB'=rownames(ht_est)) %>% 
    left_join(BA_info_fmets,by='NAME_ABB') ) %>% 
  mutate(Category=factor(Category,levels=BA_category_lvl))
annotation_row =  data.frame('Sub_pathway'=row_name$SUB_PATHWAY2,
                             'Conjugation'=factor(row_name$Category,
                                                  levels=BA_category_lvl)) 
groupcolor <- c("#8788DA","#37D3B2",'#FA86B0'); 
names(groupcolor) <- c("Primary BA", "Secondary BA",'Sum')
conju=c("#FFD14D",'#B5CDA3','#FA86B0');
names(conju)=BA_category_lvl
ann_colors=list('Sub_pathway'=groupcolor,
                'Conjugation'=conju)
summary(ht_est)

ht <- pheatmap(as.matrix(ht_est),
               cellheight=15,
               cellwidth=30,
               angle_col = '315',
               fontsize_row = 15,fontsize_col=15,
               treeheight_row = 25,
               show_rownames = T,
               show_colnames = T,
               
               border_color='white',
               cluster_rows=T,
               cluster_cols=F,
               annotation_row = annotation_row,
               annotation_colors=ann_colors,
               color= colorRamp2(c(-2,0,2),c("#006796","white","#E64B35")),
               
               row_names_side = "right", name='Beta/SD',
               display_numbers = as.matrix(q_anno), 
               fontsize_number=20,
               number_color='white',
               legend=T,
               heatmap_legend_param = list(legend_height=unit(2.5,"cm"),
                                           legend_direction='vertical',
                                           title_position="topleft",
                                           labels_gp = gpar(fontsize = 10),#legend label
                                           title_gp = gpar(fontsize = 10,
                                                           fontface = "bold"))
); ht

pdf('results/GEE/All_lipids_BAs_heatmap.pdf',height=15,width=10)
ht
dev.off()

#scatter plot with line -----
# data pre
data <- fmets_direct_BA %>% filter(time==0) %>% 
  select(all_of(c('CLIENT_SAMPLE_ID2',mets_id))) %>% 
  left_join(metadata[,c('sno','time','Time',lipid_var)],
            by=c('CLIENT_SAMPLE_ID2'='sno')) 

setnames(data,'CLIENT_SAMPLE_ID2','sno')
data[,mets_id] <- map_df(data[,mets_id],LOG10)
setnames(data,mets_id,BA_info_fmets$NAME_ABB)

gee_fmets_baseline %>% filter(as.numeric(q_fdr<0.05))
df <- gee_fmets_baseline %>% 
  mutate(q_fdr=p_3dg_anno(q_fdr))

lipid_var=c('BMI','Triglycerides','TC_HDL');
lipid_var_nm=c('BMI','Triglycerides','TC/HDL ratio')

# plot-BMI
lipid_id=lipid_var[1]; 
lipid_nm=lipid_var_nm[1]
show_list=c("LCA","3HCOA")

p0 <- c()
for (n in c(1:length(show_list))){
  ba=show_list[n]; 
  out=lipid_id;
  out_nm=lipid_nm
  qd=df %>% filter(NAME_ABB %in% ba & lipid_var %in% out)
  
  if (n==1){
      p=ggplot(data,
             aes_string(x=ba,y=out))+
      geom_point(size=1,color='#4870b7')+
      geom_smooth(method=lm, size=1,alpha=0.3,
                  fill='#F98F34',color='#F98F34')+
      labs(x=ba,y='',title=bquote(paste(.(lipid_nm),~'(kg/m'^.(2),')')),
           subtitle=bquote(paste(italic('q')~.(qd$q_fdr))))+
      theme_classic()+
      theme(axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            legend.position='none',
            plot.title=element_text(size=15),
            plot.subtitle=element_text(size=15),
            axis.text.y=element_text(angle=90,size=15),
            axis.text.x=element_text(size=15));p
    p0[[ba]] <- p
  }


  if (n==2){
    p=ggplot(tem %>% filter(!!sym(ba)>=qq2[1] & !!sym(ba)<=qq2[2]),
             aes_string(x=ba,y=out))+
      geom_point(size=1,color='#4870b7')+
      geom_smooth(method=lm, size=1,alpha=0.3,
                  fill='#F98F34',color='#F98F34')+
      labs(x=ba,y='',subtitle=bquote(paste(italic('q')~.(qd$q_fdr))))+
      theme_classic()+
      theme(axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            legend.position='none',
            plot.subtitle=element_text(size=15),
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=15));p
    p0[[ba]] <- p
  }
}



# plot-TG
lipid_id=lipid_var[2]; 
lipid_nm=lipid_var_nm[2]
show_list=c('UDCA','UCA')

p1 <- c()
for (n in c(1:length(show_list))){
  ba=show_list[n]; 
  out=lipid_id;
  out_nm=lipid_nm
  qd=df %>% filter(NAME_ABB %in% ba & lipid_var %in% out)
  
  if (n==1){
    p=ggplot(data,
             aes_string(x=ba,y=out))+
      geom_point(size=1,color='#4870b7')+
      geom_smooth(method=lm, size=1,alpha=0.3,
                  fill='#F98F34',color='#F98F34')+
      labs(x=ba,y='',title=bquote(paste(.(lipid_nm),~'(mg/dL)')),
           subtitle=bquote(paste(italic('q')~.(qd$q_fdr))))+
      theme_classic()+
      theme(axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            legend.position='none',
            plot.title=element_text(size=15),
            plot.subtitle=element_text(size=15),
            axis.text.y=element_text(angle=90,size=15),
            axis.text.x=element_text(size=15));p
    p1[[ba]] <- p
  }
  
  
  if (n==2){
    p=ggplot(tem %>% filter(!!sym(ba)>=qq2[1] & !!sym(ba)<=qq2[2]),
             aes_string(x=ba,y=out))+
      geom_point(size=1,color='#4870b7')+
      geom_smooth(method=lm, size=1,alpha=0.3,
                  fill='#F98F34',color='#F98F34')+
      labs(x=ba,y='',subtitle=bquote(paste(italic('q')~.(qd$q_fdr))))+
      theme_classic()+
      theme(axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            legend.position='none',
            plot.subtitle=element_text(size=15),
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=15));p
    p1[[ba]] <- p
  }
}


# plot-TC/HDLc ratio
lipid_id=lipid_var[3]; 
lipid_nm=lipid_var_nm[3]
show_list=c('UDCA','3HCOA')

p2 <- c()
for (n in c(1:length(show_list))){
  ba=show_list[n]; 
  out=lipid_id;
  out_nm=lipid_nm
  qd=df %>% filter(NAME_ABB %in% ba & lipid_var %in% out)
  if (n==1){
    p=ggplot(data,
             aes_string(x=ba,y=out))+
      geom_point(size=1,color='#4870b7')+
      geom_smooth(method=lm, size=1,alpha=0.3,
                  fill='#F98F34',color='#F98F34')+
      labs(x=ba,y='',title=bquote(paste(.(lipid_nm),~'(mg/dL)')),
           subtitle=bquote(paste(italic('q')~.(qd$q_fdr))))+
      theme_classic()+
      theme(axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            legend.position='none',
            plot.title=element_text(size=15),
            plot.subtitle=element_text(size=15),
            axis.text.y=element_text(angle=90,size=15),
            axis.text.x=element_text(size=15));p
    p2[[ba]] <- p
  }
  
  
  if (n==2){
    p=ggplot(tem %>% filter(!!sym(ba)>=qq2[1] & !!sym(ba)<=qq2[2]),
             aes_string(x=ba,y=out))+
      geom_point(size=1,color='#4870b7')+
      geom_smooth(method=lm, size=1,alpha=0.3,fill='#F98F34',color='#F98F34')+
      labs(x=ba,y='',subtitle=bquote(paste(italic('q')~.(qd$q_fdr))))+
      theme_classic()+
      theme(axis.title.x=element_text(size=15),
            axis.title.y=element_blank(),
            legend.position='none',
            plot.subtitle=element_text(size=15),
            axis.text.y=element_blank(),
            axis.text.x=element_text(size=15));p
    p2[[ba]] <- p
  }
}


pdf('results/GEE/Line_scatter_BAs_lipids.pdf',width=10,height=15,onefile=F)
egg::ggarrange(p1,p2,p3,nrow=3,heights=c(1,1,1))

dev.off()


