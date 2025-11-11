# date: 11/11/2025
# author:  Perera-Cruz O.

# abbreviated code for cMD characterization
# analysis between regionally-averaged surface brain data and independent variables (age, cognition, plasma biomarkers)
# only contains main parts of the code, for clarity

# library

library(dplyr)
library(readr)
library(lessR)
library(table1)
library(foreach)
library(tidyverse)
library(psych)
library(ggpubr)
library(ggplot2)
library(corrplot)
library(ppcor)
library(WGCNA)
library(survival)
library(ggpmisc)
library(broom)
library(AICcmodavg)
library(lme4)
library(car)
library(patchwork)
library(jtools)
library(lme4)
library(geepack)
library(emmeans)
library(car)
library(ggVennDiagram)

# independent correlations

pcor.test(df_SA_BBHI$Age, df_SA_BBHI$RAVLT_total_raw, df_SA_BBHI[,c("Educ_yrs","Sex")], method="spearman")
pcor.test(df_SA_BBHI$Age, df_SA_BBHI$RAVLT_delayed_raw, df_SA_BBHI[,c("Educ_yrs","Sex")], method="spearman")
pcor.test(NoNA_df$PTAU, NoNA_df$NFL, NoNA_df[,c("Age","Sex")], method="spearman")
pcor.test(df_PTAU$Age, df_PTAU$PTAU, df_PTAU[,c("Sex")], method="spearman")
pcor.test(df_NFL$Age, df_NFL$NFL, df_NFL[,c("Sex")], method="spearman")
pcor.test(df_SA_BBHI$RAVLT_delayed_raw, df_SA_BBHI$RAVLT_total_raw, df_SA_BBHI[,c("Educ_yrs","Sex","Age")], method="spearman")
pcor.test(df_PTAU$PTAU, df_PTAU$RAVLT_total_raw, df_PTAU[,c("Sex","Educ_yrs","Age")], method="spearman")
pcor.test(df_NFL$RAVLT_delayed_raw, df_NFL$NFL, df_NFL[,c("Sex","Educ_yrs","Age")], method="spearman")
pcor.test(df_PTAU$PTAU, df_PTAU$RAVLT_delayed_raw, df_PTAU[,c("Sex","Educ_yrs","Age")], method="spearman")
pcor.test(df_NFL$RAVLT_total_raw, df_NFL$NFL, df_NFL[,c("Sex","Educ_yrs","Age")], method="spearman")
pcor.test(df_NFL$RAVLT_total_raw, df_NFL$NFL, df_NFL[,c("Sex","Educ_yrs","Age")], method="spearman")

# correlation analyses and Venn diagram on all-bmk sub-sample

all_bmk<-subset(database, !is.na(Ptau_BC) & !is.na(Log10_hsCRP) & !is.na(Log10_NFL) & !is.na(ApoE4)) %>%
  filter(Study=="BBHI")

pcor.test(all_bmk$Real.Age, all_bmk$PACC, all_bmk[,c("Educ_yrs","Sex")], method="pearson")
pcor.test(all_bmk$Real.Age, all_bmk$Ptau_BC, all_bmk[,c("Sex")], method="pearson")
pcor.test(subset(all_bmk, !is.na(Log10_NFL))$Real.Age,
          subset(all_bmk, !is.na(Log10_NFL))$Log10_NFL,
          subset(all_bmk, !is.na(Log10_NFL))[,c("Sex")], method="pearson")
pcor.test(subset(all_bmk, !is.na(Log10_hsCRP))$Real.Age,
          subset(all_bmk, !is.na(Log10_hsCRP))$Log10_hsCRP,
          subset(all_bmk, !is.na(Log10_hsCRP))[,c("Sex")], method="pearson")
pcor.test(subset(all_bmk, !is.na(Log10_hsCRP))$PACC,
          subset(all_bmk, !is.na(Log10_hsCRP))$Log10_hsCRP,
          subset(all_bmk, !is.na(Log10_hsCRP))[,c("Sex","Real.Age","Educ_yrs")], method="pearson")
pcor.test(subset(all_bmk, !is.na(Log10_NFL))$PACC,
          subset(all_bmk, !is.na(Log10_NFL))$Log10_NFL,
          subset(all_bmk, !is.na(Log10_NFL))[,c("Sex","Real.Age","Educ_yrs")], method="pearson")
pcor.test(subset(all_bmk, !is.na(Ptau_BC))$PACC,
          subset(all_bmk, !is.na(Ptau_BC))$Ptau_BC,
          subset(all_bmk, !is.na(Ptau_BC))[,c("Sex","Real.Age","Educ_yrs")], method="pearson")
pcor.test(all_bmk$Ptau_BC, all_bmk$Log10_NFL, all_bmk[,c("Sex","Real.Age")])
pcor.test(all_bmk$Ptau_BC, all_bmk$Log10_hsCRP, all_bmk[,c("Sex","Real.Age")])
pcor.test(all_bmk$Log10_NFL, all_bmk$Log10_hsCRP, all_bmk[,c("Sex","Real.Age")])

wide_bmk<-bmk_data %>%
  mutate(value=1) %>%
  pivot_wider(
    names_from=bmk, values_from=value, values_fill = 0
  )

p<-ggVennDiagram(list(
  pTau181 = wide_bmk$id[wide_bmk$Ptau == 1],
  NfL = wide_bmk$id[wide_bmk$NFL == 1],
  hsCRP = wide_bmk$id[wide_bmk$hsCRP == 1],
  APOE4 = wide_bmk$id[wide_bmk$APOE4 ==1]
)) + scale_fill_gradient(low = "white", high = "skyblue")

p$layers[[3]]$aes_params$size <- 6  # region numbers
p$layers[[4]]$aes_params$size <- 6  # category labels

# risk profile and age group comparison

bmk_data <- database %>%
  mutate(Ptau=ifelse(Study=="SA",NA,Ptau)) %>%
  pivot_longer(
    cols = c("Ptau", "NFL", "hsCRP","APOE4"),
    names_to = "bmk",
    values_to = "bmk_value"
  ) %>%
  select(id, Real.Age,bmk, bmk_value,MMSE,PACC,Sex, Educ_yrs) %>%
  filter(!is.na(bmk_value)) %>%
  mutate(Sex=as.numeric(ifelse(Sex==1,0,1)))

summary(aov(Real.Age~bmk, bmk_data))                        # F + pval
lm(PACC~bmk+Real.Age+Sex+Educ_yrs, bmk_data)                # estimate (b) + pval
glm(Sex~bmk, bmk_data, family=binomial)                     # estimate (z) + pval
Anova(glm(ApoE4~Age_group+Sex, database, family=binomial), type=3)      # Chisq + pval
Anova(lm(Log10_hsCRP~Age_group+Sex, database), type=3)      # F + pval

# cMD and CTh vs Age plotting

lh_plots<-list()
rh_plots<-list()
bilat_plots<-list()
new_df<-df_all[,c("id","Sex","Real.Age","ApoE4","NFL_group","hsCRP_g2","Ptau_75","Study")]
new_df$Sex<-as.character(new_df$Sex)
new_df$ApoE4<-as.character(new_df$ApoE4)
new_df$NFL_group<-as.character(new_df$NFL_group)
new_df$hsCRP_g2<-as.character(new_df$hsCRP_g2)
new_df$Ptau_75<-as.character(new_df$Ptau_75)
new_df_ptau<-filter(new_df, Study=="BBHI")
ROIs<-list(c("lh.medialorbitofrontal","lh_medialorbitofrontal_thickness",
             "rh.medialorbitofrontal","rh_medialorbitofrontal_thickness"),
           c("lh.middletemporal","lh_middletemporal_thickness",
             "rh.middletemporal","rh_middletemporal_thickness"),
           c("lh.posteriorcingulate","lh_posteriorcingulate_thickness",
             "rh.posteriorcingulate","rh_posteriorcingulate_thickness"),
           c("lh.rostralanteriorcingulate","lh_rostralanteriorcingulate_thickness",
             "rh.rostralanteriorcingulate","rh_rostralanteriorcingulate_thickness"),
           c("lh.pericalcarine","lh_pericalcarine_thickness",
             "rh.pericalcarine","rh_pericalcarine_thickness"),
           c("lh.precentral","lh_precentral_thickness",
             "rh.precentral","rh_precentral_thickness"))
groups<-c("ApoE4","NFL_group","hsCRP_g2","Ptau_75")

for (var in ROIs){
  cMD_ROI_L<-var[1]
  CTh_ROI_L<-var[2]
  cMD_ROI_R<-var[3]
  CTh_ROI_R<-var[4]
  ROI_name<-sub("lh.","", var[1])
  ROI_name <- paste0(toupper(substr(ROI_name, 1, 1)), substr(ROI_name, 2, nchar(ROI_name)))
  
  for (grp in groups){
    if (grp=="Ptau_75"){
      
      new_df_ptau[[paste0(cMD_ROI_L,"_std")]]<-as.numeric(scale(subset(df_all, Study=="BBHI")[[cMD_ROI_L]]))
      new_df_ptau[[paste0(CTh_ROI_L,"_std")]]<-as.numeric(scale(subset(df_all, Study=="BBHI")[[CTh_ROI_L]]))
      new_df_ptau[[paste0(cMD_ROI_R,"_std")]]<-as.numeric(scale(subset(df_all, Study=="BBHI")[[cMD_ROI_R]]))
      new_df_ptau[[paste0(CTh_ROI_R,"_std")]]<-as.numeric(scale(subset(df_all, Study=="BBHI")[[CTh_ROI_R]]))
      
      new_df_ptau[[paste0(ROI_name,".bilat")]]<-as.numeric((scale(subset(df_all, Study=="BBHI")[[cMD_ROI_L]])
                                                            +scale(subset(df_all, Study=="BBHI")[[cMD_ROI_R]]))/2)
      new_df_ptau[[paste0(ROI_name,"_bilat")]]<-as.numeric((scale(subset(df_all, Study=="BBHI")[[CTh_ROI_L]])
                                                            +scale(subset(df_all, Study=="BBHI")[[CTh_ROI_R]]))/2)
      
      cMD_lm_lh<-lm(new_df_ptau[[paste0(cMD_ROI_L,"_std")]]~new_df_ptau$Sex)
      CTh_lm_lh<-lm(new_df_ptau[[paste0(CTh_ROI_L,"_std")]]~new_df_ptau$Sex)
      cMD_lm_rh<-lm(new_df_ptau[[paste0(cMD_ROI_R,"_std")]]~new_df_ptau$Sex)
      CTh_lm_rh<-lm(new_df_ptau[[paste0(CTh_ROI_R,"_std")]]~new_df_ptau$Sex)
      cMD_lm_bilat<-lm(new_df_ptau[[paste0(ROI_name,".bilat")]]~new_df_ptau$Sex)
      CTh_lm_bilat<-lm(new_df_ptau[[paste0(ROI_name,"_bilat")]]~new_df_ptau$Sex)
      
      new_df_ptau[[paste0(cMD_ROI_L,"_res")]]<-as.numeric(resid(cMD_lm_lh))
      new_df_ptau[[paste0(CTh_ROI_L,"_res")]]<-as.numeric(resid(CTh_lm_lh))
      new_df_ptau[[paste0(cMD_ROI_R,"_res")]]<-as.numeric(resid(cMD_lm_rh))
      new_df_ptau[[paste0(CTh_ROI_R,"_res")]]<-as.numeric(resid(CTh_lm_rh))
      new_df_ptau[[paste0(ROI_name,".bilat_res")]]<-as.numeric(resid(cMD_lm_bilat))
      new_df_ptau[[paste0(ROI_name,"_bilat_res")]]<-as.numeric(resid(CTh_lm_bilat))
      
    }else{
      
      new_df[[paste0(cMD_ROI_L,"_std")]]<-as.numeric(scale(df_all[[cMD_ROI_L]]))
      new_df[[paste0(CTh_ROI_L,"_std")]]<-as.numeric(scale(df_all[[CTh_ROI_L]]))
      new_df[[paste0(cMD_ROI_R,"_std")]]<-as.numeric(scale(df_all[[cMD_ROI_R]]))
      new_df[[paste0(CTh_ROI_R,"_std")]]<-as.numeric(scale(df_all[[CTh_ROI_R]]))
      
      new_df[[paste0(ROI_name,".bilat")]]<-as.numeric((scale(df_all[[cMD_ROI_L]])
                                                       +scale(df_all[[cMD_ROI_R]]))/2)
      new_df[[paste0(ROI_name,"_bilat")]]<-as.numeric((scale(df_all[[CTh_ROI_L]])
                                                       +scale(df_all[[CTh_ROI_R]]))/2)
      
      cMD_lm_lh<-lm(new_df[[paste0(cMD_ROI_L,"_std")]]~new_df$Sex)
      CTh_lm_lh<-lm(new_df[[paste0(CTh_ROI_L,"_std")]]~new_df$Sex)
      cMD_lm_rh<-lm(new_df[[paste0(cMD_ROI_R,"_std")]]~new_df$Sex)
      CTh_lm_rh<-lm(new_df[[paste0(CTh_ROI_R,"_std")]]~new_df$Sex)
      cMD_lm_bilat<-lm(new_df[[paste0(ROI_name,".bilat")]]~new_df$Sex)
      CTh_lm_bilat<-lm(new_df[[paste0(ROI_name,"_bilat")]]~new_df$Sex)
      
      new_df[[paste0(cMD_ROI_L,"_res")]]<-as.numeric(resid(cMD_lm_lh))
      new_df[[paste0(CTh_ROI_L,"_res")]]<-as.numeric(resid(CTh_lm_lh))
      new_df[[paste0(cMD_ROI_R,"_res")]]<-as.numeric(resid(cMD_lm_rh))
      new_df[[paste0(CTh_ROI_R,"_res")]]<-as.numeric(resid(CTh_lm_rh))
      new_df[[paste0(ROI_name,".bilat_res")]]<-as.numeric(resid(cMD_lm_bilat))
      new_df[[paste0(ROI_name,"_bilat_res")]]<-as.numeric(resid(CTh_lm_bilat))
      
    }
    
    formula_cMD<-paste0(var[1],"~Real.Age*factor(",grp,")+Sex")
    model_cMD<-lm(as.formula(formula_cMD), df_all)
    coef_text_cMD <- paste0("Slope", ifelse(summary(model_cMD)$coefficients[5,1]<0.001," < 0.001",
                                  paste0(" = ", round(summary(model_cMD)$coefficients[5,1], 3))),      # slope for Age:grp
                        "\nAdj R-squared", ifelse(summary(model_cMD)$adj.r.squared<0.001, " < 0.001",
                                  paste0(" = ", (round(summary(model_cMD)$adj.r.squared,3)))),         # global model adj-R
                        "\np-value", ifelse(summary(model_cMD)$coefficients[5,4]<0.05,"<0.05",
                                  paste0(" = ", round(summary(model_cMD)$coefficients[5,4], 3)))       # p-val for Age:grp
                        )
    formula_CTh<-paste0(var[2],"~Real.Age*factor(",grp,")+Sex")
    model_CTh<-lm(as.formula(formula_CTh), df_all)
    coef_text_CTh <- paste0("Slope", ifelse(summary(model_CTh)$coefficients[5,1]<0.001," < 0.001",
                                            paste0(" = ", round(summary(model_CTh)$coefficients[5,1], 3))),      # slope for Age:grp
                            "\nAdj R-squared", ifelse(summary(model_CTh)$adj.r.squared<0.001, " < 0.001",
                                                      paste0(" = ", (round(summary(model_CTh)$adj.r.squared,3)))),         # global model adj-R
                            "\np-value", ifelse(summary(model_CTh)$coefficients[5,4]<0.05,"<0.05",
                                                paste0(" = ", round(summary(model_CTh)$coefficients[5,4], 3)))       # p-val for Age:grp
    )

    plot_lh<-ggplot(data=subset(new_df, !is.na(new_df[[grp]])), aes(x=Real.Age)) +
      geom_point(aes(y=.data[[paste0(cMD_ROI_L,"_res")]], color="cMD"), alpha=0.2) +
      geom_smooth(aes(y=.data[[paste0(cMD_ROI_L,"_res")]], color="cMD", linetype=!!sym(grp)), method=loess) +
      #stat_cor(aes(y=.data[[paste0(cMD_ROI_L,"_res")]], color="cMD"), size=5, label.y=3) +
      #stat_poly_eq(aes(y=.data[[paste0(cMD_ROI_L,"_std")]],
       #               label = paste(after_stat(adj.rr.label),
        #                            after_stat(p.value.label), sep="~~~"),
      #color="cMD"), formula=formula, size=5, label.y=0.15, label.x=0.71, parse=T) +
      geom_point(aes(y=.data[[paste0(CTh_ROI_L,"_res")]], color="CTh"), alpha=0.2) +
      geom_smooth(aes(y=.data[[paste0(CTh_ROI_L,"_res")]], color="CTh", linetype=!!sym(grp)), method=loess) +
      #stat_cor(aes(y=.data[[paste0(CTh_ROI_L,"_res")]], color="CTh"), size=5, label.y=3.7) +
      theme_minimal() + scale_color_manual(name="Brain metric",
                                           values=c("cMD"="#f18f01","CTh"="#00296b")) +
      labs(
        #title=paste("cMD and CTh across age for", ROI_name, "cortex (LH)"),
         #  subtitle=paste("Standardized variables corrected for Sex and grouped by", grp),
           x="Age (yrs)", y="Brain metric (z-scores)",  shape=grp, linetype=grp) + ylim(-4,5) +
         annotate("text", x=(ifelse(grp=="ApoE4",51,47)), y=-3, label=coef_text_cMD, size=5, color="#f18f01") +
         annotate("text", x=(ifelse(grp=="ApoE4",73,62)), y=-3, label=coef_text_CTh, size=5, color="#00296b")

    formula_cMD<-paste0(var[3],"~Real.Age*factor(",grp,")+Sex")
    model_cMD<-lm(as.formula(formula_cMD), df_all)
    coef_text_cMD <- paste0("Slope", ifelse(summary(model_cMD)$coefficients[5,1]<0.001," < 0.001",
                                            paste0(" = ", round(summary(model_cMD)$coefficients[5,1], 3))),      # slope for Age:grp
                            "\nAdj R-squared", ifelse(summary(model_cMD)$adj.r.squared<0.001, " < 0.001",
                                                      paste0(" = ", (round(summary(model_cMD)$adj.r.squared,3)))),         # global model adj-R
                            "\np-value", ifelse(summary(model_cMD)$coefficients[5,4]<0.05,"<0.05",
                                                paste0(" = ", round(summary(model_cMD)$coefficients[5,4], 3)))       # p-val for Age:grp
    )
    formula_CTh<-paste0(var[4],"~Real.Age*factor(",grp,")+Sex")
    model_CTh<-lm(as.formula(formula_CTh), df_all)
    coef_text_CTh <- paste0("Slope", ifelse(summary(model_CTh)$coefficients[5,1]<0.001," < 0.001",
                                            paste0(" = ", round(summary(model_CTh)$coefficients[5,1], 3))),      # slope for Age:grp
                            "\nAdj R-squared", ifelse(summary(model_CTh)$adj.r.squared<0.001, " < 0.001",
                                                      paste0(" = ", (round(summary(model_CTh)$adj.r.squared,3)))),         # global model adj-R
                            "\np-value", ifelse(summary(model_CTh)$coefficients[5,4]<0.05,"<0.05",
                                                paste0(" = ", round(summary(model_CTh)$coefficients[5,4], 3)))       # p-val for Age:grp
    )

    plot_rh<-ggplot(data=subset(new_df, !is.na(new_df[[grp]])), aes(x=Real.Age)) +
      geom_point(aes(y=.data[[paste0(cMD_ROI_R,"_res")]], color="cMD"), alpha=0.2) +
      geom_smooth(aes(y=.data[[paste0(cMD_ROI_R,"_res")]], color="cMD", linetype=!!sym(grp)), method=loess) +
      geom_point(aes(y=.data[[paste0(CTh_ROI_R,"_res")]], color="CTh"), alpha=0.2) +
      geom_smooth(aes(y=.data[[paste0(CTh_ROI_R,"_res")]], color="CTh", linetype=!!sym(grp)), method=loess) +
      theme_minimal() + scale_color_manual(name="Brain metric",
                                           values=c("cMD"="#f18f01","CTh"="#00296b")) +
      labs(x="Age (yrs)", y="Brain metric (z-scores)",  shape=grp, linetype=grp) + ylim(-4,5) +
      annotate("text", x=(ifelse(grp=="ApoE4",51,47)), y=-3, label=coef_text_cMD, size=5, color="#f18f01") +
      annotate("text", x=(ifelse(grp=="ApoE4",73,62)), y=-3, label=coef_text_CTh, size=5, color="#00296b")
    
    if (grp=="Ptau_75"){
      
      formula_cMD<-paste0(paste0(ROI_name,".bilat"),"~Real.Age*factor(",grp,")+Sex")
      model_cMD<-lm(as.formula(formula_cMD), new_df_ptau)
      coef_text_cMD <- paste0(#"Slope", ifelse(summary(model_cMD)$coefficients[5,1]<0.001," < 0.001",
        #               paste0(" = ", round(summary(model_cMD)$coefficients[5,1], 3))),      # slope for Age:grp
        "\nR-adj", ifelse(summary(model_cMD)$adj.r.squared<0.001, " < 0.001",
                          paste0(" = ", (round(summary(model_cMD)$adj.r.squared,3)))),         # global model adj-R
        "\nP-val", ifelse(summary(model_cMD)$coefficients[5,4]<0.05,"<0.05",
                          paste0(" = ", round(summary(model_cMD)$coefficients[5,4], 3)))       # p-val for Age:grp
      )
      formula_CTh<-paste0(paste0(ROI_name,"_bilat"),"~Real.Age*factor(",grp,")+Sex")
      model_CTh<-lm(as.formula(formula_CTh), new_df_ptau)
      coef_text_CTh <- paste0(#"Slope", ifelse(summary(model_CTh)$coefficients[5,1]<0.001," < 0.001",
        #               paste0(" = ", round(summary(model_CTh)$coefficients[5,1], 3))),      # slope for Age:grp
        "\nR-adj", ifelse(summary(model_CTh)$adj.r.squared<0.001, " < 0.001",
                          paste0(" = ", (round(summary(model_CTh)$adj.r.squared,3)))),         # global model adj-R
        "\nP-val", ifelse(summary(model_CTh)$coefficients[5,4]<0.05,"<0.05",
                          paste0(" = ", round(summary(model_CTh)$coefficients[5,4], 3)))       # p-val for Age:grp
      )
      
      plot_bilat<-ggplot(data=subset(new_df_ptau, !is.na(new_df_ptau[[grp]])), aes(x=Real.Age)) +
        geom_point(aes(y=.data[[paste0(ROI_name,".bilat_res")]], color="cMD"), alpha=0.2) +
        geom_smooth(aes(y=.data[[paste0(ROI_name,".bilat_res")]], color="cMD", linetype=!!sym(grp)), method=loess) +
        geom_point(aes(y=.data[[paste0(ROI_name,"_bilat_res")]], color="CTh"), alpha=0.2) +
        geom_smooth(aes(y=.data[[paste0(ROI_name,"_bilat_res")]], color="CTh", linetype=!!sym(grp)), method=loess) +
        theme_minimal() + scale_color_manual(name="Brain metric",
                                             values=c("cMD"="#f18f01","CTh"="#00296b")) +
        labs( x="Age (yrs)", y="Brain metric (z-scores)",  shape=grp, linetype=grp) + ylim(-4,5) +
        annotate("text", x=(ifelse(grp=="ApoE4",51,47)), y=-3, label=coef_text_cMD, size=8, color="#f18f01") +
        annotate("text", x=(ifelse(grp=="ApoE4",73,62)), y=-3, label=coef_text_CTh, size=8, color="#00296b") +
        theme(text=element_text(size=12), axis.text=element_text(size=18), axis.title=element_text(size=21))
      
    }else{
      
      formula_cMD<-paste0(paste0(ROI_name,".bilat"),"~Real.Age*factor(",grp,")+Sex")
      model_cMD<-lm(as.formula(formula_cMD), new_df)
      coef_text_cMD <- paste0(#"Slope", ifelse(summary(model_cMD)$coefficients[5,1]<0.001," < 0.001",
        #               paste0(" = ", round(summary(model_cMD)$coefficients[5,1], 3))),      # slope for Age:grp                
        "\nR-adj", ifelse(summary(model_cMD)$adj.r.squared<0.001, " < 0.001",
                          paste0(" = ", (round(summary(model_cMD)$adj.r.squared,3)))),         # global model adj-R
        "\nP-val", ifelse(summary(model_cMD)$coefficients[5,4]<0.05,"<0.05",
                          paste0(" = ", round(summary(model_cMD)$coefficients[5,4], 3)))       # p-val for Age:grp
      )
      formula_CTh<-paste0(paste0(ROI_name,"_bilat"),"~Real.Age*factor(",grp,")+Sex")
      model_CTh<-lm(as.formula(formula_CTh), new_df)
      coef_text_CTh <- paste0(#"Slope", ifelse(summary(model_CTh)$coefficients[5,1]<0.001," < 0.001",
        #               paste0(" = ", round(summary(model_CTh)$coefficients[5,1], 3))),      # slope for Age:grp                
        "\nR-adj", ifelse(summary(model_CTh)$adj.r.squared<0.001, " < 0.001",
                          paste0(" = ", (round(summary(model_CTh)$adj.r.squared,3)))),         # global model adj-R
        "\nP-val", ifelse(summary(model_CTh)$coefficients[5,4]<0.05,"<0.05",
                          paste0(" = ", round(summary(model_CTh)$coefficients[5,4], 3)))       # p-val for Age:grp
      )
      
      plot_bilat<-ggplot(data=subset(new_df, !is.na(new_df[[grp]])), aes(x=Real.Age)) +
        geom_point(aes(y=.data[[paste0(ROI_name,".bilat_res")]], color="cMD"), alpha=0.2) +
        geom_smooth(aes(y=.data[[paste0(ROI_name,".bilat_res")]], color="cMD", linetype=!!sym(grp)), method=loess) +
        geom_point(aes(y=.data[[paste0(ROI_name,"_bilat_res")]], color="CTh"), alpha=0.2) + 
        geom_smooth(aes(y=.data[[paste0(ROI_name,"_bilat_res")]], color="CTh", linetype=!!sym(grp)), method=loess) + 
        theme_minimal() + scale_color_manual(name="Brain metric",
                                             values=c("cMD"="#f18f01","CTh"="#00296b")) +
        labs( x="Age (yrs)", y="Brain metric (z-scores)",  shape=grp, linetype=grp) + ylim(-4,5) +
        annotate("text", x=(ifelse(grp=="ApoE4",51,47)), y=-3, label=coef_text_cMD, size=8, color="#f18f01") +
        annotate("text", x=(ifelse(grp=="ApoE4",73,62)), y=-3, label=coef_text_CTh, size=8, color="#00296b") +
        theme(text=element_text(size=12), axis.text=element_text(size=18), axis.title=element_text(size=21))
      
    }
    
    lh_plots[[paste0(ROI_name,"_",grp)]]<-plot_lh
    rh_plots[[paste0(ROI_name,"_",grp)]]<-plot_rh
    bilat_plots[[paste0(ROI_name,"_",grp)]]<-plot_bilat
    
    ggsave(filename=paste0("~/oriol/Uni/Doctorat/Projectes/TFM Expanded/Resultats/",grp,"_",ROI_name,".tiff"),
           plot=plot_bilat,
           width=(2.72*3.2), height=(1.93*3.2), units="in", dpi=300, device="tiff", bg="white")
    
  }
}

