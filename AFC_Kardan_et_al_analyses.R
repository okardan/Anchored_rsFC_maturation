########################
# Title: Assessing neurocognitive maturation in early adolescence based on baby and adult functional brain landscapes
# Contact: Omid Kardan omidk@med.umich.edu
# This script compiles all non-brain variables from the ABCD curated files(R5.0) and uses them in the regressions and the mediations
# with AFC 
# 

# Requires downloaded ABCD data tables: abcd_p_demo.csv, abcd_y_lt.csv, abcd_ant01.csv, ph_p_pds.csv, ph_y_pds.csv (from Release 5.0) and 
# abcd_mrinback02.csv (from Release 4.0)

# Also requires abcd_sub_event_list.csv:
# create a file named "abcd_sub_event_list.csv" with two columns. one column containing all ABCD subids repeated in 5 rows per sub (one row for each year) and 2nd column with years labeled
# as 'eventname' with levels: 'baseline_year_1_arm_1', '1_year_follow_up_y_arm_1', '2_year_follow_up_y_arm_1', '3_year_follow_up_y_arm_1', '4_year_follow_up_y_arm_1'

# Also requires mod_sumrs_baseline_rest_mc2024_allruns.csv and mod_sumrs_2Y_rest_mc2024_allruns.csv which contain the brain metrics and are 
# produced by the Matlab script calc_AFC_Kardan_et_al.m

library(dplyr)
library(ggplot2)
library(tidyverse)
library(ppcor)
library(fastDummies)
library(sjPlot)
library(corrplot)
library(lme4)
library(plyr)
library(mediation)
library(ggExtra)


setwd("MyDir")
dat_list <- read.csv('abcd_sub_event_list.csv') # a csv file containing all ABCD subids repeated in 4 rows per sub (one row for each year with years labeled as 'eventname')
##################### demog and race/ethnicity #############################
datafordemo0 <- read.csv('MyDir/abcd_p_demo.csv')
datafordemo <- merge(dat_list,datafordemo0[,c('src_subject_id','eventname','demo_comb_income_v2_l',
                                              'demo_prnt_ed_v2_l','demo_prtnr_ed_v2_l',
                                              'race_ethnicity','demo_sex_v2')],
                     by = c('src_subject_id','eventname'), all.x = TRUE)
datafordemo$Income <- datafordemo$demo_comb_income_v2_l
datafordemo$Income[datafordemo$Income==999] <- NA
datafordemo$Income[datafordemo$Income==777] <- NA
datafordemo$Income_cat = factor( datafordemo$Income, levels= 1:10, 
                                 labels = c("5000", "8500", "14000", "20500", "30000",
                                            "42500", "62500", "87500", "150000", "200000") )
datafordemo$HighestEdParent <- datafordemo$demo_prnt_ed_v2_l 
datafordemo$HighestEdParent <- as.numeric(datafordemo$HighestEdParent)
datafordemo$HighestEdParent[datafordemo$HighestEdParent == 999] <- NA
datafordemo$HighestEdParent[datafordemo$HighestEdParent == 777] <- NA
datafordemo$HighestEdPartner <- datafordemo$demo_prtnr_ed_v2_l
datafordemo$HighestEdPartner <- as.numeric(datafordemo$HighestEdPartner)
datafordemo$HighestEdPartner[datafordemo$HighestEdPartner == 999] <- NA
datafordemo$HighestEdPartner[datafordemo$HighestEdPartner == 777] <- NA
#retain the highest education out of the two parents/partners
datafordemo$HighestEd <- pmax(datafordemo$HighestEdParent, datafordemo$HighestEdPartner, na.rm =TRUE)
datafordemo$Male_bin = ifelse(datafordemo$demo_sex_v2 == 1 | datafordemo$demo_sex_v2 == 3,1,0) # 3 intsex_male and no intsex_fem
#dummy code race
datafordemo <- datafordemo %>% 
  dummy_cols(ignore_na = TRUE, select_columns = c("race_ethnicity"))
colnames(datafordemo)[15] = c("White") 
colnames(datafordemo)[16] = c("Black")
colnames(datafordemo)[17] = c("Hispanic")
colnames(datafordemo)[18] = c("Asian")
colnames(datafordemo)[19] = c("Other")

datafordem <- datafordemo[,c('src_subject_id','eventname','Income','HighestEd','Male_bin',
                             'White','Black','Hispanic','Asian','Other')]

##################### age and site and relatives and height and puberty #
dataforage0 <- read.csv('MyDir/abcd_y_lt.csv')
dataforage <- merge(dat_list,dataforage0[,c('src_subject_id','eventname','site_id_l',
                                            'interview_age','rel_family_id')],
                    by = c('src_subject_id','eventname'), all.x = TRUE)

dataforHeigth0 <- read.csv('MyDir/abcd_ant01.csv')
dataforHeight <- merge(dat_list,dataforHeigth0[,c('src_subject_id','eventname','anthroheightcalc')],
                       by = c('src_subject_id','eventname'), all.x = TRUE)


dataforPubDev_p <- read.csv('MyDir/ph_p_pds.csv') # pds_p_ss_male_category_2 and pds_p_ss_female_category_2
dataforPubDev_p <- dataforPubDev_p %>% mutate(pds_p_ss = coalesce(pds_p_ss_female_category_2, pds_p_ss_male_category_2))

dataforPubDev_y <- read.csv('MyDir/ph_y_pds.csv')
dataforPubDev_y <- dataforPubDev_y %>% mutate(pds_y_ss = coalesce(pds_y_ss_female_category_2, pds_y_ss_male_cat_2))

dataforPubDev0 <- merge(dataforPubDev_p[,c("subid","eventname","pds_p_ss")], dataforPubDev_y[,c("subid","eventname","pds_y_ss")], by = c("subid","eventname"))
dataforPubDev0$pds_ss <- rowMeans(dataforPubDev0[,c("pds_p_ss","pds_y_ss")], na.rm = TRUE)

dataforPubDev <- merge(dat_list,dataforPubDev0[,c('subid','eventname','pds_ss')],
                       by = c('subid','eventname'), all.x = TRUE)

##################### N-back #
cog0 <- read.csv('MyDir/abcd_mrinback02.csv')
cog0$nback0_acc <- cog0$The.rate.of.correct.responses.to.0.back.stimuli.during.run.1.and.run.2
cog0$nback2_acc <- cog0$The.rate.of.correct.responses.to.2.back.stimuli.during.run.1.and.run.2
cog_nbk <- merge(dat_list,cog0[,c('src_subject_id','eventname','nback0_acc','nback2_acc')],
                 by = c('src_subject_id','eventname'), all.x = TRUE)

########## combine all data 
df_total_withcog <- datafordem %>% 
  full_join(dataforage, by = c('src_subject_id','eventname'))%>%
  full_join(dataforHeight, by = c('src_subject_id','eventname'))%>%
  full_join(dataforPubDev, by = c('src_subject_id','eventname'))%>%
  full_join(cog_nbk, by = c('src_subject_id','eventname'))
write.csv(df_total_withcog,'MyDir/df_total_withcog2024.csv')  


##### spread the ses and save as seperate years
df_total_withcog <- read.csv('MyDir/df_total_withcog2024.csv') 

dy <- df_total_withcog %>% mutate(Income=as.numeric(Income),HighestEd=as.numeric(HighestEd), Male_bin=as.numeric(Male_bin), White=as.numeric(White),
                                  Black=as.numeric(Black), Hispanic=as.numeric(Hispanic), Asian=as.numeric(Asian), Other=as.numeric(Other))

dy <- dy %>% group_by(src_subject_id) %>% mutate_at(vars(Income, HighestEd, Male_bin, White, Black, Hispanic, Asian, Other),
                                                    ~replace_na(.,mean(.,na.rm=TRUE))) # average the SES measures if they are available for multiple years
write.csv(dy[dy$eventname == 'baseline_year_1_arm_1',],'MyDir/df_filled_demog_withcog_Y0.csv') 
write.csv(dy[dy$eventname == '1_year_follow_up_y_arm_1',],'MyDir/df_filled_demog_withcog_Y1.csv') 
write.csv(dy[dy$eventname == '2_year_follow_up_y_arm_1',],'MyDir/df_filled_demog_withcog_Y2.csv') 

##########################                                    ################################
##########################                                  ####################################
##########################                     ##################################################
# read the csv files made above

df_filled_demog_withcog_Y0 <- read.csv('MyDir/df_filled_demog_withcog_Y0.csv')
df_filled_demog_withcog_Y2 <- read.csv('MyDir/df_filled_demog_withcog_Y2.csv')
df_filled_demog_withcog_Y2$rel_family_id <- df_filled_demog_withcog_Y0$rel_family_id

# read the csv files made in the Matlab script

wbdiff_b <- read.csv('MyDir/mod_sumrs_baseline_rest_mc2024_allruns.csv')
wbdiff_all_b <- merge(wbdiff_b,df_filled_demog_withcog_Y0[,c('subid','interview_age','Income','HighestEd',
                                                             'Male_bin','White','Black','Hispanic','Asian','Other',
                                                             'site_id_l','rel_family_id','anthroheightcalc','pds_ss',
                                                             'nback0_acc','nback2_acc')],
                      by.x = 'subids', by.y = 'subid', all.x = TRUE)  # 

wbdiff_all_b <- wbdiff_all_b[!is.nan(wbdiff_all_b$wbdiff_adult),] 

wbdiff_y2 <- read.csv('MyDir/mod_sumrs_2Y_rest_mc2024_allruns.csv')
wbdiff_all_y2 <- merge(wbdiff_y2,df_filled_demog_withcog_Y2[,c('subid','interview_age','Income','HighestEd',
                                                               'Male_bin','White','Black','Hispanic','Asian','Other',
                                                               'site_id_l','rel_family_id','anthroheightcalc','pds_ss',
                                                               'nback0_acc','nback2_acc')],
                       by.x = 'subids', by.y = 'subid', all.x = TRUE)  # 

wbdiff_all_y2 <- wbdiff_all_y2[!is.nan(wbdiff_all_y2$wbdiff_adult),] 


#############################
#############################
#                 Main Results 
#############################
#############################

dat1 <- wbdiff_all_b %>% mutate(mat_score = wbdiff_adult - wbdiff_baby, age = interview_age/12, age_q = (scale(interview_age))^2)
dat1age <- dat1[!is.na(dat1$age) &  !is.na(dat1$pFD) & !is.na(dat1$mat_score) & dat1$nruns>1,] 
dat2 <- wbdiff_all_y2 %>% mutate(mat_score = wbdiff_adult - wbdiff_baby, age = interview_age/12, age_q = (scale(interview_age))^2)
dat2age <- dat2[!is.na(dat2$age) &  !is.na(dat2$pFD) & !is.na(dat2$mat_score)& dat2$nruns>1,] 
dat1age$NBK <- rowMeans(dat1age[,c('nback0_acc','nback2_acc')], na.rm=T)
dat2age$NBK <- rowMeans(dat2age[,c('nback0_acc','nback2_acc')], na.rm=T)

dat_both <- merge(dat1age[,c('subids','rel_family_id' ,'mat_score','pFD','nruns','Male_bin','Income','pds_ss','age','anthroheightcalc','HighestEd','site_id_l','White',
                             'Black','Hispanic','Asian','Other','NBK')],
                  dat2age[,c('subids','rel_family_id','mat_score','pFD','nruns','Male_bin','Income','pds_ss','age','anthroheightcalc','HighestEd','site_id_l','White',
                             'Black','Hispanic','Asian','Other','NBK')],by = 'subids',all = TRUE)  # n = 3292

dat1age<- dat1age %>% mutate(agegroup="Y0")
dat2age <- dat2age %>% mutate(agegroup="Y2")

dat_both$delta_mat_score <- dat_both$mat_score.y - dat_both$mat_score.x
dat_both$delta_NBK <- dat_both$NBK.y - dat_both$NBK.x


lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

###  Results section 3.2. Anchored index of functional maturation predicts cognitive performance across and within youth

model0 <- lmer(NBK.x ~  1 + mat_score.x + pFD.x  +nruns.x +  +Male_bin.x + White.x +
                 Black.x + Hispanic.x + Asian.x+ (1|site_id_l.x), data = dat_both)
lm.beta.lmer(model0)
summary(model0)  # n = 6489, beta = .063
dat_y0_complete <- dat_both %>% drop_na(NBK.x, mat_score.x, pFD.x, nruns.x, Male_bin.x, White.x, Black.x, Hispanic.x, Asian.x)
length(dat_y0_complete$subids[dat_y0_complete$White.x == 1])/6489 # 57.8% white
length(dat_y0_complete$subids[dat_y0_complete$Black.x == 1])/6489 # 12.0% Black
length(dat_y0_complete$subids[dat_y0_complete$Hispanic.x == 1])/6489 # 18.2% Hispanic
length(dat_y0_complete$subids[dat_y0_complete$Asian.x == 1])/6489 # 2.0% Asian
length(dat_y0_complete$subids[dat_y0_complete$Other.x == 1])/6489 # 10.0% other
length(dat_y0_complete$subids[dat_y0_complete$Male_bin.x == 1])/6489 # 48.5% male


model2 <- lmer(NBK.y ~  1 + mat_score.y + pFD.y  +nruns.y  +Male_bin.y + White.y +
                 Black.y + Hispanic.y + Asian.y+ (1|site_id_l.y), data = dat_both)
lm.beta.lmer(model2)
summary(model2)  # n = 5089, beta = .099
dat_y2_complete <- dat_both %>% drop_na(NBK.y, mat_score.y, pFD.y, nruns.y, Male_bin.y, White.y, Black.y, Hispanic.y, Asian.y)
length(dat_y2_complete$subids[dat_y2_complete$White.y == 1])/5089 # 56.7% white
length(dat_y2_complete$subids[dat_y2_complete$Black.y == 1])/5089 # 12.0% Black
length(dat_y2_complete$subids[dat_y2_complete$Hispanic.y == 1])/5089 # 19.3% Hispanic
length(dat_y2_complete$subids[dat_y2_complete$Asian.y == 1])/5089 # 1.7% Asian
length(dat_y2_complete$subids[dat_y2_complete$Other.y == 1])/5089 # 10.3% other
length(dat_y2_complete$subids[dat_y2_complete$Male_bin.y == 1])/5089 # 50.7% male

model02 <- lmer(delta_NBK ~  1 + delta_mat_score + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both)

lm.beta.lmer(model02)  # not additionally adjusted for baseline N-back acc
summary(model02)  # n = 3200, beta = .070
model02_bl <- lmer(delta_NBK ~  1 + NBK.x + delta_mat_score + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                     Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both)

lm.beta.lmer(model02_bl) # additionally adjusted for baseline N-back acc
summary(model02)  # n = 3200, beta = .071
model02_bl <- lmer(delta_NBK ~  1 + NBK.x + mat_score.x + delta_mat_score + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                     Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both)

lm.beta.lmer(model02_bl) # additionally adjusted for baseline N-back acc and baseline AFC
summary(model02)  # n = 3200, beta = .094
dat_both_complete <- dat_both %>% drop_na(delta_NBK, delta_mat_score, pFD.x, pFD.y, nruns.x, nruns.y, Male_bin.x, White.x, Black.x, Hispanic.x, Asian.x)
length(dat_both_complete$subids[dat_both_complete$White.x == 1])/3200 # 61.3% white
length(dat_both_complete$subids[dat_both_complete$Black.x == 1])/3200 # 9.9% Black
length(dat_both_complete$subids[dat_both_complete$Hispanic.x == 1])/3200 # 18.0% Hispanic
length(dat_both_complete$subids[dat_both_complete$Asian.x == 1])/3200 # 1.5% Asian
length(dat_both_complete$subids[dat_both_complete$Other.x == 1])/3200 # 9.3% other
length(dat_both_complete$subids[dat_both_complete$Male_bin.x == 1])/3200 # 49.4% male


t.test(dat_y2_complete$mat_score.y, y = dat_y0_complete$mat_score.x,pair=FALSE) # Welch t = 5.62, p<.001
t.test(dat_both_complete$mat_score.y, y = dat_both_complete$mat_score.x,pair=TRUE) # paired t = 9.24, p<.001

## making Figure 4
library(ggExtra)
dat1age<- dat1age %>% mutate(dataset="1")
dat2age <- dat2age %>% mutate(dataset="2")

dat1age_temp <- dat1age %>% dplyr::select(NBK, mat_score, dataset)
dat2age_temp <- dat2age %>% dplyr::select(NBK, mat_score, dataset)
dat12age <- rbind(dat1age_temp, dat2age_temp)


p12 <- ggplot(dat12age,aes(x = mat_score, y = (NBK)*100/1, color= dataset, group=dataset)) + 
  geom_point(alpha=.20, size = 4) +
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#E7B800","#c6a2c6"), aesthetics = "color") +
  theme_classic(base_size = 16) +  theme_classic(base_size = 16) + 
  xlim(c(-.5,.6)) +
  ylim(c(0,100)) +
  theme(legend.position = "none") +
  labs(y = "N-back Accuracy (%)", x = 'Anchored rsFC maturation score', color="", title = "") 

p12
ggsave(file = "Fig4A_AFC_Y0_Y2.png", ggMarginal(p12, groupColour = TRUE, groupFill = TRUE), width = 7, height = 7)

m1<- mean(dat_both$delta_mat_score[!is.na(dat_both$delta_mat_score) & !is.na(dat_both$delta_NBK)])
m2<- mean(dat_both$delta_NBK[!is.na(dat_both$delta_mat_score) & !is.na(dat_both$delta_NBK)]*100)
mm <- as.data.frame(rbind(c(m1,m2),c(0,0)))
p1 <- ggplot(data=dat_both) + 
  aes(x = delta_mat_score, y = delta_NBK*100) + 
  geom_point(alpha=.06, size = 1.5,color = "#00AFBB") +
  geom_point(data = mm,aes(x = V1,
                           y = V2, size=3),shape = c(4,3), size=c(3,1), color = c("red","black"), inherit.aes = FALSE)+
  geom_smooth(method = "lm", aes(x = delta_mat_score, y = delta_NBK*100), size=1.5, color = "#00AFBB") +
  theme_classic(base_size = 16) + coord_cartesian(xlim = c(-.45,0.6), ylim = c(-30,30))+
  labs(x = "Change in anchored rsFC maturation score", y = "Change in N-back accuracy (%)",  title = " ") 
p1
ggsave(file = 'Fig4B_delta_AFC_scatter.png',ggMarginal(p1), dpi = 600, width = 7, height = 7)

###### Results section 3.3. AFC predicts cognitive performance above and beyond family resources
model0 <- lmer(mat_score.x ~  1 + pFD.x  +nruns.x +Male_bin.x + White.x + 
                 Black.x + Hispanic.x + Asian.x+ Income.x + HighestEd.x +  (1|site_id_l.x), data = dat_both)
lm.beta.lmer(model0)
summary(model0)  

model2 <- lmer(mat_score.y ~  1 + pFD.y  +nruns.y +Male_bin.y + White.y + 
                 Black.y + Hispanic.y + Asian.y + Income.y + HighestEd.y + (1|site_id_l.y), data = dat_both)
lm.beta.lmer(model2)
summary(model2)

model02 <- lmer(delta_mat_score ~  1 + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                  Black.x + Hispanic.x + Asian.x + Income.x + HighestEd.x  
                + (1|site_id_l.x), data = dat_both)

lm.beta.lmer(model02)
summary(model02) 

boot_func0 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(NBK.x ~  1 + mat_score.x +  pFD.x + +nruns.x + 
                              Male_bin.x + White.x +  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  lmwithage_resampled <- lmer(NBK.x ~  1 + mat_score.x  + pFD.x + nruns.x + 
                                Income.x + HighestEd.x +
                                Male_bin.x + White.x +
                                Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - unname(bet2[1])
}

boot_coefs0 <- replicate(1000, boot_func0(dat_both))
(pval = 1- sum(boot_coefs0>0)/length(boot_coefs0))  # p = 0.107 for change in beta
median(boot_coefs0) # 0.0068


boot_func2 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(NBK.y ~  1 + mat_score.y +  pFD.y + +nruns.y + 
                              Male_bin.y + White.y +  Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  lmwithage_resampled <- lmer(NBK.y ~  1 + mat_score.y  + pFD.y + nruns.y + 
                                Income.y + HighestEd.y +
                                Male_bin.y + White.y +
                                Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - unname(bet2[1])
}

boot_coefs2 <- replicate(1000, boot_func2(dat_both))
(pval = 1- sum(boot_coefs2>0)/length(boot_coefs2))  # p = 0.347 for change in beta
median(boot_coefs2) # 0.0018


boot_func <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(delta_NBK ~  1 + delta_mat_score +  pFD.x + pFD.y +nruns.x + nruns.y +
                              Male_bin.x + White.x +  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  lmwithage_resampled <- lmer(delta_NBK ~  1 + delta_mat_score +  pFD.x + pFD.y +nruns.x + nruns.y +
                                Income.x + HighestEd.x +
                                Male_bin.x + White.x +
                                Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - unname(bet2[1])
}

boot_coefs <- replicate(1000, boot_func(dat_both))
(pval = 1- sum(boot_coefs>0)/length(boot_coefs))  # p = 0.174 for change in beta
median(boot_coefs)  # 0.0045


###### Results section 3.4. AFC as a statistical mediator of age's association with cognitive task performance
model0 <- lmer(mat_score.x ~  1 + pFD.x  +nruns.x +Male_bin.x + White.x + 
                 Black.x + Hispanic.x + Asian.x+ Income.x + HighestEd.x + anthroheightcalc.x + age.x + pds_ss.x + (1|site_id_l.x), data = dat_both)
lm.beta.lmer(model0)
summary(model0) 
model00 <- lmer(NBK.x ~  1 + mat_score.x +  pFD.x + +nruns.x + Income.x + HighestEd.x + + anthroheightcalc.x + age.x + pds_ss.x +
                  Male_bin.x + White.x +  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both)
lm.beta.lmer(model00)
summary(model00) 

boot_func0 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(NBK.x ~  1 + mat_score.x +  pFD.x + +nruns.x + Income.x + HighestEd.x +
                              Male_bin.x + White.x +  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  lmwithage_resampled <- lmer(NBK.x ~  1 + mat_score.x  + pFD.x + nruns.x + Income.x + HighestEd.x
                              + anthroheightcalc.x + age.x + pds_ss.x +
                                Male_bin.x + White.x +
                                Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - unname(bet2[1])
}

boot_coefs0 <- replicate(1000, boot_func0(dat_both))
(pval = 1- sum(boot_coefs0>0)/length(boot_coefs0))  # p = 0.003 for change in beta of AFC predicting N-back
median(boot_coefs0) # 0.007

# Mediation in figure 5A
mediation_func0 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmapath_resampled <- lmer(mat_score.x ~  1 + age.x + anthroheightcalc.x  + pds_ss.x + pFD.x + +nruns.x + Income.x + HighestEd.x +
                              Male_bin.x + White.x +  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  lmbpath_resampled <- lmer(NBK.x ~  1 + mat_score.x + age.x + anthroheightcalc.x +  pds_ss.x + pFD.x + nruns.x + Income.x + HighestEd.x +
                              Male_bin.x + White.x + Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  lmnomed_resampled <- lmer(NBK.x ~  1 + age.x + anthroheightcalc.x +  pds_ss.x + pFD.x + nruns.x + Income.x + HighestEd.x +
                              Male_bin.x + White.x + Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  
  a <- lm.beta.lmer(lmapath_resampled)
  b <-  lm.beta.lmer(lmbpath_resampled)
  c <- lm.beta.lmer(lmnomed_resampled)
  c(unname(a[1]), unname(b[1]), unname(c[1]), unname(b[2]))
}

med_coefs0 <- replicate(1000, mediation_func0(dat_both))
abccp0 <- as.data.frame(t(med_coefs0))
abccp0$ind <- abccp0$V1 * abccp0$V2
abccp0$dir <- abccp0$V3 - abccp0$V4
(mean(abccp0$V1)) # .0374
(pnorm(mean(abccp0$V1)/sd(abccp0$V1),lower.tail = FALSE)) # p = .004
(mean(abccp0$V2)) # .048
(pnorm(mean(abccp0$V2)/sd(abccp0$V2),lower.tail = FALSE)) # p = .0003
(mean(abccp0$V3))  #  .170
(pnorm(mean(abccp0$V3)/sd(abccp0$V3),lower.tail = FALSE)) # p = .0000
(mean(abccp0$V4))  # .168
(pnorm(mean(abccp0$V4)/sd(abccp0$V4),lower.tail = FALSE)) # p = .00000
(mean(abccp0$ind))  # .002
(pnorm(mean(abccp0$ind)/sd(abccp0$ind),lower.tail = FALSE))  # p = 0.0207 for a*b
(mean(abccp0$dir))  # .002
(pnorm(mean(abccp0$dir)/sd(abccp0$dir),lower.tail = FALSE))  # p = 0.0267 for c-c'


# Mediation in figure 5B
mediation_func2 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmapath_resampled <- lmer(mat_score.y ~  1 + age.y + anthroheightcalc.y  + pds_ss.y + pFD.y + +nruns.y + Income.y + HighestEd.y +
                              Male_bin.y + White.y +  Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  lmbpath_resampled <- lmer(NBK.y ~  1 + mat_score.y + age.y + anthroheightcalc.y +  pds_ss.y + pFD.y + nruns.y + Income.y + HighestEd.y +
                              Male_bin.y + White.y + Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  lmnomed_resampled <- lmer(NBK.y ~  1 + age.y + anthroheightcalc.y +  pds_ss.y + pFD.y + nruns.y + Income.y + HighestEd.y +
                              Male_bin.y + White.y + Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  
  a <- lm.beta.lmer(lmapath_resampled)
  b <-  lm.beta.lmer(lmbpath_resampled)
  c <- lm.beta.lmer(lmnomed_resampled)
  c(unname(a[1]), unname(b[1]), unname(c[1]), unname(b[2]))
}

med_coefs2 <- replicate(1000, mediation_func2(dat_both))
abccp2 <- as.data.frame(t(med_coefs2))
abccp2$ind <- abccp2$V1 * abccp2$V2
abccp2$dir <- abccp2$V3 - abccp2$V4
(mean(abccp2$V1)) # .0257
(pnorm(mean(abccp2$V1)/sd(abccp2$V1),lower.tail = FALSE)) # p = .049
(mean(abccp2$V2)) # .096
(pnorm(mean(abccp2$V2)/sd(abccp2$V2),lower.tail = FALSE)) # p = .0000
(mean(abccp2$V3))  #  .080
(pnorm(mean(abccp2$V3)/sd(abccp2$V3),lower.tail = FALSE)) # p = .0000
(mean(abccp2$V4))  # .079
(pnorm(mean(abccp2$V4)/sd(abccp2$V4),lower.tail = FALSE)) # p = .00000
(mean(abccp2$ind))  # .0024
(pnorm(mean(abccp2$ind)/sd(abccp2$ind),lower.tail = FALSE))  # p = 0.0618 for a*b
(mean(abccp2$dir))  # .0026
(pnorm(mean(abccp2$dir)/sd(abccp2$dir),lower.tail = FALSE))  # p = 0.0514 for c-c'


# analyses with change in age and change in AFC
dat_both$delta_pds <- dat_both$pds_ss.y - dat_both$pds_ss.x
dat_both$mean_pds <- (dat_both$pds_ss.y + dat_both$pds_ss.x)/2
dat_both$delta_height <- dat_both$anthroheightcalc.y - dat_both$anthroheightcalc.x
dat_both$mean_height <- (dat_both$anthroheightcalc.y + dat_both$anthroheightcalc.x)/2
dat_both$delta_age <- dat_both$age.y - dat_both$age.x
(mean(dat_both$delta_age,na.rm = TRUE)*12)
(sd(dat_both$delta_age,na.rm = TRUE)*12)
dat_both$mean_age <- (dat_both$age.y + dat_both$age.x)/2

model02 <- lmer(delta_mat_score ~  1 + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                  Black.x + Hispanic.x + Asian.x + Income.x + HighestEd.x + pds_ss.x + delta_pds + age.x + delta_age + anthroheightcalc.x + delta_height 
                + (1|site_id_l.x), data = dat_both)

lm.beta.lmer(model02)
summary(model02)  
tab_model(model02, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE")

# Mediation model in Figure 6

mediation_func02 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmapath_resampled <- lmer(delta_mat_score ~  1 + delta_age + age.x + mat_score.x + NBK.x + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                              Black.x + Hispanic.x + Asian.x + Income.x + HighestEd.x + pds_ss.x + delta_pds  + anthroheightcalc.x + delta_height
                            + (1|site_id_l.x), data = resampled_data)
  lmbpath_resampled <- lmer(delta_NBK ~  1 + delta_mat_score + delta_age + age.x + mat_score.x +NBK.x+ pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                              Black.x + Hispanic.x + Asian.x + Income.x + HighestEd.x + pds_ss.x + delta_pds  + anthroheightcalc.x + delta_height
                            +(1|site_id_l.x), data = resampled_data)
  lmnomed_resampled <- lmer(delta_NBK ~  1 + delta_age + age.x + mat_score.x + NBK.x + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x + 
                              Black.x + Hispanic.x + Asian.x + Income.x + HighestEd.x + pds_ss.x + delta_pds  + anthroheightcalc.x + delta_height
                            + (1|site_id_l.x), data = resampled_data)
  
  a <- lm.beta.lmer(lmapath_resampled)
  b <-  lm.beta.lmer(lmbpath_resampled)
  c <- lm.beta.lmer(lmnomed_resampled)
  c(unname(a[1]), unname(b[1]), unname(c[1]), unname(b[2]))
}

med_coefs02 <- replicate(1000, mediation_func02(dat_both))
abccp02 <- as.data.frame(t(med_coefs02))
abccp02$ind <- abccp02$V1 * abccp02$V2
abccp02$dir <- abccp02$V3 - abccp02$V4
(mean(abccp02$V1)) # .0358
(pnorm(mean(abccp02$V1)/sd(abccp02$V1),lower.tail = FALSE)) # p = .017
(mean(abccp02$V2)) # .0756
(pnorm(mean(abccp02$V2)/sd(abccp02$V2),lower.tail = FALSE)) # p = .0000
(mean(abccp02$V3))  #  -0.000
(pnorm(mean(abccp02$V3)/sd(abccp02$V3),lower.tail = FALSE)) # p = .507
(mean(abccp02$V4))  # -0.003
(pnorm(mean(abccp02$V4)/sd(abccp02$V4),lower.tail = FALSE)) # p = .572
(mean(abccp02$ind))  # .0027
(pnorm(mean(abccp02$ind)/sd(abccp02$ind),lower.tail = FALSE))  # p = 0.0375 for a*b
(mean(abccp02$dir))  # .0027
(pnorm(mean(abccp02$dir)/sd(abccp02$dir),lower.tail = FALSE))  # p = 0.040 for c-c'


#############################
#############################
# Posthoc analyses and supplementary sensitivity analyses 
#############################
#############################
####### Results section 3.5. Does anchoring to baby networks matter?

# modularity
dat_both_modul <- merge(dat1age[,c('subids','rel_family_id' ,'mat_score','Qs_3','wbdiff_baby','wbdiff_adult'
                                   ,'pFD','nruns','Male_bin','Income','pds_ss','age','anthroheightcalc','HighestEd','site_id_l','White',
                                   'Black','Hispanic','Asian','Other','NBK')],
                        dat2age[,c('subids','rel_family_id','mat_score','Qs_3','wbdiff_baby','wbdiff_adult',
                                   'pFD','nruns','Male_bin','Income','pds_ss','age','anthroheightcalc','HighestEd','site_id_l','White',
                                   'Black','Hispanic','Asian','Other','NBK')],by = 'subids',all = TRUE)  # n = 3292

dat_both_modul$delta_mat_score <- dat_both_modul$mat_score.y - dat_both_modul$mat_score.x
dat_both_modul$delta_wbdiff_adult <- dat_both_modul$wbdiff_adult.y - dat_both_modul$wbdiff_adult.x
dat_both_modul$delta_wbdiff_baby <- dat_both_modul$wbdiff_baby.y - dat_both_modul$wbdiff_baby.x
dat_both_modul$delta_Q <- dat_both_modul$Qs_3.y - dat_both_modul$Qs_3.x
dat_both_modul$delta_NBK <- dat_both_modul$NBK.y - dat_both_modul$NBK.x

lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}

model02 <- lmer(delta_NBK ~  1 + delta_Q + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x +
                  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both_modul)
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}
lm.beta.lmer(model02)
summary(model02)  # n = 3200, beta = -.002, t = -0.121, p = .904
tab_model(model02, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("change in N-back Acc"))

model0 <- lmer(NBK.x ~  1 + Qs_3.x + pFD.x + nruns.x +  +Male_bin.x + White.x + 
                 Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both_modul)
lm.beta.lmer(model0)
summary(model0)  # n = 6489, beta = -.020, p = .097
tab_model(model0, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("change in N-back Acc"))

model2 <- lmer(NBK.y ~  1 + Qs_3.y +  pFD.y + + nruns.y +Male_bin.y + White.y +
                 Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = dat_both_modul)
lm.beta.lmer(model2)
summary(model2)  # n = 5089, beta = .002,  p = .885
tab_model(model2, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("change in N-back Acc"))


# within - between for only adult nets

model02 <- lmer(delta_NBK ~  1 + delta_wbdiff_adult + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x +
                  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both_modul)
lm.beta.lmer(model02)
summary(model02)  # n = 3200, beta = .006, p = .754
tab_model(model02, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("change in N-back Acc"))

model0 <- lmer(NBK.x ~  1 + wbdiff_adult.x + pFD.x + nruns.x +  +Male_bin.x + White.x +
                 Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both_modul)
lm.beta.lmer(model0)
summary(model0)  # n = 6489, beta = .014, p = .314
tab_model(model0, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("N-back Acc at Y0"))

model2 <- lmer(NBK.y ~  1 + wbdiff_adult.y +  pFD.y + nruns.y +Male_bin.y + White.y +
                 Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = dat_both_modul)
lm.beta.lmer(model2)
summary(model2)  # n = 5089, beta = .032,  p = .033
tab_model(model2, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("N-back Acc at Y2"))

# comparing adult-only and AFC
boot_func22 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(NBK.y ~  1 + mat_score.y +  pFD.y  +nruns.y +  +Male_bin.y + White.y +
                              Black.y + Hispanic.y + Asian.y+ (1|site_id_l.x), data = resampled_data)
  lmwithage_resampled <- lmer(NBK.y ~  1 + wbdiff_adult.y +  pFD.y + nruns.y +Male_bin.y + White.y +
                                Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - unname(bet2[1])
}

boot_coefs22 <- replicate(1000, boot_func22(dat_both_modul))
median(boot_coefs22)
(pval = 1- sum(boot_coefs22>0)/length(boot_coefs22))  # p = .000 for change in beta


# within - between for only baby nets

model02 <- lmer(delta_NBK ~  1 + delta_wbdiff_baby + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x +
                  Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both_modul)
lm.beta.lmer(model02)
summary(model02)  # n = 3200, beta = -.076, p = .000
tab_model(model02, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("change in N-back Acc"))

model0 <- lmer(NBK.x ~  1 + wbdiff_baby.x + pFD.x + nruns.x +  +Male_bin.x + White.x +
                 Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = dat_both_modul)
lm.beta.lmer(model0)
summary(model0)  # n = 6489, beta = -.061, p = .000
tab_model(model0, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("N-back Acc at Y0"))

model2 <- lmer(NBK.y ~  1 + wbdiff_baby.y +  pFD.y + nruns.y +Male_bin.y + White.y +
                 Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = dat_both_modul)
lm.beta.lmer(model2)
summary(model2)  # n = 5089, beta = -.079,  p = .000
tab_model(model2, show.df = TRUE, show.stat = TRUE, string.stat = "t-Statistic",show.se = FALSE, string.se = "SE",dv.labels = c("N-back Acc at Y2"))

# comparing baby-only and AFC
boot_func00 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(NBK.x ~  1 + mat_score.x +  pFD.x  +nruns.x +  +Male_bin.x + White.x +
                              Black.x + Hispanic.x + Asian.x+ (1|site_id_l.x), data = resampled_data)
  lmwithage_resampled <- lmer(NBK.x ~  1 + wbdiff_baby.x +  pFD.x + nruns.x +Male_bin.x + White.x +
                                Black.x + Hispanic.x + Asian.x + (1|site_id_l.x), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - (-1)*unname(bet2[1])
}

boot_coefs00 <- replicate(1000, boot_func00(dat_both_modul))
median(boot_coefs00)
(pval = 1- sum(boot_coefs00>0)/length(boot_coefs00))  # p = .450 for change in beta

boot_func22 <- function(data){
  resampled_data = data %>% 
    sample_frac(1, replace=TRUE)
  
  lmnoage_resampled <- lmer(NBK.y ~  1 + mat_score.y +  pFD.y  +nruns.y +  +Male_bin.y + White.y +
                              Black.y + Hispanic.y + Asian.y+ (1|site_id_l.x), data = resampled_data)
  lmwithage_resampled <- lmer(NBK.y ~  1 + wbdiff_baby.y +  pFD.y + nruns.y +Male_bin.y + White.y +
                                Black.y + Hispanic.y + Asian.y + (1|site_id_l.y), data = resampled_data)
  bet1 <- lm.beta.lmer(lmnoage_resampled)
  bet2 <-  lm.beta.lmer(lmwithage_resampled)
  unname(bet1[1]) - (-1)*unname(bet2[1])
}

boot_coefs22 <- replicate(1000, boot_func22(dat_both_modul))
median(boot_coefs22)
(pval = 1- sum(boot_coefs22>0)/length(boot_coefs22))  # p = .061 for change in beta



# supplementary 1 (exclude based on family id)
#############################
#############################
dat_both_famexclude <- dat_both %>% distinct(rel_family_id.x, .keep_all = TRUE)
model02 <- lmer(delta_NBK ~  1 + delta_mat_score + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x +
                  Black.x + Hispanic.x + Asian.x+ (1|site_id_l.x), data = dat_both_famexclude)
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}
lm.beta.lmer(model02)
summary(model02)  # n = 2739, beta = .055
dat_both_famexclude_complete <- dat_both_famexclude %>% drop_na(delta_NBK, delta_mat_score, pFD.x, pFD.y, nruns.x, nruns.y, Male_bin.x, White.x, Black.x, Hispanic.x, Asian.x)


dat_y0_famexclude <- dat_both %>% distinct(rel_family_id.x, .keep_all = TRUE)
model0 <- lmer(NBK.x ~  1 + mat_score.x + pFD.x  +nruns.x  +Male_bin.x + White.x +
                 Black.x + Hispanic.x + Asian.x+ (1|site_id_l.x), data = dat_y0_famexclude)
lm.beta.lmer(model0)
summary(model0)  # n = 5641, beta = .064
dat_y0_famexclude_complete <- dat_y0_famexclude %>% drop_na(NBK.x, mat_score.x, pFD.x, nruns.x, Male_bin.x, White.x, Black.x, Hispanic.x, Asian.x)


dat_y2_famexclude <- dat_both %>% distinct(rel_family_id.y, .keep_all = TRUE)
model2 <- lmer(NBK.y ~  1 + mat_score.y + pFD.y  +nruns.y  +Male_bin.y + White.y +
                 Black.y + Hispanic.y + Asian.y+ (1|site_id_l.y), data = dat_y2_famexclude)
lm.beta.lmer(model2)
summary(model2)  # n = 4494, beta = .095
dat_y2_famexclude_complete <- dat_y2_famexclude %>% drop_na(NBK.y, mat_score.y, pFD.y, nruns.y, Male_bin.y, White.y, Black.y, Hispanic.y, Asian.y)


t.test(dat_y2_famexclude_complete$mat_score.y, y = dat_y0_famexclude_complete$mat_score.x,pair=FALSE) # Welch t = 5.18, p<.001
t.test(dat_both_famexclude$mat_score.y, y = dat_both_famexclude$mat_score.x,pair=TRUE) # paired t = 8.67, p<.001
######################### Figure S1

dat1age_temp <- dat1age %>% dplyr::select(NBK, mat_score, dataset,rel_family_id)
dat2age_temp <- dat2age %>% dplyr::select(NBK, mat_score, dataset,rel_family_id)
dat1exclude_temp <- dat1age_temp %>% distinct(rel_family_id, .keep_all = TRUE)
dat2exclude_temp <- dat2age_temp %>% distinct(rel_family_id, .keep_all = TRUE)
dat12exclude <- rbind(dat1exclude_temp, dat2exclude_temp)


p12 <- ggplot(dat12exclude,aes(x = mat_score, y = (NBK)*100/1, color= dataset, group=dataset)) + 
  geom_point(alpha=.20, size = 4) +
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#E7B800","#c6a2c6"), aesthetics = "color") +
  theme_classic(base_size = 16) + 
  xlim(c(-.5,.6)) +
  ylim(c(0,100)) +
  theme(legend.position = "none") +
  labs(y = "N-back Accuracy (%)", x = 'Anchored rsFC maturation score', color="", title = "") 

p12
ggsave(file = "Fig_marg_supp1.png", ggMarginal(p12, groupColour = TRUE, groupFill = TRUE), width = 7, height = 7)


m1<- mean(dat_both_famexclude$delta_mat_score[!is.na(dat_both_famexclude$delta_mat_score) & !is.na(dat_both_famexclude$delta_NBK)])
m2<- mean(dat_both_famexclude$delta_NBK[!is.na(dat_both_famexclude$delta_mat_score) & !is.na(dat_both_famexclude$delta_NBK)]*100)
mm <- as.data.frame(rbind(c(m1,m2),c(0,0)))
p1 <- ggplot(data=dat_both_famexclude) + 
  aes(x = delta_mat_score, y = delta_NBK*100) + 
  geom_point(alpha=.06, size = 1.5,color = "#00AFBB") +
  geom_point(data = mm,aes(x = V1,
                           y = V2, size=3),shape = c(4,3), size=c(3,1), color = c("red","black"), inherit.aes = FALSE)+
  geom_smooth(method = "lm", aes(x = delta_mat_score, y = delta_NBK*100), size=1.5, color = "#00AFBB") +
  theme_classic(base_size = 16) + coord_cartesian(xlim = c(-.45,0.6), ylim = c(-30,30))+
  labs(x = "Change in rsFC anchored maturation score", y = "Change in N-back accuracy (%)",  title = " ") 
p1
ggsave(file = 'delta_mat_scatter2024_supp1.png',ggMarginal(p1), dpi = 600, width = 7, height = 7)
#############################

# supplementary 2 (exclude low nback performance)
#############################
#############################
dat_both_perfexclude <- dat_both[dat_both$NBK.x > .6 & dat_both$NBK.y > .6,]
model02 <- lmer(delta_NBK ~  1 + delta_mat_score + pFD.x + pFD.y +nruns.x + nruns.y +Male_bin.x + White.x +
                  Black.x + Hispanic.x + Asian.x+ (1|site_id_l.x), data = dat_both_perfexclude)
lm.beta.lmer <- function(mod) {
  b <- fixef(mod)[-1]
  sd.x <- apply(getME(mod,"X")[,-1],2,sd)
  sd.y <- sd(getME(mod,"y"))
  b*sd.x/sd.y
}
lm.beta.lmer(model02)
summary(model02)  # n = 2983, beta = .060
dat_both_perfexclude_complete <- dat_both_perfexclude %>% drop_na(delta_NBK, delta_mat_score, pFD.x, pFD.y, nruns.x, nruns.y, Male_bin.x, White.x, Black.x, Hispanic.x, Asian.x)


model0 <- lmer(NBK.x ~  1 + mat_score.x + pFD.x  +nruns.x +  +Male_bin.x + White.x +
                 Black.x + Hispanic.x + Asian.x+ (1|site_id_l.x), data = dat_both[dat_both$NBK.x > .6,])
lm.beta.lmer(model0)
summary(model0)  # n = 6060, beta = .058
dat_y0_perfexclude_complete <- dat_both[dat_both$NBK.x > .6,] %>% drop_na(NBK.x, mat_score.x, pFD.x, nruns.x, Male_bin.x, White.x, Black.x, Hispanic.x, Asian.x)


model2 <- lmer(NBK.y ~  1 + mat_score.y + pFD.y  +nruns.y +  +Male_bin.y + White.y +
                 Black.y + Hispanic.y + Asian.y+ (1|site_id_l.y), data = dat_both[dat_both$NBK.y > .6,])
lm.beta.lmer(model2)
summary(model2)  # n = 4924, beta = .071
dat_y2_perfexclude_complete <- dat_both[dat_both$NBK.y > .6,] %>% drop_na(NBK.y, mat_score.y, pFD.y, nruns.y, Male_bin.y, White.y, Black.y, Hispanic.y, Asian.y)


t.test(dat_y2_perfexclude_complete$mat_score.y, y = dat_y0_perfexclude_complete$mat_score.x,pair=FALSE) # Welch t = 5.98, p<.001
t.test(dat_both_perfexclude$mat_score.y, y = dat_both_perfexclude$mat_score.x,pair=TRUE) # paired t = 9.91, p<.001
######################### Figure S2

dat1age<- dat1age %>% mutate(dataset="1")
dat2age <- dat2age %>% mutate(dataset="2")

dat1age_temp <- dat1age %>% dplyr::select(NBK, mat_score, dataset)
dat2age_temp <- dat2age %>% dplyr::select(NBK, mat_score, dataset)
dat1exclude_temp <- dat1age_temp[dat1age_temp$NBK > .6,]
dat2exclude_temp <- dat2age_temp[dat2age_temp$NBK > .6,]
dat12exclude <- rbind(dat1exclude_temp, dat2exclude_temp)


p12 <- ggplot(dat12exclude,aes(x = mat_score, y = (NBK)*100/1, color= dataset, group=dataset)) + 
  geom_point(alpha=.20, size = 4) +
  geom_smooth(method = "lm")+
  scale_color_manual(values = c("#E7B800","#c6a2c6"), aesthetics = "color") +
  theme_classic(base_size = 16) + 
  xlim(c(-.5,.6)) +
  ylim(c(0,100)) +
  theme(legend.position = "none") +
  labs(y = "N-back Accuracy (%)", x = 'Anchored rsFC maturation score', color="", title = "") 

p12
ggsave(file = "Fig_marg_supp2.png", ggMarginal(p12, groupColour = TRUE, groupFill = TRUE), width = 7, height = 7)

m1<- mean(dat_both_perfexclude$delta_mat_score[!is.na(dat_both_perfexclude$delta_mat_score) & !is.na(dat_both_perfexclude$delta_NBK)])
m2<- mean(dat_both_perfexclude$delta_NBK[!is.na(dat_both_perfexclude$delta_mat_score) & !is.na(dat_both_perfexclude$delta_NBK)]*100)
mm <- as.data.frame(rbind(c(m1,m2),c(0,0)))
p1 <- ggplot(data=dat_both_perfexclude) + 
  aes(x = delta_mat_score, y = delta_NBK*100) + 
  geom_point(alpha=.06, size = 1.5,color = "#00AFBB") +
  geom_point(data = mm,aes(x = V1,
                           y = V2, size=3),shape = c(4,3), size=c(3,1), color = c("red","black"), inherit.aes = FALSE)+
  geom_smooth(method = "lm", aes(x = delta_mat_score, y = delta_NBK*100), size=1.5, color = "#00AFBB") +
  theme_classic(base_size = 16) + coord_cartesian(xlim = c(-.45,0.6), ylim = c(-30,30))+
  labs(x = "Change in rsFC anchored maturation score", y = "Change in N-back accuracy (%)",  title = " ") 
p1
ggsave(file = 'delta_mat_scatter2024_supp2.png',ggMarginal(p1), dpi = 600, width = 7, height = 7)
