setwd("C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\analysis-nda17-master_1.1\\analysis-nda17-master\\notebooks\\general")
#setwd("/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/analysis-nda17-master_1.1/analysis-nda17-master/notebooks/general/")
#colnames_abcd <- colnames(nda17)
# write.csv(colnames_abcd, file = "colnames_file.csv")
nda17 = readRDS("nda17.Rds")
library(summarytools)
library(dplyr)


#view(dfSummary(nda17$highest_education_imp, max.distinct.values = 25))

ABCD_demographics <- nda17[c("subjectid", "abcd_site","age", "female", "race_ethnicity_imp",
                             "highest_education_imp", "highest_education_cat_imp","highest_income_imp","rel_family_id")]
#should be: SubjID. ScanSite, Age, Male, Female, Race_Ethnicity, Eduation, Income 

#view(dfSummary(ABCD_demographics, max.distinct.values = 25))

cortthick <- nda17[c(9087:9154)]
cortarea <- nda17[c(9301:9368)]
cortvol <- nda17[c(9408:9475)]
subcortvol <- nda17[c(10158:10201)]
rs_corrs_networks <- nda17[c(113:281)]
rs_corrs_networkstosubcort <- nda17[c(18498:18744)]

rs_corrs_networks_uplefttri <- rs_corrs_networks[c(1:13,15:26,29:39,43:52,57:65,71:78,85:91,99:104,113:117,127:130,141:143,155:156,169)]
#residentialvars <- nda17[c(2995:3106)]

colnames(cortthick) <- c("ct_lh_bankssts", "ct_lh_caudalanteriorcingulate", 
                                                     "ct_lh_caudalmiddlefrontal", "ct_lh_cuneus",
                                                     "ct_lh_entorhinal", "ct_lh_fusiform",
                                                     "ct_lh_inferiorparietal", "ct_lh_inferiortemporal",
                                                     "ct_lh_isthmuscingulate","ct_lh_lateraloccipital",
                                                     "ct_lh_lateralorbitofrontal", "ct_lh_lingual",
                                                     "ct_lh_medialorbitofrontal", "ct_lh_middletemporal",
                                                     "ct_lh_parahippocampal", "ct_lh_paracentral",
                                                     "ct_lh_parsopercularis", "ct_lh_parsorbitalis",
                                                     "ct_lh_parstriangularis", "ct_lh_pericalcarine",
                                                     "ct_lh_postcentral", "ct_lh_posteriorcingulate",
                                                     "ct_lh_precentral", "ct_lh_precuneus",
                                                     "ct_lh_rostralanteriorcingulate", "ct_lh_rostralmiddlefrontal",
                                                     "ct_lh_superiorfrontal", "ct_lh_superiorparietal",
                                                     "ct_lh_superiortemporal", "ct_lh_supramarginal",
                                                     "ct_lh_frontalpole", "ct_lh_temporalpole",
                                                     "ct_lh_transversetemporal", "ct_lh_insula",
                                                     "ct_rh_bankssts", "ct_rh_caudalanteriorcingulate", 
                                                     "ct_rh_caudalmiddlefrontal", "ct_rh_cuneus",
                                                     "ct_rh_entorhinal", "ct_rh_fusiform",
                                                     "ct_rh_inferiorparietal", "ct_rh_inferiortemporal",
                                                     "ct_rh_isthmuscingulate","ct_rh_lateraloccipital",
                                                     "ct_rh_lateralorbitofrontal", "ct_rh_lingual",
                                                     "ct_rh_medialorbitofrontal", "ct_rh_middletemporal",
                                                     "ct_rh_parahippocampal", "ct_rh_paracentral",
                                                     "ct_rh_parsopercularis", "ct_rh_parsorbitalis",
                                                     "ct_rh_parstriangularis", "ct_rh_pericalcarine",
                                                     "ct_rh_postcentral", "ct_rh_posteriorcingulate",
                                                     "ct_rh_precentral", "ct_rh_precuneus",
                                                     "ct_rh_rostralanteriorcingulate", "ct_rh_rostralmiddlefrontal",
                                                     "ct_rh_superiorfrontal", "ct_rh_superiorparietal",
                                                     "ct_rh_superiortemporal", "ct_rh_supramarginal",
                                                     "ct_rh_frontalpole", "ct_rh_temporalpole",
                                                     "ct_rh_transversetemporal", "ct_rh_insula")
colnames(cortarea) <- c("ca_lh_bankssts", "ca_lh_caudalanteriorcingulate", 
                                                          "ca_lh_caudalmiddlefrontal", "ca_lh_cuneus",
                                                          "ca_lh_entorhinal", "ca_lh_fusiform",
                                                          "ca_lh_inferiorparietal", "ca_lh_inferiortemporal",
                                                          "ca_lh_isthmuscingulate","ca_lh_lateraloccipital",
                                                          "ca_lh_lateralorbitofrontal", "ca_lh_lingual",
                                                          "ca_lh_medialorbitofrontal", "ca_lh_middletemporal",
                                                          "ca_lh_parahippocampal", "ca_lh_paracentral",
                                                          "ca_lh_parsopercularis", "ca_lh_parsorbitalis",
                                                          "ca_lh_parstriangularis", "ca_lh_pericalcarine",
                                                          "ca_lh_postcentral", "ca_lh_posteriorcingulate",
                                                          "ca_lh_precentral", "ca_lh_precuneus",
                                                          "ca_lh_rostralanteriorcingulate", "ca_lh_rostralmiddlefrontal",
                                                          "ca_lh_superiorfrontal", "ca_lh_superiorparietal",
                                                          "ca_lh_superiortemporal", "ca_lh_supramarginal",
                                                          "ca_lh_frontalpole", "ca_lh_temporalpole",
                                                          "ca_lh_transversetemporal", "ca_lh_insula",
                                                          "ca_rh_bankssts", "ca_rh_caudalanteriorcingulate", 
                                                          "ca_rh_caudalmiddlefrontal", "ca_rh_cuneus",
                                                          "ca_rh_entorhinal", "ca_rh_fusiform",
                                                          "ca_rh_inferiorparietal", "ca_rh_inferiortemporal",
                                                          "ca_rh_isthmuscingulate","ca_rh_lateraloccipital",
                                                          "ca_rh_lateralorbitofrontal", "ca_rh_lingual",
                                                          "ca_rh_medialorbitofrontal", "ca_rh_middletemporal",
                                                          "ca_rh_parahippocampal", "ca_rh_paracentral",
                                                          "ca_rh_parsopercularis", "ca_rh_parsorbitalis",
                                                          "ca_rh_parstriangularis", "ca_rh_pericalcarine",
                                                          "ca_rh_postcentral", "ca_rh_posteriorcingulate",
                                                          "ca_rh_precentral", "ca_rh_precuneus",
                                                          "ca_rh_rostralanteriorcingulate", "ca_rh_rostralmiddlefrontal",
                                                          "ca_rh_superiorfrontal", "ca_rh_superiorparietal",
                                                          "ca_rh_superiortemporal", "ca_rh_supramarginal",
                                                          "ca_rh_frontalpole", "ca_rh_temporalpole",
                                                          "ca_rh_transversetemporal", "ca_rh_insula")
colnames(cortvol) <- c("cv_lh_bankssts", "cv_lh_caudalanteriorcingulate", 
                       "cv_lh_caudalmiddlefrontal", "cv_lh_cuneus",
                       "cv_lh_entorhinal", "cv_lh_fusiform",
                       "cv_lh_inferiorparietal", "cv_lh_inferiortemporal",
                       "cv_lh_isthmuscingulate","cv_lh_lateraloccipital",
                       "cv_lh_lateralorbitofrontal", "cv_lh_lingual",
                       "cv_lh_medialorbitofrontal", "cv_lh_middletemporal",
                       "cv_lh_parahippocampal", "cv_lh_paracentral",
                       "cv_lh_parsopercularis", "cv_lh_parsorbitalis",
                       "cv_lh_parstriangularis", "cv_lh_pericalcarine",
                       "cv_lh_postcentral", "cv_lh_posteriorcingulate",
                       "cv_lh_precentral", "cv_lh_precuneus",
                       "cv_lh_rostralanteriorcingulate", "cv_lh_rostralmiddlefrontal",
                       "cv_lh_superiorfrontal", "cv_lh_superiorparietal",
                       "cv_lh_superiortemporal", "cv_lh_supramarginal",
                       "cv_lh_frontalpole", "cv_lh_temporalpole",
                       "cv_lh_transversetemporal", "cv_lh_insula",
                       "cv_rh_bankssts", "cv_rh_caudalanteriorcingulate", 
                       "cv_rh_caudalmiddlefrontal", "cv_rh_cuneus",
                       "cv_rh_entorhinal", "cv_rh_fusiform",
                       "cv_rh_inferiorparietal", "cv_rh_inferiortemporal",
                       "cv_rh_isthmuscingulate","cv_rh_lateraloccipital",
                       "cv_rh_lateralorbitofrontal", "cv_rh_lingual",
                       "cv_rh_medialorbitofrontal", "cv_rh_middletemporal",
                       "cv_rh_parahippocampal", "cv_rh_paracentral",
                       "cv_rh_parsopercularis", "cv_rh_parsorbitalis",
                       "cv_rh_parstriangularis", "cv_rh_pericalcarine",
                       "cv_rh_postcentral", "cv_rh_posteriorcingulate",
                       "cv_rh_precentral", "cv_rh_precuneus",
                       "cv_rh_rostralanteriorcingulate", "cv_rh_rostralmiddlefrontal",
                       "cv_rh_superiorfrontal", "cv_rh_superiorparietal",
                       "cv_rh_superiortemporal", "cv_rh_supramarginal",
                       "cv_rh_frontalpole", "cv_rh_temporalpole",
                       "cv_rh_transversetemporal", "cv_rh_insula")



#view(dfSummary(ABCD_demographics$abcd_site, max.distinct.values = 25))
#view(dfSummary(ABCD_demographics$rel_family_id, max.distinct.values = 25))


###create dummy variables for demographic variables: (female, race.ethnicity, household.income, highest.education_v2, married)
female_dummy <- model.matrix( ~ female - 1, data= ABCD_demographics)
female <- data.frame(female_dummy)
colnames(female) <- c("Male","Female")
race.ethnicity_dummy <- model.matrix( ~ race_ethnicity_imp - 1, data = ABCD_demographics)
head(race.ethnicity_dummy)
race.ethnicity <- data.frame(race.ethnicity_dummy)
colnames(race.ethnicity) <- c("RE_Hispanic", "RE_White", "RE_Black", "RE_Asian", "RE_Other")

highest.education_dummy <- model.matrix( ~ highest_education_cat_imp - 1, data = ABCD_demographics)
head(highest.education_dummy)
high.edu_cat <- data.frame(highest.education_dummy)
colnames(high.edu_cat) <- c("Edu-LessHighHS","Edu-HSDip","Edu-SomeCollege", "Edu-Bachelor","Edu-PostGrad")

household.income_dummy <- model.matrix( ~ highest_income_imp - 1, data = ABCD_demographics)
head(household.income_dummy)
household.income <- data.frame(household.income_dummy)
colnames(household.income) <- c("Inc-Greater100K","Inc-Less50K","Inc-50to100K")

scannersite_dummy <- model.matrix( ~ abcd_site - 1, data = ABCD_demographics)
scannersite <- data.frame(scannersite_dummy)


ABCD_Demo <- cbind.data.frame(ABCD_demographics$subjectid, scannersite, ABCD_demographics$age, female, 
                              race.ethnicity, high.edu_cat, ABCD_demographics$highest_education_imp, household.income, 
                              ABCD_demographics$rel_family_id)
colnames(ABCD_Demo)[c(1:23, 36,40)] <- c("SubjID", "ScanSite_01", "ScanSite_02", "ScanSite_03", "ScanSite_04", "ScanSite_05",
                                      "ScanSite_06", "ScanSite_08", "ScanSite_09", "ScanSite_10", "ScanSite_11", "ScanSite_12",
                                      "ScanSite_13", "ScanSite_14", "ScanSite_15", "ScanSite_16", "ScanSite_17", "ScanSite_18",
                                      "ScanSite_19", "ScanSite_20", "ScanSite_21", "ScanSite_22", "Age", "HighEdu_Num", "FamilyID")

ABCD_Brain_Struct <- cbind.data.frame(cortthick,cortarea,cortvol,subcortvol)

#ct idx 1-68, ca idx 69-132, cv idx 137-204, subcortv idx 205-248

#setwd("C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Preprocessing\\Release1.1\\")
setwd("/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Preprocessing/Release1.1/")
ABCD_BrainStruct_2 <- mutate_all(ABCD_Brain_Struct, function(x) as.numeric(as.character(x)))
write.csv(ABCD_Demo, "ABCD_Demo.csv")
write.csv(ABCD_BrainStruct_2, "ABCD_Brain_Struct.csv", row.names = FALSE, na = "0")
#view(dfSummary(ABCD_BrainStruct_2))


attach(nda17)
 #remember to detach
Youth_Ratings <- nda17[,c("src_subject_id", "pmq_y_ss_mean", "fes_y_ss_fc_pr", "psb_y_ss_mean", 
                           "crpbi_y_ss_parent", "crpbi_y_ss_caregiver", 
                           "srpf_y_ss_ses", "srpf_y_ss_iiss", "srpf_y_ss_dfs", "physical_activity1_y")]
Parent_Ratings <- nda17[,c("src_subject_id", "meim_p_ss_exp", "meim_p_ss_com", "meim_p_ss_total", 
                          #"via_p_ss_hc", "via_p_ss_amer", 
                          "nsc_p_ss_mean_3_items", "fes_p_ss_fc_pr",
                          "macv_p_ss_fs", "macv_p_ss_fo", "macv_p_ss_isr", "macv_p_ss_fr", "macv_p_ss_r",
                          "psb_p_ss_mean")]
Neighborhood_Variables <- nda17[,c("src_subject_id", "reshist_addr1_d1a", "reshist_addr1_walkindex",
                                  "reshist_addr1_grndtot", "reshist_addr1_p1tot", "reshist_addr1_p1vlnt",
                                  "reshist_addr1_drugtot", "reshist_addr1_drgsale", "reshist_addr1_drgposs", 
                                  "reshist_addr1_mjsale", "reshist_addr1_dui",
                                  "reshist_addr1_adi_edu_l", "reshist_addr1_adi_edu_h", "reshist_addr1_adi_work_c",
                                  "reshist_addr1_adi_income", "reshist_addr1_adi_in_dis", "reshist_addr1_adi_home_v",
                                  "reshist_addr1_adi_rent", "reshist_addr1_adi_mortg", "reshist_addr1_adi_home_o",
                                  "reshist_addr1_adi_crowd", "reshist_addr1_adi_unemp", "reshist_addr1_adi_pov",
                                  "reshist_addr1_adi_b138", "reshist_addr1_adi_sp", "reshist_addr1_adi_ncar",
                                  "reshist_addr1_adi_ntel", "reshist_addr1_adi_nplumb", "reshist_addr1_adi_wsum",
                                  "reshist_addr1_adi_perc", "reshist_addr1_popdensity", "reshist_addr1_no2",
                                  "reshist_addr1_pm25", "reshist_addr1_proxrd")]

Cog_Scores <- nda17[,c("nihtbx_picvocab_agecorrected", "nihtbx_flanker_agecorrected",
                       "nihtbx_list_agecorrected", "nihtbx_cardsort_agecorrected",
                       "nihtbx_pattern_agecorrected", "nihtbx_picture_agecorrected",
                       "nihtbx_reading_agecorrected", "nihtbx_fluidcomp_agecorrected",
                       "nihtbx_cryst_agecorrected", "nihtbx_totalcomp_agecorrected")]
                       #"cash_choice_task", 
                       #"lmt_scr_efficiency", "pea_wiscv_tss")]
view(dfSummary(All_Data_clean[c(22:50)], max.distinct.values = 25))                       
detach(nda17)

#ABCD_Brain_Data <- read.csv("ABCD_Brain_Struct.csv")
#ABCD_Demo <- read.csv("ABCD_Demo.csv")
#ABCD_Brain_Data[ABCD_Brain_Data == 0] = NA

All_Data_1 <- cbind.data.frame(ABCD_Demo, ABCD_Brain_Struct, Youth_Ratings, Parent_Ratings, Neighborhood_Variables, Cog_Scores)
All_Data_clean <- All_Data_1[complete.cases(All_Data_1),] #3150 cases

All_Data_2 <- cbind.data.frame(ABCD_Demo, rs_corrs_networks_uplefttri, Youth_Ratings, Parent_Ratings, Neighborhood_Variables, Cog_Scores)
All_Data_clean_rsfmri <- All_Data_2[complete.cases(All_Data_2),]

test1 <- ABCD_Demo[complete.cases(ABCD_Demo),] #4521
test2 <- ABCD_Brain_Struct[complete.cases(ABCD_Brain_Struct),] #4258
test3 <- Youth_Ratings[complete.cases(Youth_Ratings),] #4186
test4 <- Parent_Ratings[complete.cases(Parent_Ratings),]#4233 
test5 <- Neighborhood_Varibles[complete.cases(Neighborhood_Variables),] #4037
test6 <- Cog_Scores[complete.cases(Cog_Scores),]#4296

#write.csv(All_Data_clean, "ABCD_Preprocessed_BrainDemoParentYouthEnvCog_NoMissing.csv", row.names = FALSE)
write.csv(All_Data_clean_rsfmri, "ABCD_Preprocesses_rsfmriDemo_NoMissing.csv", row.names = FALSE)

### THis section on is to log transform relevant variables
library(psych)
library(ggplot2)
library(lme4)
library(stargazer)
library(gvlma)
library(extrafont)
library(lattice)
library(summarytools)
library(dplyr)
library(ade4)
library(factoextra)
library(magrittr)
library(FactoMineR)

ABCD_Data <- read.csv("C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Preprocessing\\Release1.1\\ABCD_Preprocessed_BrainDemoParentYouthEnvCog_NoMissing.csv")
XScores <- read.csv("C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Analysis\\CanCorrs\\Release1.1\\PrepForMediation\\XScores_v4_Visualization.csv")
YScores <- read.csv("C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Analysis\\CanCorrs\\Release1.1\\PrepForMediation\\YScores_v4_Visualization.csv")
setwd("C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Preprocessing\\Release1.1\\")
hist(ABCD_Data$reshist_addr1_d1a)
Res_density_log <- log1p(ABCD_Data$reshist_addr1_d1a)


UCR_GrandTotal_log <- log1p(ABCD_Data$reshist_addr1_grndtot)
UCR_TotOffenses_log <- log1p(ABCD_Data$reshist_addr1_p1tot)
UCR_ViolCrime_log <- log1p(ABCD_Data$reshist_addr1_p1vlnt)
UCR_DrugAbuse_log <- log1p(ABCD_Data$reshist_addr1_drugtot)
UCR_DrugSale_log <- log1p(ABCD_Data$reshist_addr1_drgsale)
UCR_DrugPoss_log <- log1p(ABCD_Data$reshist_addr1_drgposs)
UCR_DUI_log <- log1p(ABCD_Data$reshist_addr1_dui)

UCRCrime <- data.frame(UCR_GrandTotal_log, UCR_TotOffenses_log, UCR_ViolCrime_log, 
                       UCR_DrugAbuse_log, UCR_DrugSale_log, UCR_DrugPoss_log, UCR_DUI_log)

hist(UCR_GrandTotal_log)
hist(UCR_TotOffenses_log)
hist(UCR_ViolCrime_log)
hist(UCR_DrugAbuse_log)
hist(UCR_DrugSale_log)
hist(UCR_DrugPoss_log)

ProxRoads_log <- log(ABCD_Data$reshist_addr1_proxrd)
hist(ProxRoads_log)

res.pca_UCR <- PCA(UCRCrime, scale.unit = TRUE, graph = T)
fviz_screeplot(res.pca_UCR, addlabels =TRUE)
fviz_pca_var(res.pca_UCR, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
fviz_contrib(res.pca_UCR, choice = "var", axes = 1, top = 20)
fviz_contrib(res.pca_UCR, choice = "var", axes = 2, top = 20)

UCR_Crime_PC <- res.pca_UCR$ind$coord[,1]

ADI_Vars <- ABCD_Data[c(322:340)]
ADI_Vars[paste0(names(ADI_Vars), '_log')] <- lapply(ADI_Vars, function(x) log1p(x))

ADI_Vars_Corrected <- ADI_Vars[c(20, 2:9, 29:36)]

res.pca_ADI <- PCA(ADI_Vars_Corrected, scale.unit = TRUE, graph = T)
fviz_screeplot(res.pca_ADI, addlabels =TRUE)
fviz_pca_var(res.pca_ADI, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), 
             repel = TRUE # Avoid text overlapping
)
fviz_contrib(res.pca_ADI, choice = "var", axes = 1, top = 20)
fviz_contrib(res.pca_ADI, choice = "var", axes = 2, top = 20)

ADI_PCs <- res.pca_ADI$ind$coord[,c(1,2)]
colnames(ADI_PCs) <- c("ADI_PC1", "ADI_PC2")

cor(ABCD_Data$reshist_addr1_adi_perc, res.pca_ADI$ind$coord[,1])
plot(ABCD_Data$reshist_addr1_adi_perc, res.pca_ADI$ind$coord[,1])

Log_Vars <- data.frame(Res_density_log,ProxRoads_log, UCR_Crime_PC, ADI_PCs)

ABCD_Data_new <- cbind.data.frame(ABCD_Data, Log_Vars)
write.csv(ABCD_Data_new, "ABCD_Preprocessed_NoMissing_LogFixed.csv", row.names = FALSE)

#following two lines are just to check if there are any NAs in the data
indx <- apply(ABCD_Data_new, 2, function(x) any(is.na(x) ))
colnames(indx)

