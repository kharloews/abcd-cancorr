%% Script file to run multiple regressions on ABCD Data
%% Change script to fit the OS you are using
OS_for_analysis = 2 % 1 for Mac, 2 for Windows
if OS_for_analysis == 1
    preproc_path = '/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Preprocessing/'
    filepath_to_save = '/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Preprocessing/Residuals/'
else preproc_path = 'C:\Users\Berman Lab\Dropbox\OldLaptop\DropBox_Documents\ABCD\Preprocessing\Release1.1\'
    filepath_to_save = 'C:\Users\Berman Lab\Dropbox\OldLaptop\DropBox_Documents\ABCD\Preprocessing\Release1.1\Residuals\'
end
%% Load in data preprocessed using R scripts
All_Data = readtable([preproc_path 'ABCD_Preprocessed_BrainDemoParentYouthEnvCog_NoMissing.csv']);
ABCD_Demos = All_Data(:, 2:38);
ABCD_BrainStructs = All_Data(:,41:288) ;
file_prefix = 'ABCD_lm_'

%% Subset Brain Data
ABCD_subs = ABCD_Demos(:,1);
DVs_Thick = ABCD_BrainStructs(:, 1:68);
DVs_Area = ABCD_BrainStructs(:, 69:136);
DVs_Vol = ABCD_BrainStructs(:, 137:204);
DVs_SubcortVol = ABCD_BrainStructs(:, 205:248);
%% Create IVs -- Factors in the regression
% This will set IVs to be Scanner Site, number 20 is not included and is
% the reference
IVs_ScannerSite = ABCD_Demos(:,1:20);

%% This runs the models for the individual areas, and also for DK regions      
for m = [1];
    Type_of_Model = m
    All_DVs = {DVs_Thick, DVs_Area, DVs_Vol, DVs_SubcortVol};
if Type_of_Model == 2
        DV_Model_names = {'CortThick','CortArea', 'CortVol', 'SubcortVol'};
        IVs = IVs_Default;        
else
       DV_Model_names = {'CortThick_SS', 'CortArea_SS','CortVol_SS', 'SubcortVol_SS'};
       IVs = IVs_ScannerSite;
end 
%for q = [3,4,5]
for q = 1:length(All_DVs);
    
Coeffs = [];
SEs = [];
tstats = [];
pvalues = [];
Residuals_Raw = [];
Residuals_Norm = [];
Fitted_Values = [];
Model_Specifics = [];
Model_Specifics_names = {'AIC'; 'AICc'; 'BIC'; 'CAIC'; 'MSE'; 'Rsquared_ord'; 'Rsquared_adj'};
DVs = All_DVs{q};
DVs_regions = DVs.Properties.VariableNames;
for  s = 1:width(DVs); 
    
    DV_of_interest = DVs(:,s);
    Model = [IVs DV_of_interest];
    lm = fitlm(Model);
    coeffs_temp = lm.Coefficients.Estimate;
    SEs_temp = lm.Coefficients.SE;
    tstats_temp = lm.Coefficients.tStat;
    pvalues_temp = lm.Coefficients.pValue;
    res_raw_temp = lm.Residuals.Raw;
    res_norm_temp = lm.Residuals.Standardized;
    fitted_temp = lm.Fitted;
    AIC_temp = lm.ModelCriterion.AIC;
    AICc_temp = lm.ModelCriterion.AICc;
    BIC_temp = lm.ModelCriterion.BIC;
    CAIC_temp = lm.ModelCriterion.CAIC;
    MSE_temp = lm.MSE;
    Rsquared_ord_temp = lm.Rsquared.Ordinary;
    Rsquared_adj_temp = lm.Rsquared.Adjusted;
    model_specifics_temp = vertcat(AIC_temp, AICc_temp, BIC_temp, CAIC_temp, MSE_temp,Rsquared_ord_temp,Rsquared_adj_temp);
    
    Coeffs = [Coeffs coeffs_temp];
    pvalues = [pvalues pvalues_temp];
    SEs = [SEs SEs_temp];
    tstats = [tstats tstats_temp];
    Residuals_Raw = [Residuals_Raw res_raw_temp];
    Residuals_Norm = [Residuals_Norm res_norm_temp];
    Standard_Betas = Coeffs ./ SEs;
    Fitted_Values = [Fitted_Values fitted_temp];
    Model_Specifics = [Model_Specifics model_specifics_temp];
    
end
s
Coeff_Names = transpose(lm.CoefficientNames);
save([filepath_to_save file_prefix DV_Model_names{q} '.mat'], 'Coeffs',...
    'pvalues', 'SEs', 'tstats', 'Residuals_Raw', 'Residuals_Norm', 'Coeff_Names', 'ABCD_subs',...
    'Standard_Betas', 'Fitted_Values','Model_Specifics','Model_Specifics_names', 'DVs_regions') 
end 
end







