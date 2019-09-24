clear all
%% Script to run Canonical Correlations between Structurals and Demographics 
%% Read in relevant files
OS_for_analysis = 2 % 1 for Mac, 2 for Windows
if OS_for_analysis == 1
    load('/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Preprocessing/Release1.1/Residuals/ABCD_lm_CortVol_SS.mat')
    ABCD_Clean = readtable('/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Preprocessing/Release1.1/ABCD_Preprocessed_BrainDemoParentYouthEnvCog_NoMissing.csv');
    addpath('/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Analysis/freesurfer_statsurf_display/')
    filepath_to_save = '/Users/_omsa-tech/Dropbox/OldLaptop/DropBox_Documents/ABCD/Analysis/CanCorrs/Release1.1/CanCorrResults/'
else load('C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Preprocessing\\Release1.1\\Residuals\\ABCD_lm_CortVol_SS.mat')
    ABCD_Clean = readtable('C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Preprocessing\\Release1.1\\ABCD_Preprocessed_NoMissing_LogFixed.csv');
    addpath('C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Analysis\\freesurfer_statsurf_display\\')
    filepath_to_save = 'C:\\Users\\Berman Lab\\Dropbox\\OldLaptop\\DropBox_Documents\\ABCD\\Analysis\\CanCorrs\\Release1.1\\CanCorrResults\\'

end


%% Set up Variables
Indices_by_Lobe = [3,11,13,16:19,23,26,27,31, 37, 45, 47, 50:53,57, 60,61, 65,...
    7,21,24,28,30,41, 55,58,62,64,...
    1,5,6,8,14,15,29,32,33, 35,39,40,42,48,49,63,66,67,...
    4,10,12,20,38,44,46,54,...
    2,9,22,25, 36,43,56,59,...
    34, 68]%This is to go from regular DK order to the one in the CanCorrs

To_Orig_Indices = [33, 59,1,51,34,35,23,36,60,52,2,53,3,37,38,...
    4,5,6,7,54,24,61,8,25,62,9,10,26,39,27,11,40,41,67,...
    42,63,12,55,43,44,28,45,64,56,13,57,14,46,47,...
    15,16,17,18,58,29,65,19,30,66,20,21,31,48,32,22,49,50, 68]%This is to go from CanCorrs to the Original DK order
ABCD_X_VOL = Residuals_Norm(:,Indices_by_Lobe);

ABCD_Demos = ABCD_Clean(:, [23,25, 27, 26,28:29,31:33,35,38,37]); % Order has been changed to reflect: Age, Female, White, Hispanic, Black, Asian, Edu (minus Bachelor, <50K, >100K
Parent_Vars = ABCD_Clean(:, [300,301,303:310]);
Youth_Vars = ABCD_Clean(:, 290:298);
Neighbor_Vars = ABCD_Clean(:, [313, 358, 359, 357, 340, 355, 342, 343, 356]);
Cog_Vars = ABCD_Clean(:, [345:351]);

ABCD_Y_VOL = table2array([ABCD_Demos Parent_Vars Youth_Vars Neighbor_Vars]);

data = horzcat(ABCD_X_VOL, ABCD_Y_VOL);
fileprefix = 'ABCD_CanCorr_CortVol_v5'
figuretitle = 'Canonical Correlation ABCD Cort Vol Demo YPNeigh Cog-- '
%% This will create a folder to save all the figures
if exist([filepath_to_save fileprefix])==0
    mkdir(filepath_to_save, fileprefix)
end 
%%
VarNames_temp = {'lh bankssts';'lh caudalanteriorcing';'lh caudalmiddlefront';...
'lh cuneus'; 'lh entorhinal'; 'lh fusiform';...
'lh inferiorpariet';'lh inferiortemp';'lh isthmuscing';...        
'lh ltr occipital'; 'lh ltr orbitofront'; 'lh lingual';...
'lh medialorbitofront'; 'lh middletemp'; 'lh parahippocampal';... 
'lh paracentral'; 'lh parsopercularis'; 'lh parsorbitalis';...
'lh parstriangularis'; 'lh pericalcarine'; 'lh postcentral';...
'lh posteriorcing'; 'lh precentral'; 'lh precuneus';...
'lh rostralanteriorcing'; 'lh rostralmiddlefront'; 'lh superiorfront';...
'lh superiorpariet'; 'lh superiortemp'; 'lh supramarginal';... 
'lh frontpole'; 'lh temppole'; 'lh transversetemp';'lh insula';...
'rh bankssts'; 'rh caudalanteriorcing'; 'rh caudalmiddlefront';...
'rh cuneus'; 'rh entorhinal'; 'rh fusiform';...
'rh inferiorpariet'; 'rh inferiortemp'; 'rh isthmuscing';...
'rh ltr occipital'; 'rh ltr orbitofront'; 'rh lingual';...
'rh medialorbitofront'; 'rh middletemp'; 'rh parahippocampal';...
'rh paracentral'; 'rh parsopercularis'; 'rh parsorbitalis';...
'rh parstriangularis'; 'rh pericalcarine'; 'rh postcentral';...
'rh posteriorcing'; 'rh precentral'; 'rh precuneus';...
'rh rostralanteriorcing'; 'rh rostralmiddlefront'; 'rh superiorfront';...
'rh superiorpariet'; 'rh superiortemp'; 'rh supramarginal';...
'rh frontpole'; 'rh temppole'; 'rh transversetemp';'rh insula';...
'Age'; 'Female';... % Demo
'White'; 'Hispanic'; 'Black'; 'Asian';... % Race
'LT HS Edu'; 'HS Diploma Edu'; 'Some College Edu';'Post Grad Edu';...%SES
'LT 50k H Inc'; '100k H Inc';
'MEIM Exploration P'; 'MEIM Comm & Attach P';...%Parent
'Neigh Safety Prot P'; 'Fam Env Conflict P'; 'MACV Fam Support P';...
'MACV Fam Oblig P'; 'MACV Indep SelfReliance P'; 'MACV Family as Ref P';...
'MACV Religion P'; 'ProSocial Beh P';...
'Parental Monitoring Quest Y'; 'Fam Environ Scale Conflict Y'; 'ProSocial Beh Y';...%Youth
'CRPBI Accept Parent Y'; 'CRPBI Accept CG T'; 'SRPF School Env Y';...
'SRPF School Involv Y'; 'SRPF School Diseng Y';'Physical Activity';...
'Walkability Idx'; 'ADIPC1'; 'ADIPC2'; 'Crime_PC';... %Neigh Est
'ADI Percentile';...
'Pop Density'; 'NO2'; 'PM 2.5'; 'Prox to Roads'};... %Satellite measure

XVariableNames = 'Desikan-Killani Cort Areas'
YVariableNames = 'Subject Reported Data'
VarNames = VarNames_temp([3,11,13,16:19,23,26,27,31, 37, 45, 47, 50:53,57, 60,61, 65,...
    7,21,24,28,30,41, 55,58,62,64,...
    1,5,6,8,14,15,29,32,33, 35,39,40,42,48,49,63,66,67,...
    4,10,12,20,38,44,46,54,...
    2,9,22,25, 36,43,56,59,...
    34, 68, 69:108]);

xindices = 1:68;
yindices = 69:108;
T = data(:,xindices); 
H = data(:,yindices); 
XX=zscore(T);
YY=zscore(H);
[Xcoef Ycoef r1 Xscore Yscore stats]=canoncorr(XX,YY);
[Xcoef_orig Ycoef_orig r1_orig Xscore_orig Yscore_orig stats_orig]=canoncorr(XX,YY);
X_loading = corr(XX,Xscore);
Y_loading = corr(YY,Yscore);

mean_loads_all_X = [];
mean_loads_all_Y = [];
mean_loads_all_X_vis = [];
mean_loads_all_Y_vis = [];

X_sig_loads_all = [];
X_loads_idx = [];

Y_sig_loads_all = [];
Y_loads_idx = [];


X_times_loadings_LVs = {};
Y_times_loadings_LVs = {};

X_LVs_names = {};
Y_LVs_names = {};

cont_for_loop = 1
orig_u = Xscore;
orig_v = Yscore;

for k = 1:min(size(T,2),size(H,2)); %count LVs
    k
    if cont_for_loop==0
        break 
    end
    
    allXloads = [];
    allYloads = [];
    pee=[];
    
    rr =[];
    meanCorrsX = X_loading;
    meanCorrsY = Y_loading;
    X = XX;
    Y = YY;
    numsamp = 4000;
        sampsize = size(T,1);
    for z=1:numsamp
        tempd=datasample([X Y], sampsize);
        Xs=tempd(:,1:size(T,2));
        Ys=tempd(:,size(T,2)+1:(size(T,2)+size(H,2)));
        [Xcoef Ycoef r Xscore Yscore stats]=canoncorr(Xs,Ys);

%         X_loading = corr(Xs,Xscore);
%         Y_loading = corr(Ys,Yscore);
%         
%         
%         A = X_loading(:,k).*sign(meanCorrsX(:,k));
%         B = Y_loading(:,k).*sign(meanCorrsY(:,k));
%         for p = 1:size(A)
%             if A(p)>0
%                 A(p) = X_loading(p,k);
%             end
%         end
%         for p =1:size(B)
%             if B(p)>0
%                 B(p) = Y_loading(p,k);
%             end
%         end
%         allXloads = [allXloads;A'];
%         allYloads = [allYloads;B'];
%         rr = [rr;r(k)];
%     end
                pee=[pee;stats.p];
                
                Xscore = (Xs - repmat(mean(Xs),sampsize,1))*Xcoef;
                Yscore = (Ys - repmat(mean(Ys),sampsize,1))*Ycoef;
                temp_u = (X - repmat(mean(X),sampsize,1))*Xcoef;
                temp_v = (Y - repmat(mean(Y),sampsize,1))*Ycoef;
                
                X_loading = Xcoef;
                Y_loading = Ycoef;
                
                X_loading = corr(Xs,Xscore);%**********************
                Y_loading = corr(Ys,Yscore);
                
                [rss,validityx] = corr(orig_v,temp_v); [rss2,validityy] = corr(orig_u,temp_u);
                %         validcolumnx = ceil(abs(validityx(k,1))-0.71); % or threshold set as minimum r that results in p <0.05 for df?
                %         validcolumny = ceil(abs(validityy(k,1))-0.71); % is this valid for flipping
                validcolumnx = validityx(k,k); % or threshold set as minimum r that results in p <0.05 for df?
                validcolumny = validityy(k,k);
                %%%%%$###
%                 if (validcolumnx<0.05 & validcolumny<0.05 & sign(rss(k,k))<0 & sign(rss2(k,k))<0)
%                     
%                     A = -X_loading(:,k); B = -Y_loading(:,k);
%                 else
%                     A = X_loading(:,k); B = Y_loading(:,k);
%                 end
                %%%%%%%$##
                    if sign(rss(k,k))<0 
                    
                    A = -X_loading(:,k); 
                else
                    A = X_loading(:,k);
                    end
                    if sign(rss2(k,k))<0
                    B = -Y_loading(:,k);
                    else
                     B = Y_loading(:,k);  
                    end
                    %               A = X_loading(:,k).*sign(meanCorrsX(:,k));
                    %               B = Y_loading(:,k).*sign(meanCorrsY(:,k));
                    %               %         pmax1 = find(abs(meanCorrsX(:,k))==max(abs(meanCorrsX(:,k))));
                    %               %         pmax2 = find(abs(meanCorrsY(:,k))==max(abs(meanCorrsY(:,k))));
                    %
                    %               for p = 1:size(A)
                    %                   if A(p)>0 %instead of p use a fixed index for the largest load
                    %                       A(p) = X_loading(p,k);
                    %                   end
                    %               end
                    %
                    %
                    %               for p =1:size(B)
                    %                   if B(p)>0 %instead of p use a fixed index for the largest load
                    %                       B(p) = Y_loading(p,k);
                    %                   end
                    %               end
                
                allXloads = [allXloads;A'];
                allYloads = [allYloads;B'];
                rr = [rr;r(k)];
    end
            
        %Start Omid
    loc1 = [1:size(allXloads,2)]';
    col1 = [1*ones(11,1);2*ones(11,1);3*ones(5,1);4*ones(5,1);5*ones(9,1);...
        6*ones(9,1);7*ones(4,1);8*ones(4,1);9*ones(4,1);10*ones(4,1); 11*ones(1,1); 12*ones(1,1)];
    sigcol1 = sign([(mean(allXloads)-2*std(allXloads)).*(mean(allXloads)+2*std(allXloads))]');
    Huh1 = [mean(allXloads)' loc1 col1 sigcol1];
    Hmm1 = Huh1(:,1);
    
    loc2 = [1:size(allYloads,2)]';
    col2 = [1*ones(6,1); 2*ones(6,1); 3*ones(10,1); 4*ones(9,1);...
        5*ones(9,1)];
    sigcol2 = sign([(mean(allYloads)-2*std(allYloads)).*(mean(allYloads)+2*std(allYloads))]');
    Huh2 = [mean(allYloads)' loc2 col2 sigcol2];
    Hmm2 = Huh2(:,1);
    % End Omid
    
    figure(k);
    set(gcf,'color','w');
 subplot(1,2,1);
    h = barh(mean(allXloads),.5,'FaceColor',[0 .6 .6],'EdgeColor',[0.1 .9 .9],'LineWidth',1.5);
    hold on
    mean_allXloads = mean(allXloads);
    SE_allXloads = 2*std(allXloads);
    for i = 1:numel(Hmm1)
        
        h = barh(i,Hmm1(i));
        if i == 1, hold on, end
        he = errorbar(Hmm1(i),i, 2*std(allXloads(:,i)),'horizontal','.r');
        if Huh1(i, 3) == 1
            col = [255 219 15] ./ 255;
        elseif Huh1(i, 3) == 2
            col = [255 236 131] ./ 255;
        elseif Huh1(i, 3) == 3
            col = [255 15 54] ./ 255;
        elseif Huh1(i, 3) == 4
            col = [255 131 154] ./ 255;
        elseif Huh1(i, 3) == 5
            col = [232 116 17] ./ 255;
        elseif Huh1(i, 3) == 6
            col = [232 154 88] ./ 255;
        elseif Huh1(i, 3) == 7
            col = [119 2 232] ./ 255;
        elseif Huh1(i, 3) == 8
            col = [170 107 232] ./ 255;
        elseif Huh1(i, 3) == 9
            col = [3 110 255] ./ 255;
        elseif Huh1(i, 3) == 10
            col = [118 178 255] ./ 255;
        elseif Huh1(i, 3) == 11
            col = [62 255 91] ./ 255;
        elseif Huh1(i, 3) == 12
            col = [122 255 148] ./ 255;
        else col = 'c';
        end
         if Huh1(i,4)<=0 
             col = 'w';
         end
        set(h, 'FaceColor', col)

    
end
    set(gca,'FontWeight','bold','FontSize',16)
    title(XVariableNames);
    hold on
    %errorbar(mean(allXloads),2*std(allXloads),'.r', 'horizontal');
    set(gca,'Ytick',[1:length(VarNames(xindices))]);
    set(gca,'Yticklabel',VarNames(xindices),'FontWeight','bold','FontSize',12)
    set(gca, 'YAxisLocation', 'left')
    xlabel('Canonical Loading','FontWeight','bold','FontSize',16);
    hold off
 subplot(1,2,2);
     barh(mean(allYloads),.5,'FaceColor',[0 .6 .6],'EdgeColor',[0.1 .9 .9],'LineWidth',1.5);
    %%
    hold on
    mean_allYloads = mean(allYloads);
    SE_allYloads = 2*std(allYloads);
for i = 1:numel(Hmm2)
        
        h = barh(i,Hmm2(i));
        if i == 1, hold on, end
        he = errorbar(Hmm2(i),i, 2*std(allYloads(:,i)),'horizontal','.r');
        if Huh2(i, 3) == 1
            col = [230 41 73] ./ 255;
        elseif Huh2(i, 3) == 2
            col = [240 118 36] ./ 255;
        elseif Huh2(i, 3) == 3
            col = [244 189 45] ./ 255;
        elseif Huh2(i, 3) == 4
            col = [30 212 222] ./ 255;
        elseif Huh2(i, 3) == 5
            col = [28 125 225] ./ 255;
        elseif Huh2(i, 3) == 6
            col = [232 154 88] ./ 255;
        elseif Huh2(i, 3) == 7
            col = [119 2 232] ./ 255;
        elseif Huh2(i, 3) == 8
            col = [170 107 232] ./ 255;
        elseif Huh2(i, 3) == 9
            col = [3 110 255] ./ 255;
        elseif Huh2(i, 3) == 10
            col = [118 178 255] ./ 255;
        elseif Huh2(i, 3) == 11
            col = [62 255 91] ./ 255;
        elseif Huh2(i, 3) == 12
            col = [122 255 148] ./ 255;
        else col = 'c';
        end
         if Huh2(i,4)<=0 
             col = 'w';
         end
        set(h, 'FaceColor', col)
                     
    
    end

hold off
%%
set(gca,'FontWeight','bold','FontSize',16)
    title(YVariableNames);% rr or r1?
    hold on
    %errorbar(mean(allYloads),2*std(allYloads),'horizontal','.r');
    set(gca,'Ytick',[1:length(VarNames(yindices))]);
    set(gca,'Yticklabel',VarNames(yindices),'FontWeight','bold','FontSize',12)
    set(gca, 'YAxisLocation', 'right')
    xlabel('Canonical Loading','FontWeight','bold','FontSize',16);
    hold off
    sgtitle(['Latent Variable ', num2str(k), ', F (', num2str(stats.df1(k)), ', ', num2str(stats.df2(k)), ')= ',num2str(stats.F(1)), ', R^2 = ',num2str(mean(rr).^2),', p< ', num2str(stats.pF(k)),' R = ',num2str(mean(rr))]);
    %% section to prep values for visualization
    mean_loads_all_X(:,k) = mean_allXloads;
    mean_loads_all_Y(:,k)  = mean_allYloads;
    X_sig_loads = 0;
    X_sig_loads_count = 0;
    X_sig_idx = [];
    
    Y_sig_loads = 0;
    Y_sig_loads_count = 0;
    Y_sig_idx = [];
%% X var 
    for j = 1:numel(mean_allXloads); %This replaces meanXloads by 0 so they aren't shown on visualization
        if Huh1(j,4) <= 0;
            mean_loads_all_X_vis(j,k) = 0;
        else mean_loads_all_X_vis(j,k) = mean_loads_all_X(j,k);
        end
        if abs(mean_loads_all_X_vis(j,k))>0; %This separates the positive and negative non-zero loadings 
            X_sig_loads = X_sig_loads + mean_loads_all_X_vis(j,k);
            X_sig_loads_count = X_sig_loads_count + 1;
            X_sig_idx = [X_sig_idx, j];
        end     
    end
    
    if X_sig_loads_count > 1
            X_sig_loads_mean = X_sig_loads ./ X_sig_loads_count
        else X_sig_loads_mean = 0
    end
%% Y var
        for j = 1:numel(mean_allYloads); %This replaces meanXloads by 0 so they aren't shown on visualization
        if Huh2(j,4) <= 0;
            mean_loads_all_Y_vis(j,k) = 0;
        else mean_loads_all_Y_vis(j,k) = mean_loads_all_Y(j,k);
        end
        if abs(mean_loads_all_Y_vis(j,k))>0; %This separates the positive and negative non-zero loadings 
            Y_sig_loads = Y_sig_loads + mean_loads_all_Y_vis(j,k);
            Y_sig_loads_count = Y_sig_loads_count + 1;
            Y_sig_idx = [Y_sig_idx, j];
        end     
    end
    
    if Y_sig_loads_count > 1
            Y_sig_loads_mean = X_sig_loads ./ Y_sig_loads_count
        else Y_sig_loads_mean = 0
    end

    X_times_loadings_LVs{1,k} = zscore(ABCD_X_VOL(:,X_sig_idx)) .* mean_allXloads(1,X_sig_idx);
    Y_times_loadings_LVs{1,k} = zscore(ABCD_Y_VOL(:,Y_sig_idx)) .* mean_allYloads(1,Y_sig_idx);
    
    X_sig_loads_idx{1,k} = X_sig_idx;
    Y_sig_loads_idx{1,k} = Y_sig_idx;
    
    X_LVs_sig_names{1,k} = VarNames(X_sig_idx,:);
    Y_LVs_sig_names{1,k} = VarNames(Y_sig_idx + 68,:);
    
    mean_loads_all_X_vis_orig = (mean_loads_all_X_vis(To_Orig_Indices,:));
    if sum(mean_loads_all_X_vis(:,k))==0 
        cont_for_loop = 0 %This is the counter that appears at the beginning of the for loop, that will stop you from running
        % more LV plots than necessary
    
    end
    
    %% Code to save LV bar plots figures
baseFileName = sprintf([filepath_to_save fileprefix '/' fileprefix '_LV%d.fig'],k)
savefig(baseFileName)
%% Code to visualize brain areas
Vis_Matrix =[];
Vis_Matrix(:,1) = [mean_loads_all_X_vis_orig(1:34,k)];
Vis_Matrix(:,2) = [mean_loads_all_X_vis_orig(35:68,k)];

Values = cell(1, 2);
ValuesMask = cell(1, 2);

for z = 1:2
	Values{z} = Vis_Matrix(:,z);
    for j=1:numel(Vis_Matrix(:,z))
        Vis_Matrix(j,z)
        if Vis_Matrix(j,z)==0
            ValuesMask{1,z} = [ValuesMask{1,z};false(1,1)]
        else ValuesMask{1,z} = [ValuesMask{1,z};true(1,1)]
        end
    end
end
figure(100+k)
freesurfer_statsurf_scalar(Values, ValuesMask, 'aparc', 'ScalarName', 'Canonical Loadings', 'ValueLimits', [-1 1])

brainFileName = sprintf([filepath_to_save fileprefix '/' fileprefix '_LV_Brain_%d.png'],k)
freesurfer_statsurf_savefig(brainFileName)
end

% correlation matrix
figure 
[cormat,cormat_pvals] = corr([T H]);
imagesc(cormat);
set(gca,'Ytick',[1:length(VarNames([xindices yindices]))])
set(gca,'Xtick',[1:length(VarNames([xindices yindices]))])
set(gca,'Xticklabel',VarNames([xindices yindices]),'FontWeight','bold','FontSize',9,'XTickLabelRotation', 45)
set(gca,'Yticklabel',VarNames([xindices yindices]),'FontWeight','bold','FontSize',9)
title('Correlations of Variables in CanCorr', 'FontSize', 24);
save([filepath_to_save fileprefix '.mat'],'cormat_pvals', 'mean_loads_all_X', 'mean_loads_all_Y', ...
    'X_times_loadings_LVs','Y_times_loadings_LVs',...
    'X_sig_loads_idx','Y_sig_loads_idx',...
    'X_LVs_sig_names', 'Y_LVs_sig_names','VarNames',...
    'mean_loads_all_X_vis','mean_loads_all_X_vis_orig','mean_loads_all_Y_vis',...
    'T', 'H', 'Xcoef', 'Ycoef', 'r1', 'Xscore', 'Yscore', 'stats', 'ABCD_X_VOL', 'cormat', 'ABCD_Y_VOL',...
    'Xcoef_orig', 'Ycoef_orig', 'r1_orig', 'Xscore_orig', 'Yscore_orig', 'stats_orig') 

%one at a time
%X = [ones(size(data(:,1))), H];
%[~, ~, r]=regress(T, X);


