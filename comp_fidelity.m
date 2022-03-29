clear all;
permutedata_dir = 'MeanSingle_ver2_TuningMean_TuningSingle_permute';
realdata_dir = 'MeanSingle_ver2_TuningMean_TuningSingle';
FS_no = 'all';
Sub_lists = {'01','03','05','07','09','11','13','15','17','19','21','23'};
SD_lists  = {'SD7'};
line_colors ={'b','c','r'};
nConds = 6;
ROI_lists = {'V1_thr','V2_thr','V3_thr','IPS_thr','SPL_thr','FEF_thr'};
ROI_labels = {'V1','V2','V3','IPS','SPL','FEF'};

figNo = 1; 
nChan=6; c_center=3;
%ori_unit = pi/nChan;
%ori_rad = 0:ori_unit:(pi-ori_unit);
%ori_rad = ori_rad - ori_rad(c_center);
%tt = cos(abs(ori_rad))
ori_rad = cosd(abs((30:30:180)*2-180))
% plot(ori_rad)

%% Mean Task  ====================================================.
Att_lists = {'Mean_Mean', 'Mean_Single'};

% %% compute a vector mean of the 1d reconstruction for permute data
for ss=1:length(Sub_lists)   
    for rr=1:length(ROI_lists)     
        for dd=1:length(Att_lists)
            for niteration = 1:1001
                if niteration == 1001 %% compute a vector mean of the 1d reconstruction for real data 
                    FileName = sprintf('%s/%sAvg_fullCross_Mean_%s_%s_result_%s_shift.txt', ...
                    realdata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no);
                else
                    FileName = sprintf('%s/%sAvg_fullCross_Mean_%s_%s_result_%s_shift_permute_%i.txt', ...
                    permutedata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no, niteration);
                end
                fullCross_Mean.(Att_lists{dd}).(ROI_lists{rr})= load(FileName);
                dataIn   = mean(fullCross_Mean.(Att_lists{dd}).(ROI_lists{rr}),1); 
                
                % vector 
%                 ff = dataIn .* cos(abs(ori_rad));
                ff = dataIn .* ori_rad;
                MeanVec.(Att_lists{dd})(ss,rr,niteration) = mean(ff);
                
            end
        end
    end
end

for niteration = 1:1001
    [h,p,ci,stats]= ttest(MeanVec.Mean_Mean(:,niteration));
    MeanVec.Mean_Mean_tstat(niteration) = stats.tstat;
    
    [h,p,ci,stats]= ttest(MeanVec.Mean_Single(:,niteration));
    MeanVec.Mean_Single_tstat(niteration) = stats.tstat;
end

% sort  
Tempdata= [];     
Tempdata = sort(MeanVec.Mean_Mean_tstat, 'descend');
TempRank = find(Tempdata == MeanVec.Mean_Mean_tstat(1001));
MeanVec.Mean_Mean_pval = (TempRank/10)/100;
fprintf('\n %2.4f\n', MeanVec.Mean_Mean_pval)

Tempdata= [];     
Tempdata = sort(MeanVec.Mean_Single_tstat, 'descend');
TempRank = find(Tempdata == MeanVec.Mean_Single_tstat(1001));
MeanVec.Mean_Single_pval = (TempRank/10)/100;
fprintf('\n %2.4f\n',MeanVec.Mean_Single_pval)


%% =================
% % vec vs. 0
% % Mean Mean Mean, Mean Mean Single 
% % compute tstat
for dd=1:length(Att_lists)
    for rr=1:length(ROI_lists)    
        for niteration = 1:1001
            SEM = std(MeanVec.(Att_lists{dd})(:,rr,niteration))/sqrt(length(Sub_lists));
            MeanVec_tstat.(Att_lists{dd})(niteration,rr) =  ...
            mean(MeanVec.(Att_lists{dd})(:,rr,niteration))/SEM;
        end
    end
end

%% sort ascend : Fidelity
for dd=1:length(Att_lists)   
    figure(dd),clf;
    for rr=1:length(ROI_lists)
        Tempdata = sort(MeanVec_tstat.(Att_lists{dd})(:,rr), 'descend');
        TempRank = find(Tempdata == MeanVec_tstat.(Att_lists{dd})(1001,rr));
        MeanVec_pval.(Att_lists{dd})(rr) = (TempRank/10)/100;

        subplot(1,6,rr); 
        hist(MeanVec_tstat.(Att_lists{dd})(:,rr), 50); hold on
        % plot([0, 0], [0 100] , 'k--');
        plot([MeanVec_tstat.(Att_lists{dd})(1001,rr), MeanVec_tstat.(Att_lists{dd})(1001,rr)], ...
            [0 100], 'r','linewidth',2);
    %         xlim([0 150]);  ylim([0 100]);
        title(sprintf('%s',ROI_lists{rr}));
    end
    
    fprintf('MeanVec_%s', Att_lists{dd})
    MeanVec_pval.(Att_lists{dd}) 
end

%% Single 
% Task, Tuning, Orientation (Mean Mean Mean, Mean Mean Single ...)
%% Single Task  ====================================================.
Att_lists = {'Single_Mean', 'Single_Single'};

% %% compute a vector mean of the 1d reconstruction for permute data
for ss=1:length(Sub_lists)   
    for rr=1:length(ROI_lists)     
        for dd=1:length(Att_lists)
            for niteration = 1:1001
                if niteration == 1001 %% compute a vector mean of the 1d reconstruction for real data 
                    FileName = sprintf('%s/%sAvg_fullCross_Single_%s_%s_result_%s_shift.txt', ...
                    realdata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no);
                else
                    FileName = sprintf('%s/%sAvg_fullCross_Single_%s_%s_result_%s_shift_permute_%i.txt', ...
                    permutedata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no, niteration);
                end
                fullCross_Single.(Att_lists{dd}).(ROI_lists{rr})= load(FileName);
                dataIn   = mean(fullCross_Single.(Att_lists{dd}).(ROI_lists{rr}),1); 
                
                % vector 
%                 ff = dataIn .* cos(abs(ori_rad));
                ff = dataIn .* ori_rad;
                SingleVec.(Att_lists{dd})(ss,rr,niteration) = mean(ff);
                
            end
        end
    end
end

SingleSingleMean = SingleVec.Single_Mean(:,:,1001);
SingleSingleSingle = SingleVec.Single_Single(:,:,1001);

for niteration = 1:1001
    [h,p,ci,stats]= ttest(SingleVec.Single_Mean(:,niteration));
    SingleVec.Single_Mean_tstat(niteration) = stats.tstat;
    
    [h,p,ci,stats]= ttest(SingleVec.Single_Single(:,niteration));
    SingleVec.Single_Single_tstat(niteration) = stats.tstat;
end

% sort  
Tempdata= [];     
Tempdata = sort(Single_linear.Single_Mean_tstat, 'descend');
TempRank = find(Tempdata == Single_linear.Single_Mean_tstat(1001));
Single_linear.Single_Mean_pval = (TempRank/10)/100;
fprintf('\n %2.4f\n', Single_linear.Single_Mean_pval)

Tempdata= [];     
Tempdata = sort(Single_linear.Single_Single_tstat, 'descend');
TempRank = find(Tempdata == Single_linear.Single_Single_tstat(1001));
Single_linear.Single_Single_pval = (TempRank/10)/100;
fprintf('\n %2.4f\n',Single_linear.Single_Single_pval)

%% =================
% % vec vs. 0
% % Single Single Mean, Single Single Single 
% % compute tstat
for dd=1:length(Att_lists)
    for rr=1:length(ROI_lists)    
        for niteration = 1:1001
            SEM = std(SingleVec.(Att_lists{dd})(:,rr,niteration))/sqrt(length(Sub_lists));
            SingleVec_tstat.(Att_lists{dd})(niteration,rr) =  ...
            mean(SingleVec.(Att_lists{dd})(:,rr,niteration))/SEM;
        end
    end
end

%% sort ascend : Fidelity
for dd=1:length(Att_lists)   
    figure(dd),clf;
    for rr=1:length(ROI_lists)
        Tempdata = sort(SingleVec_tstat.(Att_lists{dd})(:,rr), 'descend');
        TempRank = find(Tempdata == SingleVec_tstat.(Att_lists{dd})(1001,rr));
        SingleVec_pval.(Att_lists{dd})(rr) = (TempRank/10)/100;

        subplot(1,6,rr); 
        hist(SingleVec_tstat.(Att_lists{dd})(:,rr), 50); hold on
        % plot([0, 0], [0 100] , 'k--');
        plot([SingleVec_tstat.(Att_lists{dd})(1001,rr), SingleVec_tstat.(Att_lists{dd})(1001,rr)], ...
            [0 100], 'r','linewidth',2);
    %         xlim([0 150]);  ylim([0 100]);
        title(sprintf('%s',ROI_lists{rr}));
    end
    
    fprintf('MeanVec_%s', Att_lists{dd})
    SingleVec_pval.(Att_lists{dd}) 
end


%% plot 
%% TunignAll, Mean, Single orientation plotting
% %% Mean orientation during Mean Task(blue), Single Task(cyan)
% % Data to be plotted as a bar graph
AllMean = [mean(MeanMeanMean,1); mean(SingleSingleSingle,1)]';
Allerror = [std(MeanMeanMean)/sqrt(length(Sub_lists)); ...
    std(SingleSingleSingle)/sqrt(length(Sub_lists))]';

% Creating axes and the bar graph
figure(22),clf;
ax = axes;
h = bar(AllMean,'BarWidth',1.0);
h(1).FaceColor = 'blue';
h(2).FaceColor = 'magenta'; hold on;

% Finding the number of groups and the number of bars in each group
ngroups = size(AllMean, 1);
nbars = size(AllMean, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    % Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, AllMean(:,i), Allerror(:,i));
end
ylim([-0.02 0.1]);
set(gca,'XTick',[1:length(ROI_labels)]);
set(gca, 'XTickLabel', ROI_labels);





