clear all;
permutedata_dir = 'TuningSingle_permute';
realdata_dir = 'MeanSingle_ver2_TuningMean_TuningSingle';
FS_no = 'all';
Sub_lists = {'01','03','05','07','09','11','13','15','17','19','21','23'};
line_colors ={'b','c','r'};
ROI_lists = {'V1_thr','V2_thr','V3_thr','IPS_thr','SPL_thr','FEF_thr'};
ROI_labels = {'V1','V2','V3','IPS','SPL','FEF'};

ori_rad = cosd(abs((30:30:180)*2-180))

%% Mean Task  ====================================================
Att_lists = {'Mean_Mean', 'Mean_Single'};

% %% compute a vector mean of the 1d reconstruction for permute data
for ss=1:length(Sub_lists)   
    for rr=1:length(ROI_lists)     
        for dd=1:length(Att_lists)
            for niteration = 1:1001
                if niteration == 1001 %% compute a vector mean of the 1d reconstruction for real data 
                    FileName = sprintf('%s/%s_Mean_%s_%s_result_%s_shift.txt', ...
                    realdata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no);
                else
                    FileName = sprintf('%s/%s_Mean_%s_%s_result_%s_shift_permute_%i.txt', ...
                    permutedata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no, niteration);
                end
                fullCross_Mean.(Att_lists{dd}).(ROI_lists{rr})= load(FileName);
                dataIn   = mean(fullCross_Mean.(Att_lists{dd}).(ROI_lists{rr}),1); 
                
                % vector 
                ff = dataIn .* ori_rad;
                MeanVec.(Att_lists{dd})(ss,rr,niteration) = mean(ff);
                
            end
        end
    end
end


%% Single Task  ====================================================
Att_lists = {'Single_Mean', 'Single_Single'};
% %% compute a vector mean of the 1d reconstruction for permute data
for ss=1:length(Sub_lists)   
    for rr=1:length(ROI_lists)     
        for dd=1:length(Att_lists)
            for niteration = 1:1001
                if niteration == 1001 %% compute a vector mean of the 1d reconstruction for real data 
                    FileName = sprintf('%s/%s_Single_%s_%s_result_%s_shift.txt', ...
                    realdata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no);
                else
                    FileName = sprintf('%s/%s_Single_%s_%s_result_%s_shift_permute_%i.txt', ...
                    permutedata_dir, Sub_lists{ss}, Att_lists{dd}, ROI_lists{rr}, FS_no, niteration);
                end
                fullCross_Single.(Att_lists{dd}).(ROI_lists{rr})= load(FileName);
                dataIn   = mean(fullCross_Single.(Att_lists{dd}).(ROI_lists{rr}),1); 
                
                % vector 
                ff = dataIn .* ori_rad;
                SingleVec.(Att_lists{dd})(ss,rr,niteration) = mean(ff);
                
            end
        end
    end
end


