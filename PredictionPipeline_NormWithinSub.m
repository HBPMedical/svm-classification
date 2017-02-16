%==========================================================================
%This script is used to compare different parameter influence on the
%prediction. The pathologically proven (PP) subjects are used here. 
%==========================================================================

cd C:\DATA\ExternalDrive\Jing_work\LREN\Projects\
rmpath(genpath('C:\Program Files\MATLAB\R2014b\toolbox\stats\stats\'))

%==========================================================================
%EXPERIMENT PARAMETER 

    %correct age&tiv effects based on healthy controls ("hc") or all subjects ("all")
    correct_type='hc'; 
    %LOO indicates if the validation subject is used for the estimation of beta
    LOO=1;
    linear='agetiv'; % remove effect of age and tiv ('agetiv'); remove effect of age ('age')
     
    ResultsDir='ADsubgrouping\C3_Imaging\TivAgeLinear_NormWithinSub\';
    ResultsFile='Dijon_nonAD_MRIage_HcTivAgeLinear_NormWithinSub';

%==========================================================================
%load data

Mask= spm_vol('SVM/RESULTS/ADvsHC_ADNI_PP/GM_maskADNIdifPP/GM_025maskADNIdifPP1.nii');
mask = spm_read_vols(Mask);

% info of PP cohort
PPsubPTID = readtable('SVM/PPsubjects_data/PPsubjects_TrueLabel.xlsx');
HC = PPsubPTID.HC(cellfun(@isempty,PPsubPTID.HC)~=1,1);
AD = PPsubPTID.AD(cellfun(@isempty,PPsubPTID.AD)~=1,1);
PPSub = [HC;AD]; % HC + AD =Sub
PPlabel=[repmat(1, 1, length(HC)),repmat(2, 1,length(AD))]'; %HC is 1;AD is 2 

TrainImages = strcat('SVM/PP_preproc/image_preproc/smwc1',PPSub,'.nii');
PPinfo = readtable('SVM/PP_preproc/PPTIV_John.xls');

% info of ADNI cohort
ADNIad = readtable('ADsubgrouping\ADNI_Imaging_Data\adPtidTivAge.xls');
ADNIhc = readtable('ADsubgrouping\ADNI_Imaging_Data\cnPtidTivAge.xls');

ADNISub = [ADNIhc.PTID;ADNIad.PTID]; % HC + AD =Sub
PredictTiv = [ADNIhc.TIV;ADNIad.TIV];
PredictAge = [ADNIhc.AGE;ADNIad.AGE];
Dir = 'ADNI_preprocess\good_seg_PPnormalize_SmoothResults_SMWC1\';
PredictImages = strcat(Dir,'smwc1',ADNISub,'.nii');
PredictLabel = [ones(1, size(ADNIhc,1)),repmat(2, 1,size(ADNIad,1))]';  %HC is 1;AD is 2 

PredictID=ADNISub;

% % Dijon
% C3 = readtable('ADsubgrouping\C3_Imaging\Dijon_info_MRIage.csv');
% Dir='ADsubgrouping\C3_Imaging\Dijon_PPnorm_smwc1\Dijon_PPnorm_smwc1\';
% PredictImages = strcat(Dir,C3.path,'.nii');
% PredictTiv = C3.Tiv;
% PredictAge = C3.Age;
% PredictLabel = C3.DEMT0;
% 
% PredictID=C3.path;

% Dijon non AD 
C3 = readtable('ADsubgrouping\C3_Imaging\Dijon_nonAD.csv');
Dir='ADsubgrouping\C3_Imaging\Dijon_PPnorm_smwc1\Dijon_PPnorm_smwc1\';
PredictImages = strcat(Dir,C3.path,'.nii');
PredictTiv = C3.Tiv;
PredictAge = C3.Age;
PredictLabel = C3.DEMT0;

PredictID=C3.path;


% Bordeaux
% Bordeaux_AD = readtable('ADsubgrouping\C3_Imaging\Bordeaux_AD.csv');
% Bordeaux_HC = readtable('ADsubgrouping\C3_Imaging\Bordeaux_HC.csv');
% 
% Bordeaux_HC_MRI=strcat('ADsubgrouping\C3_Imaging\Bordeaux_HC\',Bordeaux_HC.Image);
% Bordeaux_AD_MRI=cellstr(strcat('ADsubgrouping\C3_Imaging\Bordeaux_AD\',Bordeaux_AD.Images));
% 
% PredictImages = [Bordeaux_HC_MRI;Bordeaux_AD_MRI];
% PredictTiv = [Bordeaux_HC.TIV;Bordeaux_AD.TIV];
% PredictAge = [Bordeaux_HC.Age;Bordeaux_AD.MRIage];
% PredictLabel = [repmat(1, 1, size(Bordeaux_HC,1)),repmat(2, 1,size(Bordeaux_AD,1))]';
% 
% PredictID=[num2cell(Bordeaux_HC.ID);num2cell(Bordeaux_AD.num)];
%==========================================================================
%order the subjects as the same sequence of Sub
[ismem,loc] = ismember(PPSub,PPinfo.Sub);
pps=PPinfo.Sub(loc');
TrainTiv=PPinfo.TIV(loc');
TrainAge=PPinfo.Age(loc');

for i = 1:length(TrainImages)
    volume = spm_vol(TrainImages{i});
    cc = spm_read_vols(volume).*(abs(det(volume.mat))/100^3);
    cc_v(:,i)=cc(mask(:)>0);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare covariate matrix and remove confunds effects

if isequal({'all'},{correct_type})
    l=length(PPSub);
elseif isequal({'hc'},{correct_type})
    l=length(HC);
end

%use all normal subjects to estimate the effect of tive and age
if isequal({'age'},{linear})
    cov_train=[ones(1,l);TrainAge(1:l,1)']';
elseif isequal({'agetiv'},{linear})
    cov_train=[ones(1,l);TrainTiv(1:l,1)';TrainAge(1:l,1)']';
end

B1=inv(cov_train'*cov_train)*cov_train'*cc_v(:,1:l)';

if isequal({'age'},{linear})
    B1=B1(2,:);
    volres_allnor=cc_v-B1'*TrainAge(:,1)';
%     volres_allnor=cc_v-B1'*[ones(length(PPSub),1),TrainAge(:,1)]';
elseif isequal({'agetiv'},{linear})
    B1=B1(2:3,:);
    volres_allnor=cc_v-B1'*[TrainTiv(:,1)';TrainAge(:,1)'];
%     volres_allnor=cc_v-B1'*[ones(length(PPSub),1)';TrainTiv(:,1)';TrainAge(:,1)'];
end

MEAN_sub=mean(volres_allnor);
STD_sub=std(volres_allnor);
volres_allnor=volres_allnor-repmat(MEAN_sub,size(volres_allnor,1),1);
volres_allnor=volres_allnor.*repmat((1./STD_sub),size(volres_allnor,1),1);
    
kernel_allnor=volres_allnor'*volres_allnor;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate the cost value 

C_pow=[-6:1:-2];
C = 10.^C_pow;

for j = 1:length(C)
    %use the 33 pp subjects for LOO cross validation
    for valsub =1:length(PPSub)
        % use each subject sequently as validation subject 
        train_id=1:length(PPSub)~=valsub;
        
        if LOO==1
            if sum(valsub==1:l)>0 
                % if current validation subject is a normal subject, beta is
                % estimated without using this normal subject. 
                cov_oneNL=cov_train(1:l~=valsub,:);
                beta=inv(cov_oneNL'*cov_oneNL)*cov_oneNL'*cc_v(:,1:l~=valsub)';
                
                if isequal({'age'},{linear})
                    beta=beta(2,:);
                    vol_res = cc_v-beta'*TrainAge(:,1)'; % remove the effects of age and tiv
%                     vol_res = cc_v-beta'*[ones(length(PPSub),1),TrainAge(:,1)]'; % only take the residuals
                elseif isequal({'agetiv'},{linear})
                    beta=beta(2:3,:);
                    vol_res = cc_v-beta'*[TrainTiv(:,1)';TrainAge(:,1)'];
%                     vol_res = cc_v-beta'*[ones(length(PPSub),1)';TrainTiv(:,1)';TrainAge(:,1)'];
                end
                
                MEAN_sub=mean(vol_res);
                STD_sub=std(vol_res);
                vol_res=vol_res-repmat(MEAN_sub,size(vol_res,1),1);
                vol_res=vol_res.*repmat((1./STD_sub),size(vol_res,1),1);
                
                kernel=vol_res'*vol_res;
                % current subject is used for validation
                test_kernel=kernel(valsub,train_id); 
                % the rest subjects are for training
                train_kernel=kernel(train_id,train_id); 
            else
                test_kernel=kernel_allnor(valsub,train_id); 
                train_kernel=kernel_allnor(train_id,train_id); 
            end
        else
            test_kernel=kernel_allnor(valsub,train_id); 
            train_kernel=kernel_allnor(train_id,train_id); 
        end
        m=svmtrain(PPlabel(train_id,1),[(1:length(train_kernel))' train_kernel], ['-t 4 -c ',num2str(C(j))]);
        [pred_lab,acc]=svmpredict(PPlabel(valsub,1),[(1:size(test_kernel,1))' test_kernel],m);
        accuracy(valsub,j)=pred_lab;
        clear m
    end
    %build a classifier with current beta and all subjects to record the
    %support vectors that current classifier is using
    svm = svmtrain(PPlabel,[(1:length(kernel_allnor))' kernel_allnor], ['-t 4 -c ',num2str(C(j))]);
    SV(j,1)=sum(svm.nSV);
    clear svm
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%performance of the classifier


for m = 1:length(C)
    plotY(m,1) = sum(accuracy(:,m)==PPlabel)./33;
end

fprintf('Highest accuracy is %i\n', max(plotY(:,1)));
C_min=min(C(ismember(plotY(:,1), max(plotY(:,1)))));
fprintf('C min value that provide highest accuracy is %i\n', C_min);
fprintf('Correctly classified healthy controls are %i\n', sum(accuracy(1:15,C==C_min)==1));
fprintf('Correctly classified patients are %i\n', sum(accuracy(16:33,C==C_min)==2));

svm = svmtrain(PPlabel,[(1:length(kernel_allnor))' kernel_allnor], ['-t 4 -c ',num2str(C_min)]);


plotX= C_pow;
plot(plotX,plotY)
title('accuracy with LOO method (with differ mask, volumes are corrected by age and tiv)');
xlabel('cost (log 10 scale)');
ylabel('accuracy');

figure(3)
plot(plotX,SV)
ylim_min=min(SV)-(max(SV)-min(SV))/10;
ylim_max=max(SV)+(max(SV)-min(SV))/10;
ylim([ylim_min,ylim_max]);
title('number of support vectors that are used for building the classifier');
xlabel('cost (log 10 scale)');
ylabel('number of SV');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%regions with high weights
w=(svm.sv_coef.*PPlabel(svm.SVs))'* volres_allnor(:,svm.SVs)';

spm_check_registration(Mask.fname);
grey_mask=spm_vol(Mask.fname);
[cc xyz]=spm_read_vols(grey_mask);
cc(mask>0)=w;
xyz(4,:)=1;
xyz=inv(grey_mask.mat)*xyz(:,abs(cc)>=0.0015);
int_xyz=randn(1,size(xyz,2));
spm_orthviews('AddBlobs',1,xyz,int_xyz,grey_mask.mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predict in ADNI cohort
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear cc cc_v y vol_res 

if isequal({'age'},{linear})
    cov_predict=PredictAge(:,1)'; 
%     cov_predict=[ones(length(PredictID),1)';PredictAge(:,1)']; %covariate matrix [tiv,age]
elseif isequal({'agetiv'},{linear})
    cov_predict=[PredictTiv(:,1),PredictAge(:,1)]'; %covariate matrix [tiv,age]
%     cov_predict=[ones(length(PredictID),1),PredictTiv(:,1),PredictAge(:,1)]';
end

%read in GM volumes
for i = 1:length(PredictImages)
    volume = spm_vol(PredictImages{i});
    cc = spm_read_vols(volume)*(abs(det(volume.mat))/100^3);
    cc = cc(mask(:)>0);
    
    cc=cc-B1'*cov_predict(:,i);
    res=(cc-mean(cc))/std(cc);
    
%     The score is calculated using the weights from SVM and the normalized
%     volume of a new subject
%     score(1,i)= w*res;

    test_kernel = res'*volres_allnor;
    [SVM_label,acc]=svmpredict(PredictLabel(i,1),[(1:size(test_kernel,1))' test_kernel],svm);
    pred_label(i,1)=SVM_label;
end


labels = table(PredictID,PredictLabel,pred_label);

if exist(ResultsDir)==0
    mkdir(ResultsDir)
end
writetable(labels,strcat(ResultsDir,ResultsFile,'.txt'));




