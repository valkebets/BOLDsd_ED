%% Resting-state fMRI preprocessing for BOLD signal variability analysis %%
% This script computes the preprocessing steps listed below on
% resting-state fMRI data: 
%   Realignment
%   Segmentation
%   Co-registration
%   Normalisation to MNI space
%   Smoothing
%   Regress out WM and CSF signal
%   Bandpass filtering 
%   Motion scrubbing


clear

addpath(genpath('./spm12'));
addpath(genpath('./DPARSF'))
addpath(genpath('./preproc_scripts'));

            
%% Realignment, Coregistration, New Segmentation

struct_dir = %%% define path to anatomical file
funct_dir = %%% define path to functional timeseries

preprocess(funct_dir,struct_dir,'procChain',{'realign','QC','coregister','newSegment','label'},'atlasFile',atlas_file,'QCcoef',1.5,'highpass', 0);


%% Normalisation, Smoothing        

% Normalisation
y_file = dir(fullfile(struct_dir,'Segmented','y*.nii'));
matlabbatch{1}.spm.spatial.normalise.write.subj.def = {fullfile(struct_dir,'Segmented',y_file.name)};

r_files = dir(fullfile(funct_dir,'realigned','rf*.nii'));
for i = 1:length(r_files)
    norm_files{i} = fullfile(funct_dir,'realigned',r_files(i).name);
end

matlabbatch{1}.spm.spatial.normalise.write.subj.resample = norm_files;
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb = [-78 -112 -70
    78 76 85];
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox = [3 3 3];
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = 1;
matlabbatch{1}.spm.spatial.normalise.write.woptions.prefix = 'w';
spm_jobman('run',matlabbatch);

% Smoothing
wr_files = dir(fullfile(funct_dir,'realigned','wrf*.nii'));
for i = 1:length(wr_files)
    norm_files{i} = fullfile(funct_dir,'realigned',wr_files(i).name);
end

matlabbatch{1}.spm.spatial.smooth.data = norm_files;
matlabbatch{1}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';
spm_jobman('run',matlabbatch);        


%% Regress out WM & CSF signal

masks_dir = './preproc_scripts/Masks';

% Extract average signal in white matter & CSF masks
[CSFavg,WMavg] = getMaskAverageCSF_WM(funct_dir,struct_dir,masks_dir,'swrf');
x_csf(1:length(CSFavg),1) = CSFavg;
x_wm(1:length(WMavg),1) = WMavg;
predictors = [x_csf x_wm];

% Vectorize timeseries
swr_files = dir(fullfile(funct_dir,'swrf*.nii'));
for i = 1:length(swr_files)
    Ai = spm_vol(fullfile(funct_dir,swr_files(i).name));
    A = spm_read_vols(Ai);
    A_vec = reshape(A,[Ai.dim(1)*Ai.dim(2)*Ai.dim(3) 1]);
    y(i,:) = A_vec;
end

nTimePoints = size(y,1);
nVoxels = size(y,2);
y_predicted = nan(nTimePoints,nVoxels);

% Compute linear regression using WM+CSF signal as predictors 
for v = 1:nVoxels
    ts = squeeze(y(:,v)); 
    ts(isnan(ts)) = 0;
    confounds = squeeze(predictors(:,:));
    B = LinearModel.fit(confounds,ts);
    [this_y_predicted,~] = predict(B,confounds);
    y_predicted(:,v) = this_y_predicted;
end

% Discard predicted Y
y_final = y - y_predicted;
y_final(isnan(y_final))=0;


%% Bandpass filtering
    
myTR = 2.1; % fMRI repetition time
myBand = [0.01 0.1]; % frequency band for bandpass filtering

ts_filtered = y_IdealFilter(y_final,myTR,myBand); % uses DPARSF
        

%% Motion scrubbing

artifacts_file = fullfile(funct_dir,'artifacts.mat');
[ts_scrubbed,indd] = scrubTimeSeries(ts_filtered,artifacts_file,funct_dir,0,'cut',0.5);

 