% This script computes voxelwise BOLD signal variability on preprocessed 
% fMRI timeseries
%
% Valeria Kebets, September 2021

clear

addpath(genpath('./spm12'));

template_file = %%% specify path & name of DARTEL group template or other template in MNI space

% Load template
Ai = spm_vol(template_file);
A = spm_read_vols(Ai);
mask = find(A~=0);

% For each voxel, compute standard deviation of BOLD signal across timecourse
for i = 1:size(ts_scrubbed,1)
    ts_scrubbed_sd(i,1) = std(ts_scrubbed(i,:));
end

% Constrain to voxels within template
ts_scrubbed_sd_masked = ts_scrubbed_sd(mask);

% Spatial z-scoring (within-subject)
ts_scrubbed_sd_masked_z = zscore(ts_scrubbed_sd_masked);

% Write BOLD variability map to volume
Bi = Ai;
Bi.fname = fullfile('BOLDsd_' <subj_name> '.nii');
Bi.dt = [spm_type('float32') 0];
B = zeros(size(A));
B(Aidx) = ts_scrubbed_sd_masked_z;
spm_write_vol(Bi,B);

