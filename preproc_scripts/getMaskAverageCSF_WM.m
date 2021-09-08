function [CSFavg,WMavg] = getMaskAverageCSF_WM(pathRS,pathSeg,masksPath,prefix)

% Based on Jonas Richiardi & Nora Leonardi's timeCourseExtractor script

cd(fullfile(masksPath));
CSFthreshold = 0.95; % threshold in CSF tissue probability map to consider a voxel as CSF
WMthreshold = 0.95; % same for WM

%% Warp CSF and WM masks to subject space OR reslice
Masks = {'WhiteMask_09_121x145x121.nii','CsfMask_07_121x145x121.nii'};
deform_field = dir(fullfile(pathSeg,'iy*.nii')); % from "new segment", template to native space deformation field
for j = 1:length(Masks)
    tmp = which(Masks{j});
    if ~isempty(tmp) && ~exist(fullfile(pathSeg,['w',Masks{j}]),'file'),
        fprintf('Warping mask %s to native space...\n',Masks{j});
        nii_mask2target(tmp,fullfile(pathSeg,deform_field.name));
    elseif isempty(tmp)
        error('No mask of right size found. Please reslice the masks of DPARSFA to the dimensions 121x145x121');
    elseif exist(fullfile(pathSeg,['w',Masks{j}]),'file')
        fprintf('Warped mask %s already exists\n',['w',Masks{j}]);
    end
end
% find overlap between mask and segmentation
%overlapMask(paths.sseg) % saves c2 and c3 with prefix 'w'
overlapMaskNative(pathSeg,masksPath) % saves c2 and c3 with prefix 'w'


%% Extract CSF & WM

fprintf('Extracting CSF...');
[segFN_noPath, foosubdirs] = spm_select('List',pathSeg,'^wc3.*\.nii$'); % new file with overlap with mask
sSeg_fn{1} = fullfile(pathSeg,segFN_noPath); clear segFN_noPath foosubdirs
fprintf('Extracting WM...');
[segFN_noPath, foosubdirs] = spm_select('List',pathSeg,'^wc2.*\.nii$');
sSeg_fn{2} = fullfile(pathSeg,segFN_noPath);
[fMeanFN_noPath, foosubdirs] = spm_select('List',pathRS,['^' prefix '.*\.nii$']);
if isempty(fMeanFN_noPath)
    [fMeanFN_noPath, foosubdirs] = spm_select('List',pathRS,['^' prefix '.*\.img$']);
end
fMeanFN_noPath = fMeanFN_noPath(1,:); % take first file
fMean_fn = fullfile(pathRS,fMeanFN_noPath);
CSFV = mapVolumeToVolume(sSeg_fn{1},fMean_fn);
CSFmaskLidx = CSFV > CSFthreshold;
WMV = mapVolumeToVolume(sSeg_fn{2},fMean_fn); % map seg to mean voxel space
WMmaskLidx = WMV > WMthreshold;


%% Get list of filenames of the images in temporal order & read volumes
scansIdx={[]};
fVolsFNlist = getImageFNinAcqOrder(pathRS,prefix,'nii',scansIdx{1}); %%doesn't work with some datasets file format
fVolsFNlist_fullpath = cellfun(@(x) fullfile(pathRS,x), fVolsFNlist,'UniformOutput',false);

% If doesn't work with data, then uncomment the next 4 ligns et comments the previous 2
%ffiles = dir(fullfile(pathRS,[prefix '*.nii']));
%for iter = 1:length(ffiles)
%    fVolsFNlist_fullpath{iter,1} = fullfile(pathRS,ffiles(iter,1).name);
%end

% read all headers and files in temporal order into a 4-D array
V0i = spm_vol(fVolsFNlist_fullpath);
V0idx = 1:length(V0i);
% clear and preallocate to avoid fragmenting
clear V0;
V0 = zeros(V0i{1}.dim(1),V0i{1}.dim(2),V0i{1}.dim(3),length(V0i),'single');

% read all volumes
fprintf('Reading all volumes\n');
for iter = 1:length(V0i),
    V0(:,:,:,iter) = spm_read_vols(V0i{iter});
end


%% Compute average CSF & WM signal

disp([' Computing average CSF signal...']);
% average signal of all CSF voxels corresponding to this block
[Ci,Cj,Ck] = ind2sub(size(V0(:,:,:,1)),find(CSFmaskLidx));
%CSFsignal=V0(Ci,Cj,Ck,V0idx);
nCSFvoxels = numel(Ci); % also, sum(CSFmaskLidx)
tcCSF = zeros(1,1,1,size(V0,4));
%[csfx1,csfx2,csfx3]=ndgrid(1:szCSFvoxels(1),1:szCSFvoxels(2),1:szCSFvoxels(3));
for idx = 1:numel(Ci)
    tcCSF = tcCSF+V0(Ci(idx),Cj(idx),Ck(idx),:);
end
CSFavg = squeeze(tcCSF/nCSFvoxels)';
clear tcCSF Ci Cj Ck

disp([' Computing average WM signal...']);
% average signal of all WM voxels corresponding to this block
[Ci,Cj,Ck] = ind2sub(size(V0(:,:,:,1)),find(WMmaskLidx));
%CSFsignal=V0(Ci,Cj,Ck,V0idx);
nWMvoxels = numel(Ci); % also, sum(CSFmaskLidx)
tcWM = zeros(1,1,1,size(V0,4));
%[csfx1,csfx2,csfx3]=ndgrid(1:szCSFvoxels(1),1:szCSFvoxels(2),1:szCSFvoxels(3));
for idx = 1:numel(Ci)
    tcWM = tcWM + V0(Ci(idx),Cj(idx),Ck(idx),:);
end
WMavg = squeeze(tcWM/nWMvoxels)';

end
