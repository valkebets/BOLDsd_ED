function overlapMaskNative(pathseg,pathMasks)
% find overlap between segmentation and mask in native space, threshold and
% save
% v1.0 Nov 2012 Nora Leonardi

fprintf('Creating new template masks for CSF and WM...');
for i=1:2
    if i==1
        % path to CSF segmentation file
        [segFN_noPath, ~] = spm_select('List',pathseg,'^c3.*\.nii$');
        % path to CSF mask
        [MaskFN_noPath, ~] = spm_select('List',pathMasks,'^CsfMask.*\.nii$');
    else
        % path to WM segmentation file
        [segFN_noPath, ~] = spm_select('List',pathseg,'^c2.*\.nii$');
        % path to WM mask
        [MaskFN_noPath, ~] = spm_select('List',pathMasks,'^WhiteMask.*\.nii$');
    end
    
    sSeg_fn=fullfile(pathseg,segFN_noPath);
    Vseg_i=spm_vol(sSeg_fn);
    [pth,nam,ext] = spm_fileparts(Vseg_i.fname);
    if exist([pth,'/w',nam,ext],'file'), fprintf('Mask %d already exists. ',i);
    else
        % read seg file
        Vseg=spm_read_vols(Vseg_i);
        % read mask file
        Mask_fn=fullfile(pathMasks,MaskFN_noPath); 
        Vmask_i=spm_vol(Mask_fn);
        Vmask=spm_read_vols(Vmask_i);

        % find overlap and threshold
        Vseg(~logical(Vmask))=0; % remove all CSF parts not in mask
        Vseg(Vseg<0.4)=0; % remove all parts with proba below 0.4 (c.f. FCON1000 preproc scripts)
        Vseg_i.fname=[pth,'/w',nam,ext];
        spm_write_vol(Vseg_i,Vseg);
    end
end
fprintf('Done\n');