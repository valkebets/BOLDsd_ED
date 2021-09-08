%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3944 $)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.preproc.data = '<UNDEFINED>';
matlabbatch{1}.spm.spatial.preproc.output.GM = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.output.WM = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.output.CSF = [1 1 1];
matlabbatch{1}.spm.spatial.preproc.output.biascor = 1;
matlabbatch{1}.spm.spatial.preproc.output.cleanup = 0;
matlabbatch{1}.spm.spatial.preproc.opts.tpm = {
                                               fullfile(spm('Dir'),'tpm','grey.nii')
                                               fullfile(spm('Dir'),'tpm','white.nii')
                                               fullfile(spm('Dir'),'tpm','csf.nii')
                                               };
matlabbatch{1}.spm.spatial.preproc.opts.ngaus = [2
                                                 2
                                                 2
                                                 4];
matlabbatch{1}.spm.spatial.preproc.opts.regtype = 'mni';
matlabbatch{1}.spm.spatial.preproc.opts.warpreg = 1;
matlabbatch{1}.spm.spatial.preproc.opts.warpco = 25;
matlabbatch{1}.spm.spatial.preproc.opts.biasreg = 0.0001;
matlabbatch{1}.spm.spatial.preproc.opts.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.opts.samp = 3;
matlabbatch{1}.spm.spatial.preproc.opts.msk = {''};
