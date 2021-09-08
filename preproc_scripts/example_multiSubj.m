% basic multi-subject example for preprocessing

path = 'E:\AnnArbor_b_classic';

dirs = dir(fullfile(path, 's*'));
folders = cell(size(dirs));

for i = 1:size(dirs,1)
    folders{i} = fullfile(path,dirs(i).name);
end



for i = 1:size(folders,1)
    fprintf('%s%d%s%d\n','Processing subject number ', i, '. Total number of subjects is ', size(folders,1));
    % setup paths for this subject's anatomical data
    structPath = fullfile(folders{i},'anat');
    % ditto for functional data
    functPath = fullfile(folders{i}, 'func');
    % check sanity and start
    if(exist(structPath,'dir') && exist(functPath,'dir'))
        preprocess(functPath,structPath,'procChain',...
            {'realign','QC','coregister','newSegment','label'},...
            'atlasFile','mydata/atlases/Hn30r83.nii');
    end
end
