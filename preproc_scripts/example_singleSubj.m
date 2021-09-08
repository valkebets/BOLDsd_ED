% basic single-subject example for preprocessing
clc;

p.atlasPath='/Users/Richiardi/_skool/works/biosignals/MRI/data/atlas_AAL_correctLR/AAL116_correctLR.nii';

p.rootPath='/Users/Richiardi/temp/testProc/c03';
p.structPath = fullfile(p.rootPath,'ana');
p.functPath = fullfile(p.rootPath,'resting');

preprocess(p.functPath,p.structPath,'procChain',{'realign','QC','coregister','newSegment','label'},...
		   'QCcoef',1.5,'highpass', 0,'atlasFile',p.atlasPath);
