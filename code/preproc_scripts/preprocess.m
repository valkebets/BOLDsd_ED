function preprocess(functPath, structPath, varargin)
% Preprocessing pipeline for fMRI connectivity analysis (functional
% realignment, functional-structural coregistration, segmentation,
% normalisation, atlasing)
%
% IN:
%   functPath: string, full path to raw functional data
%   structPath: string, full path to raw structural data
%   varargin:
%       1)  the processing chain to use, e.g.
%       'procChain',{'realign','QC','coregister','newSegment','label'}
%       2) the name of the atlas, e.g.
%       'atlasFile','myAtlas.nii'
%       3) other arguments for quality control, e.g.
%       'QCcoef', 1.5
%       4) high-pass filtering, e.g.
%       'highpass', 0
% 
% USAGE
%   Please see the "example_*.m" scripts for an example call.
%
% REQUIREMENTS
% - SPM12
% 
% INSTALLATION
% - SPM12 must be on your matlab path
%
% ACKNOWLEDGEMENTS
% This file borrows liberally from the following sources
% - Rik Henson's example script for SPM 5 (MRC Cognition and Brain Sciences Unit)
% http://en.wikibooks.org/wiki/SPM-Example_batch_script
% - Various functions in the IBASPM toolbox, in particular auto_labelling
% from Yasser Aleman-Gomez and colleages at the Cuba Neurosciences Center
% - highpassFilter.m is art_despike from Paul Mazaika
% (http://www.stanford.edu/~mazaika/)
% - editfilenames.m is from SPM (http://www.fil.ion.ucl.ac.uk/spm/)
%
% If you use this script in published work please cite
% for cognitive data: Richiardi et al., Decoding brain states from fMRI
%   connectivity graphs, NeuroImage (56), 2011, pp. 616--626
% for clinical data: Richiardi et al., Classifying minimally-disabled
%   multiple sclerosis patients from resting-state functional connectivity,
%   NeuroImage (62), 2012, pp. 2021--2033
%
%
% VERSION HISTORY
% v1.0 Jonas Richiardi / Medical Image Processing Laboratory
% - initial release for SPM5
% v2.0 Manuel Wuethrich
% - SPM8 support
% - use spm jobmanager
% - support for newSegment and unifiedSegmentation
% - port IBASPM auto_labelling to SPM8
% - basic QC capability
% v2.0.1 Jonas Richiardi
% - cross-platformish
% - path settings fixed
% v2.0.1b Jonas Richiardi
% - removed ~ ignore code and changed to 4-output forme for filepart for 
% backwards compatibility with pre-7.9 versions e.g. as in Greedy SCC cluster
% - corrected path setting code
% v2.0.2 June 2011 Jonas Richiardi
% - supports smoothing
% - supports coregistration direction setting (FtoS or StoF)
% v2.0.3 June 2011 JR
% - supports custom atlases
% v2.1 Oct 2011 JR
% - fixes a critical bug introduced in v2.0 - wrong interpolation
% parameters were used in new_Auto_Labelling.m
% v2.1.1 Dec 2011 JR
% - stricter dir-and-file checking in pp_loadVolumes
% v2.1.2 Dec 2011 JR
% - add version tagging to tasksDone structure to help traceability
% - bugfix and better existence type checking in pp_loadVolume
% - support for Matlab 7.13 (fileparts syntax)
% v2.1.3 Nov 2012 JR
% - fixed typo in classic segment code path
% v2.1.4 Oct 2013 JR / Stanford University
% - bit of cleanup code and doc, atlas file mandatory
% v2.1.5 July 2015 JR 
% - fixed bug when no atlasing present (thanks Martin Ndengera)
% v3.0 July 2015 JR
% - SPM12 support (realign syntax, segment_job.m... updated)
%

DEBUGMODE=false; % set this to true to help diagnosing problems

%% check and set necessary paths
if exist('spm.m','file')~=2
    error('SPM seems not to be on your matlab path.');
else
    spmLoc=spm('Dir');
    myPath=path();
    % add 'config' SPM dir to the path
    spmConfigLoc=[spmLoc filesep 'config'];
    if exist(spmConfigLoc,'dir')~=7
        error(['SPM installation seems to be missing a config dir at '...
            spmConfigLoc ', quitting.']);
    else
        addpath(spmConfigLoc);
    end
end

if (DEBUGMODE==true)
    path
end

%% initialize and load data
warning('Slice-timing correction should be added to the chain for long TRs');
disp('** STARTING PROCESSING');
Tstart = clock;
spm('defaults', 'fMRI');
spm_jobman('initcfg');
pp_loadVolumes;
% generate full path of this file
thisLocation=which('preprocess.m');
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, quux1, quux2, quux3] = fileparts(thisLocation);
else
    [jobsParentPath, ~, ~] = fileparts(thisLocation);
end


%% realign
if  tasksTodo.realign
    % -------------------- realign -------------------------------------
    mkdir(alignFolder);
    jobs = {fullfile(jobsParentPath,'jobs','align_job.m')};
    %spm_jobman('serial', jobs, '', functFiles');
    spm_jobman('run', jobs{1}, functFiles);
    
    movefile(fullfile(functPath,'rf*'),alignFolder);
    movefile(fullfile(functPath,'rp*'),alignFolder);
    movefile(fullfile(functPath,'mean*'),alignFolder);
    
    copyfile(jobs{:}, jobFolder);
    
    
    %---------------- update tasksTodo -------------------------
    tasksDone.realign = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% QC
if tasksTodo.QC

    % -------------------- load ----------------------------------------
    transData = load(fullfile(alignFolder, ['rp_' functFilenames(1,1:end-3) 'txt']));
    
    functSize = size(spm_read_vols(spm_vol(functFiles{1})));
    intens = zeros(size(functFiles,1),1);
    intensTop100 = zeros(size(functFiles,1),1);
    angle = zeros(size(transData,1),1);
    transl = zeros(size(transData,1),1);
    maxTransl = zeros(size(functFiles,1),1);
    corners = cell(8,1);
    
    for k = 0:7
        mask = dec2bin(k,3) == '1';
        corners{k+1} = (mask .* functSize)';
    end
    
    for i = 1:size(functFiles,1)
        transMatrix = spm_matrix(transData(i,:));
        quat = dcm2qua(transMatrix(1:3,1:3));
        angle(i) = 2*acos(quat(1))/(2*pi)*360;
        transl(i) = norm(transMatrix(1:3,4));
        for k = 1:8
            dist = norm(corners{k} - transMatrix(1:3,:) * [corners{k};1]);
            if dist > maxTransl(i), maxTransl(i) = dist; end
        end
        
        image = spm_read_vols(spm_vol(functFiles{i}));
        intens(i) = mean(image(:));
        array = sort(image(:),'descend');
        intensTop100(i) = mean(array(1:100));
    end    
   
    intens = (intens-mean(intens))/std(intens);
    intensTop100 = (intensTop100-mean(intensTop100))/std(intensTop100); 
    
    % ---------- compute outliers --------------------------------------
    x = (1:size(intens,1))';
    c = polyfit(x,intens,2);
    f = polyval(c,x);
    
    % ---------  on mean intensity -------------------------------------
    sortIntens = sort(intens);
    lQ = sortIntens(round(0.25*size(intens,1)));
    uQ = sortIntens(round(0.75*size(intens,1)));
    IQR = uQ - lQ;
    min = f + lQ - QCcoef*IQR;
    max = f + uQ + QCcoef*IQR;
    
    indexInt = [find(intens > max); find(intens < min)];
    
%     figure;
%     plot(intens,'b'); hold on;
%     plot(max,'g');
%     plot(min,'g');
%     plot(indexInt,intens(indexInt),'ro');hold off;
    
%    ---------- compute outliers --------------------------------------
%     x = (1:size(intensTop100,1))';
%     c = polyfit(x,intensTop100,2);
%     f = polyval(c,x);
    
    %----------- on mean intensity of top 100 voxels ------------------
%     sortIntensTop100 = sort(intensTop100);
%     lQ = sortIntensTop100(round(0.25*size(intensTop100,1)));
%     uQ = sortIntensTop100(round(0.75*size(intensTop100,1)));
%     IQR = uQ - lQ;
%     min = f + lQ - 1.5*IQR;
%     max = f + uQ + 1.5*IQR;
%     
%     indexTop100 = [find(intensTop100 > max); find(intensTop100 < min)];
%     
%     
%     
%     
%     
%     figure;
%     plot(intensTop100,'b'); hold on;
%     plot(max,'g');
%     plot(min,'g');
%     plot(indexTop100,intensTop100(indexTop100),'ro');hold off;
%     
    % ------------ apply ---------------------------------------------
    mask = zeros(size(intens,1),1);
    mask(indexInt) = 1;
    
    if highpass
        mkdir(QCFolder);
        alignFilenames = spm_select('List',alignFolder,['^r.*\.' volExt '$']);
        alignFiles=cell(size(alignFilenames,1),1);
        QCFiles=cell(size(alignFilenames,1),1);
        for f=1:size(alignFilenames,1)
            alignFiles{f}=fullfile(alignFolder,alignFilenames(f,:));
            QCFiles{f}=fullfile(QCFolder,alignFilenames(f,:));
        end

        for i = size(functFiles,1):-1:1
           if mask(i) == 0, copyfile(alignFiles{i},QCFolder);
           else QCFiles(i) = []; end
        end

        highpassFilter(char(QCFiles), 2, 0);

        for f=1:size(QCFiles,1)
            delete(QCFiles{f});
        end
    end
    
    
    
    % -------------- save artifacts ----------------------------------
    artifacts = struct('rotation', {angle}, 'translation', {transl},'maxTranslation',{maxTransl},...
        'intensity', {intens}, 'intensityTop100', {intensTop100},'mask',mask);
    save(fullfile(alignFolder, 'artifacts'), 'artifacts');
    
    %---------------- update tasksTodo -------------------------
    tasksDone.QC = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% coregister: estimate
if  tasksTodo.coregister
    jobs = {fullfile(jobsParentPath,'jobs','coreg_job.m')};
    inputs = cell(3, 1);
    if strcmp(coregDirection,'StoF')
        % coreg struct to func (change headers of struct file)
        inputs{1} = {fullfile(alignFolder, ['mean' functFilenames(1,:)])};
        inputs{2} =  cellstr(structFile); % Coreg: Estimate: Source Image - cfg_files
        inputs{3} = {''};
    else
        % coreg funct to struct (change headers of func files)
        inputs{1} =  cellstr(structFile); % Coreg: Estimate: target Image
        inputs{2} = {fullfile(alignFolder, ['mean' functFilenames(1,:)])};
        tmp_others=cellstr(spm_select('List',alignFolder,['^rf.*\.' volExt '$']));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
    end
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregister = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% coregister: estimate and and reslice
if  tasksTodo.coregReslice
    warning(preprocess:untested12,'CoregReslice Not tested with SPM12');
    jobs = {fullfile(jobsParentPath,'jobs','coregReslice_job.m')};
    inputs = cell(4, 1);
    if strcmp(coregDirection,'FtoS')
        inputs{1} = cellstr(structFile);  % ref (fixed) -> struct
        inputs{2} = {fullfile(alignFolder, ['mean' functFilenames(1,:)])}; % source (moved) -> func
        tmp_others=cellstr(spm_select('List',alignFolder,['^rf.*\.' volExt '$']));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
        inputs{4} = 'r'; % prefix
    elseif strcmp(coregDirection,'FtoF')
        inputs{1} = cellstr(refMeanFunc);  % ref (fixed) -> meanfunc
        inputs{2} = {fullfile(alignFolder, ['mean' functFilenames(1,:)])}; % source (moved) -> func
        tmp_others=cellstr(spm_select('List',alignFolder,['^rf.*\.' volExt '$']));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
        inputs{4} = 'r'; % prefix
    else
        error(['unknown coregDirection: ' coregDirection]);
    end
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregReslice = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% smooth
if  tasksTodo.smooth
    warning(preprocess:untested12,'Smooth not tested with SPM12');
    jobs = {fullfile(jobsParentPath,'jobs','smooth_job.m')};
    inputs = cell(3, 1);
    % select all funcs
    tmp_data=cellstr(spm_select('List',alignFolder,['^' smoothPrefix '.*\.' volExt '$']));
    tmp_data{end+1}=['mean' functFilenames(1,:)];
    tmp_data=cellfun(@(x) fullfile(alignFolder,x),tmp_data,'UniformOutput',false); % prepend dir
    inputs{1} = tmp_data;   % functional data
    inputs{2} = smoothFWHM; % smoothing kernel specs
    inputs{3} = ['s' num2str(smoothFWHM(1))];
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregister = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% Classic segmentation and normalization
if tasksTodo.segment == CLASSICSEGMENT
    warning(preprocess:untested12,'ClassicSegment Not tested with SPM12');
    structVol = spm_vol(structFile);
    templVol = spm_vol(fullfile(spm('Dir'),'templates','T1.nii'));
    
    % -------------- segmentation --------------------------------
    fprintf('%s','* Segmentation... ');
    mkdir(segFolder);
    defaults.segment.write.wrt_cor = 0;
    
    spm_segment(structVol,templVol,defaults.segment);
    movefile(fullfile(structPath,'c*'), segFolder);
    fprintf('%s\n','Segmentation done.');
    
    % -------------- normalization -----------------------------
    fprintf('%s','* Normalisation... ');
    mkdir(normFolder);
    
    spm_normalise(templVol, structVol);
    matfile = fullfile(structPath,[structFilename(1:end-4) '_sn.mat']);
    spm_write_sn(structVol.fname, matfile);
    
    movefile(fullfile(structPath,'w*'), normFolder);
    movefile(fullfile(structPath,'*sn.mat'), normFolder);
    
    fprintf('%s\n','Normalisation done.');
    
    % ------------ update tasksTodo -------------------------------
    tasksDone.segment = CLASSICSEGMENT;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
    
end

%% New segmentation and normalization
if tasksTodo.segment == NEWSEGMENT
    % --------- do segmentation ---------------------------------
    fprintf('%s','* Segmentation... ');
    mkdir(segFolder);
    
    jobs = {fullfile(jobsParentPath,'jobs','newSegment_job.m')};
    spm_jobman('serial', jobs, '', cellstr(structFile));
    
    movefile(fullfile(structPath,'i*'), segFolder); % move inverse warping field
    movefile(fullfile(structPath,'y*'), segFolder); % move warping field
    movefile(fullfile(structPath,'c*'), segFolder); % move tissue class files
    movefile(fullfile(structPath,'*seg8.mat'), segFolder);
    copyfile(jobs{:}, jobFolder);
    
    % ----------- normalization ----------------------------------
    fprintf('%s','* Normalization... ');
    mkdir(normFolder);
    
    jobs = {fullfile(jobsParentPath,'jobs','deform_job.m')};
    inputs = cell(3, 1);
    inputs{1} = {fullfile(segFolder,['y_' structFilename])};
    inputs{2} = {structFile};
    inputs{3} = {normFolder};
    spm_jobman('serial', jobs, '', inputs{:});
    
    copyfile(jobs{:}, jobFolder);
    
    % ------------ update tasksTodo -------------------------------
    tasksDone.segment = NEWSEGMENT;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end


%% IBASPM-style Labelling
if 	tasksTodo.label
    fprintf('%s','* Atlasing... ');
    
    segmentFiles=cellstr(spm_select('FPList',segFolder,['^c.*\.' volExt '$']));
    
    %list = spm_select('List',segFolder,['^c.*\.' volExt '$']);
    %segmentFiles = cell(size(list,1),1);
    %for i = 1:size(list,1)
    %    segmentFiles{i}= fullfile(segFolder, list(i,:));
    %end
    
    if tasksDone.segment == NEWSEGMENT, deform = fullfile(segFolder, ['iy_' structFilename]);
    elseif tasksDone.segment == CLASSICSEGMENT, deform = fullfile(normFolder, [structFilename(1:end-4) '_sn.mat']); end
    
    new_Auto_Labelling(segmentFiles, atlasFile, deform, atlasFolder, []);
    fprintf('%s\n','Atlasing done.');
    
    %---------------- update tasksTodo -------------------------
    tasksDone.label = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% end
Ttotal=etime(clock, Tstart);
disp(['** DONE PROCESSING. Total time: ' num2str(Ttotal/60,'%3.1f') ' min.']);

