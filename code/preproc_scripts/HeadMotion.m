function [Excl,TemporalMask, varargout] = HeadMotion(pathRP,removeFirstNscans,FDthreshold,FDthreshold2)
%Check Head Motion
% adapted from DPARSFA_run.m

% OUTPUT
% Excl: 
% temporal_mask: nVolumes x 2. First col is frame displacement computed
%   with Power's instantaneous motion on a sphere, thresholded at 0.5.
%   Volumes to discard are marked 0. Second col is FD with Van Dijk's method
%   thresholded at 0.1. Volumes to discard according to this are marked 0.
% if third output is required, the HeadMotion vector will be output: 
%   MaxRP (1:6): vector, maximum absolute motion over x,y,z translation and 
%       alpha,beta, gamma rotation
%   MeanRP (7:12): vector, mean absolute motion over x,y,z translation and 
%       alpha,beta, gamma rotation
%   MeanRMS (13): scalar, temporal average RMS translation in mm
%   MeanFD_VanDijk (14): scalar, mean frame displacement computed with Van
%       Dijk's RMS translation
%   MeanFD_Power (15): scalar, mean frame displacement computed with
%   Power's displacement on surface of a sphere. Usual threshold 0.5
%   NumberFD_Power_05 (16): integer scalar, number of frames over threshold
%   PercentFD_Power_05 (17): scalar 0-1, percentage of frames over
%       threshold
%   NumberFD_Power_02 (18): integer scalar, number of frames over threshold
%   PercentFD_Power_02 (19): scalar, percentage of frames over threshold

if nargin < 3
 FDthreshold=0.5; % Power
end
if nargin < 4
 FDthreshold2=0.1; % Van Dijk, micro movements
end

% max(abs(Tx)), max(abs(Ty)), max(abs(Tz)), max(abs(Rx)), max(abs(Ry)), max(abs(Rz)),
% mean(abs(Tx)), mean(abs(Ty)), mean(abs(Tz)), mean(abs(Rx)), mean(abs(Ry)), mean(abs(Rz)),
% mean RMS, mean relative RMS (mean FD_VanDijk), 
% mean FD_Power, Number of FD_Power>0.5, Percent of FD_Power>0.5, Number of FD_Power>0.2, Percent of FD_Power>0.2

% add subject path
mot = dir(fullfile(pathRP,'rp*.txt')); deg=180/pi; % SPM
if isempty(mot), mot=dir(fullfile(pathRP,'*.par')); deg=1; end % FSL
%RP=load('./projects/MS_dynamics/data/rest_mc.1D');
RP = load(fullfile(pathRP,mot(1).name));
RP=RP(removeFirstNscans+1:end,:);
% FSL: mm mm mm deg deg deg (counter clock wise)
% SPM: mm mm mm rad rad rad
RP(:,4:6)=RP(:,4:6)*deg; % FSL already in degrees, SPM not

% maximum motion
MaxRP = max(abs(RP));
% mean motion
MeanRP = mean(abs(RP));

%Calculate FD Van Dijk (Van Dijk, K.R., Sabuncu, M.R., Buckner, R.L., 2012. The influence of head motion on intrinsic functional connectivity MRI. Neuroimage 59, 431-438.)
RPRMS = sqrt(sum(RP(:,1:3).^2,2)); % root mean squared displacement in 3D
RPRMS2 = sqrt(sum(RP(:,4:6).^2,2));
%max(RPRMS)
%max(RPRMS2)
MeanRMS = mean(RPRMS);

% number of displacements > .1 mm in 3D btw adjacent vols
FD_VanDijk = abs(diff(RPRMS));
FD_VanDijk = [0;FD_VanDijk];
MeanFD_VanDijk = mean(FD_VanDijk);
% test=diff(RP(:,1:3)); test=[0 0 0; test]; mm=sqrt(sum(test.^2,2));

%Calculate FD Power (Power, J.D., Barnes, K.A., Snyder, A.Z., Schlaggar, B.L., Petersen, S.E., 2012. Spurious but systematic correlations in functional connectivity MRI networks arise from subject motion. Neuroimage 59, 2142-2154.) 
RPDiff=diff(RP); 
RPDiff=[zeros(1,6);RPDiff];
RPDiffSphere=RPDiff;
RPDiffSphere(:,4:6)=RPDiffSphere(:,4:6)*50*pi/180; % displacement on surface of sphere with radius 50 mm, arc length = theta(rad) * r
FD_Power=sum(abs(RPDiffSphere),2); % instantaneous head motion
MeanFD_Power = mean(FD_Power);

eval(['NumberFD_Power_' strrep(num2str(FDthreshold),'.','') ' = length(find(FD_Power>FDthreshold));']) % used by Power
eval(['PercentFD_Power_' strrep(num2str(FDthreshold),'.','') ' = length(find(FD_Power>FDthreshold)) / length(FD_Power);'])
% NumberFD_Power_02 = length(find(FD_Power>0.2));
% PercentFD_Power_02 = length(find(FD_Power>0.2)) / length(FD_Power);

HeadMotion = [MaxRP,MeanRP,MeanRMS,MeanFD_VanDijk,MeanFD_Power,...
    eval(['NumberFD_Power_' strrep(num2str(FDthreshold),'.','')]),...
    eval(['PercentFD_Power_' strrep(num2str(FDthreshold),'.','')])];%,NumberFD_Power_05,PercentFD_Power_05,NumberFD_Power_02,PercentFD_Power_02];
disp(['Mean FD_Power ' num2str(mean(FD_Power),2) ', mean FD_VanDijk ' num2str(mean(FD_VanDijk),2)])

%Write the Head Motion as .csv
% fid = fopen([AutoDataProcessParameter.DataProcessDir,filesep,'RealignParameter',filesep,FunSessionPrefixSet{iFunSession},'HeadMotion.csv'],'w');
% fprintf(fid,'Subject ID\tmax(abs(Tx))\tmax(abs(Ty))\tmax(abs(Tz))\tmax(abs(Rx))\tmax(abs(Ry))\tmax(abs(Rz))\tmean(abs(Tx))\tmean(abs(Ty))\tmean(abs(Tz))\tmean(abs(Rx))\tmean(abs(Ry))\tmean(abs(Rz))\tmean RMS\tmean relative RMS (mean FD_VanDijk)\tmean FD_Power\tNumber of FD_Power>0.5\tPercent of FD_Power>0.5\tNumber of FD_Power>0.2\tPercent of FD_Power>0.2\n');
% for i=1:AutoDataProcessParameter.SubjectNum
%     fprintf(fid,'%s\t',AutoDataProcessParameter.SubjectID{i});
%     fprintf(fid,'%e\t',HeadMotion(i,:));
%     fprintf(fid,'\n');
% end
% fclose(fid);

Range=3:-0.25:0;
tmp=repmat(HeadMotion(1:6),length(Range),1);
Ex=repmat(Range',1,6);
Excl=find(any(tmp>Ex,2));
if ~isempty(Excl), 
    Mot=find(tmp(Excl(1),:)>Ex(Excl(1),:)); % index in HeadMotion
    if Mot(1)>3, unit='deg'; else unit='mm'; end
    disp(['Max head motion > ' num2str(Range(Excl(1))) ' ' unit ':  ' num2str(HeadMotion(1:6),2)]);
end

% for i=1:length(Range)
%     ExcludingCriteria=Range(i);
%     BigHeadMotion=find(HeadMotion(1:6)>ExcludingCriteria, 1);
%     if ~isempty(BigHeadMotion)
%         Excl(i)=1;
%         sprintf('\nExclusion Criteria: %2.1fmm and %2.1f degree in max head motion\n%s\n\n\n',ExcludingCriteria,ExcludingCriteria);
%     end
% end

TemporalMask=ones(size(RP,1),2);
Index= FD_Power > FDthreshold;
TemporalMask(Index,1)=0;
Index= FD_VanDijk > FDthreshold2;
TemporalMask(Index,2)=0;

if nargout==3
    varargout{1}=HeadMotion;
end

end
