function [TSscrub,indd] = scrubTimeSeries(TSblock,artifactsFile,funcPath,removeFirstNscans,ScrubbingMethod,FDthreshold,debug)

% Adapted from Nora Leonardi's timeCourseExtraction script

if ~exist('FDthreshold','var')
    FDthreshold = 0.5;
end

%% Motion estimates
fprintf('Checking motion parameters for excessive motion...');
[~,TemporalMask] = HeadMotion(funcPath,removeFirstNscans,FDthreshold);
load(artifactsFile); 
intens = artifacts.intensity(removeFirstNscans+1:end);

DVARS = abs(diff(intens)); DVARS=[0;DVARS]; % in principle RMS on voxels, excluding those with no full coverage over scan
    
% or apply later
% Power chose .5% BOLD change
y = intens'; x = 1:length(y); 
[p,err] = polyfit(x,y,2);
[y_fit,delta] = polyval(p,x,err);
TooLarge = zeros(size(x));
TooLarge(y>y_fit+3*delta) = 1;
TooLarge(y<y_fit-3*delta) = 1;
idxOut = find(TooLarge);

b = 1;
scansIdx{b} = 1:length(y); % variable from block selection, here not used
c = ismember(scansIdx{b},idxOut);
TemporalMask = [TemporalMask,ones(size(TemporalMask,1),1)];
TemporalMask(scansIdx{b}(c),3)=0;
TM = all(TemporalMask(:,[1,3])'); % any of 3 indices says to exclude (VanDijk currently not used)
% exclude previous 1 + after 1
indd = find(TM==0);
% disp([num2str(length(indd)) ' scans to exclude (adding 1 before, 1 after): ' num2str(indd)]);
% indd=[indd,indd-1,indd+1]; indd=unique(indd);
disp([num2str(length(indd)) ' scans to exclude (adding 1 before, 2 after): ' num2str(indd)]);
indd = [indd,indd-1,indd+1,indd+2]; indd=unique(indd);
if any(diff(indd)==2), % ???? would need diff twice?
    fprintf('Adding scan %d because single scan in between excluded scans\n', num2str(indd(diff(indd)==2)+1));
    indd = [indd,indd(diff(indd)==2)+1]; 
end % one scan remains between excluded scans
valid = (indd>=1 & indd<=length(TM));
indd = indd(valid);
TM(indd)=0;

save(fullfile(funcPath,['ExcMotion' strrep(num2str(FDthreshold),'.','') '.mat']), 'indd');

if exist('debug','var') && debug
    figure;plot(artifacts.intensity); hold on; plot(indd+5,artifacts.intensity(indd+5),'xr'); 
    plot(y_fit,'g');plot(y_fit+3*delta,'g--');plot(y_fit-3*delta,'g--');
    title('artifacts.intensity fit and all time points to exclude');
    savefig(fullfile(funcPath,['ExcMotion' strrep(num2str(FDthreshold),'.','')]));
    close gcf;
end


%% Scrubbing
% ScrubbingMethod - cut, nearest, linear, spline, pchip
if ~all(TM) && ~isempty(ScrubbingMethod) % then at least 1 to exclude
        fprintf('Scrubbing functional data, using %s interpolation\n',ScrubbingMethod);
        TSscrub = TSblock(:,TM); %'cut'
        if ~strcmpi(ScrubbingMethod,'cut')
            xi = 1:length(TM);
            x = xi(TM);
            TSscrub = interp1(x,TSscrub',xi,ScrubbingMethod,'extrap')'; % extrapolate since otherwise returns NaN if x value out of range
            % yi=interp1(x,Y,xi): find yi, fct Y, at points xi               
        end
else
    TSscrub = TSblock;
end