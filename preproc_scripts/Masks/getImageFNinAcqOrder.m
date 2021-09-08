function fList=getImageFNinAcqOrder(fpath,myP,myE,scanIdx,varargin)
% Return a list of files found in fpath in the order
% of the acquisition sequence.
% 
%
% IN    fpath:  path to where the functional data are stored
%       myP:    prefix of the files to list (typically according to SPM
%           convention - e.g. 'rf'
%       myE:    extension of the file to list, typically 'nii'
%       scanIdx: vector of scan indices to select, [] if all
% 
% ASSUMPTIONS:
% Works only for CHUV-rest or EMOrest resting-state data with
% volume files naming convention
% rfAAAAAAA-BBBBBB-CCCCC-CCCCC-1.nii
% where:
% AAAAAAA:  session / user (?) identifier
% BBBBBB:   scan time, format HHMMSS
% CCCCC:    volume sequence number
% e.g. 'rf2664514-163313-00001-00001-1.nii'
% or 'rfIRMF19870224FROAU-0009-00257-000257-01.nii'
%
% v1.0 Sep 2009 Jonas Richiardi
% - initial release, deals with CHUV-rest (MS) and EMOrest databases
% v1.1 Jul 2010 Jonas Richiardi
% - more flexible, also deals with CMSTrest type DB with inconsistent
% number of fields (missing BBBBBB in schema above)
% v1.2 May 2011 JR
% - still more flexible: reverts to reading from disk order if cannot
% parse file format.
% v1.2.1 Dec 2011 JR
% - support for all-dash format
% v1.3 Jan 2012 JR
% - support for skipping parsing

tryParsing=process_options(varargin,'tryParsing',true);

% get a dir listing corresponding to desiderata
dirL=dir(fullfile(fpath,[ myP '*.' myE]));

% check whether we have a 4-field or 5-field file name
foo=dirL(1).name;
dashIdx=strfind(foo,'-');
underIdx=strfind(foo,'_');
nFields_dash=numel(dashIdx)+1;
nFields_under=numel(underIdx)+1;

if tryParsing==true % check if we can really do it
    if nFields_dash==1 && nFields_under==1
        error('No dashes or underscores in filename, cannot parse');
    elseif nFields_dash>1 && nFields_under>1
        warning('getImageFNinAcqOrder:GivingUpParsing',['! Both dashes and ' ...
            ' underscores in filename, reverting to on-disk order. This ' ...
            'ordering may be OS dependent']);
        tryParsing=false;
    else
        tryParsing=true;
    end
end

if tryParsing==true
    if nFields_dash>1
        parseDashes=true;
        parseUnder=false;
        if nFields_dash==5
            parseFmtString='%s%s%n%n%s';
        elseif nFields_dash==4
            parseFmtString='%s%n%n%s';
        else
            error(['Unknown file name format, cannot guess sequence: ' foo]);
        end
    elseif nFields_under>1
        parseUnder=true;
        parseDashes=false;
        if nFields_under==6
            parseFmtString='%s%s%s%s%n%n';
        elseif nFields_under==3
            parseFmtString='%s%s%n';
        else
            error(['Unknown file name format, cannot guess sequence: ' foo]);
        end
    end
end


% build a numeric array of sequence numbers
nFiles=numel(dirL);
fList=cell(nFiles,1);
if tryParsing==true
    sn=zeros(nFiles,1); % number of images in the acquisition sequence
    if parseDashes
        for f=1:nFiles
            fList{f}=dirL(f).name;
            foo=textscan(fList{f},parseFmtString,'Delimiter','-');
            sn(f)=foo{3};
        end
    elseif parseUnder
        for f=1:nFiles
            fList{f}=dirL(f).name;
            foo=textscan(fList{f},parseFmtString,'Delimiter','_');
            sn(f)=foo{end};
        end
    end
    % sort list according to sequence number
    [sn_srt,sn_srtIdx]=sort(sn,1,'ascend');
    fList=fList(sn_srtIdx);
else
    for f=1:nFiles
        fList{f}=dirL(f).name;
    end
end

% subset if needed
if ~isempty(scanIdx)
    % check order of scan indices
    if (~all(sort(scanIdx)==scanIdx))
        warning('getImageFNinAcqOrder:weirdScanIdx',...
            ['Scan indices are not in increasing order, make sure you know' ...
            ' what you are doing']);
    end
    fList=fList(scanIdx);
end