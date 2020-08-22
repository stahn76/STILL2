%
%
% STILL2 RSEEG postproc PSD
%
% Sangtae Ahn (sangtae_ahn@med.unc.edu)
% Frohlich Lab.
%
% first written by 1/23/2017
% revised by 2/6/2017 : changed psd method pmtm -> pwelch
% 

clear all
close all
clc

%% Load dataset
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\toolbox\eeglab13_6_5b');
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\STILL2');

eeglab;
pop_editoptions( 'option_savetwofiles', 1,'option_single', 0);

myPath='D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\Data\STILL2\EGI\';
cd(myPath);
subStruct = dir;
subStruct = subStruct(cellfun(@any,strfind({subStruct.name},'P_0')));
numSubs = length(subStruct);
session_subStruct = dir([myPath '/' subStruct(1).name]);
session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
nSessions = length(session_subStruct);

% initial parameters
num_events=8;
name={'EO','EC'};


close all
clear EEG;

EO=cell(numSubs,4);
EC=cell(numSubs,4);
iEO=1;
iEC=2;
PxxEO=cell(numSubs,4);
PxxEC=cell(numSubs,4);
srate=250;


%% RSEEG (EO+EC)


for iSub = 1
    subId = subStruct(iSub).name;
    session_subStruct = dir([myPath '/' subId]);
    session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
    numSessions = length(session_subStruct);
    
    parfor iSession =2
        
        sessionId{iSession} = session_subStruct(iSession).name;
        fileStructDir{iSession} = dir([ subId '/' sessionId{iSession} '/*fil.mff']);
        fileStruct{iSession} = fileStructDir{iSession}(cellfun(@any,strfind({fileStructDir{iSession}.name},'EC')));
        if isempty(fileStruct{iSession})
            fileStruct{iSession} = fileStructDir{iSession}(cellfun(@any,strfind({fileStructDir{iSession}.name},'EO')));
        end
        fileId{iSession} = fileStruct{iSession}.name;
  
        
        % Power Spectral Density (PSD) estimate via the Thomson multitaper
        % or
        % Power Spectral Density (PSD) estimate via Welch's method
        
        
        % EO
        EO{iSub,iSession} = pop_loadset('filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '_' char(name(iEO)) '_pir.set']);
        
        
        for iEpoch = 1:size(EO{iSub,iSession}.data,3)
            [PxxEO{iSub,iSession}(:,:,iEpoch)] = pmtm(EO{iSub,iSession}.data(:,:,iEpoch)',1.25,[1:0.5:50],EO{iSub,iSession}.srate);
%             [PxxEO{iSub,iSession}(:,:,iEpoch)] = pwelch(EO{iSub,iSession}.data(:,:,iEpoch)',srate*2,srate/4,srate*2,srate);
        end
        PxxEO{iSub,iSession}=mean(PxxEO{iSub,iSession},3);
%         
%         figure;
%         plot(10*log(PxxEO{1,2}(:,17)));
        
        % EC
        EC{iSub,iSession} = pop_loadset('filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '_' char(name(iEC)) '_pir.set']);
        for iEpochEC = 1:size(EC{iSub,iSession}.data,3)
            [PxxEC{iSub,iSession}(:,:,iEpochEC)] = pwelch(EC{iSub,iSession}.data(:,:,iEpochEC)',srate*2,srate/4,srate*2,srate);
        end
        PxxEC{iSub,iSession}=mean(PxxEC{iSub,iSession},3);
%          
%         

        
        
    end
    
end


