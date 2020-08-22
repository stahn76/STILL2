%
%
% STILL2 RSEEG Preprocessing for 2nd pre-processing
% 
% Sangtae Ahn (sangtae_ahn@med.unc.edu)
% Frohlich Lab.
%
% first written by 12/19/2016
% 
%

clear all
close all
clc

%% Load dataset
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\toolbox\eeglab13_6_5b');
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\STILL2\RSEEG');

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
name={'_EO','_EC'};

close all

%% RSEEG (EO+EC)
% POOL = parpool('local',4); % open Parallel, use matlabpool in Matlab2013a (Killdevil)

for iSub = 2:numSubs
    subId = subStruct(iSub).name;
    session_subStruct = dir([myPath '/' subId]);
    session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
    numSessions = length(session_subStruct);
    
    for iSession = 1:4
        
        sessionId = session_subStruct(iSession).name;
        fileStructDir = dir([ subId '/' sessionId '/*fil.mff']);
        fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'EC')));
        if isempty(fileStruct)
            fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'EO')));
        end
        fileId = fileStruct.name;
        
        % load data after ICA
        EEG = pop_loadset('filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_pir.set']);
        
%         % check event information
%         if length(EEG.event) ~= num_events;
%             msgbox('== Error occured due to event information ==');
%             return;
%         end
            
%        for i=1:length(EEG.event)
%             EEG.event(i).latency=int64(EEG.event(i).latency);
%        end
%         
        
        
        EEG = pop_rmdat( EEG, {'DIN2'},[0 120] ,0);
%         size(EEG.data,2)
        EEG = eeg_checkset(EEG);
        
        nPath='D:\Dropbox (Frohlich Lab)\HumanStudies\STILL2\';
        pop_saveset(EEG,'filepath',[nPath],'filename',[fileId(1:end-4)  '_pir_EO.set']);
        
        end
        
        
        
        
        

end



