%
%
% STILL2 RSEEG Preprocessing 2nd
%
% Sangtae Ahn (sangtae_ahn@med.unc.edu)
% Frohlich Lab.
%
% first written by 1/21/2016
%
%

clear all
close all
clc

%% Load dataset
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\toolbox\eeglab13_6_5b');
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\STILL2\ODDBA');

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
num_events=200;

close all

%% ODDBA

for iSub = 22
    subId = subStruct(iSub).name;
    session_subStruct = dir([myPath '/' subId]);
    session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
    numSessions = length(session_subStruct);
    
    for iSession = 1 : 4
        
        sessionId = session_subStruct(iSession).name;
        fileStructDir = dir([ subId '/' sessionId '/*fil.mff']);
        fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'ODDBA')));
        fileId = fileStruct.name;
        
        % load data after ICA
        EEG = pop_loadset('filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_pi.set']);
        
        % check event information
        if length(EEG.event) ~= num_events;
            msgbox('== Error occured due to event information ==');
            return;
        end
        
        % select components manually
        pop_selectcomps(EEG, [1:EEG.rank] );
        pop_eegplot( EEG, 0, 1, 1);
        input('Press any key to continue...');
        
        % find rejected components and filtering the data
        ind_sel=find(EEG.reject.gcompreject == 1); % find the index for rejected
        EEG = pop_subcomp( EEG, ind_sel, 0);
        pop_saveset(EEG,'filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_pir.set']);
        EEG = eeg_checkset( EEG );
        
        EEG = pop_epoch(EEG,{'DIN1' 'DIN2'},[-0.4 0.6]);
        pop_saveset(EEG,'filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_pire.set']);
        
        pop_eegplot(EEG, 1, 1, 1);
        input('Press any key to continue...');
        pop_saveset(EEG,'filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_pirer.set']);
        
    end
end

%%

