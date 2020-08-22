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

for iSub = 1
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
        
        % epoching and reject bad epohcs
        type = {'DIN2','DIN3'};
        
        RSEEG=EEG;
        
        for iType = 1:length(type)
         
            [EEG,event_indices]=pop_selectevent(RSEEG,'type',type(iType));
            
            event_indices=sort([event_indices event_indices+1]);
            
            EEG.event=EEG.event(event_indices);
            
            data=[];
            for i=1:2:length(EEG.event)-1 % 2 is to find the pairs of DIN2 and DIN1
                temp=EEG.data(:,floor(EEG.event(i).latency):floor(EEG.event(i+1).latency)); % floor to make INT
                data=[data temp];
            end
            
            EEG.data=data;
            EEG.xmin=0;
            EEG.xmax=length(EEG.data)/EEG.srate;
            EEG.times=1:length(EEG.data);
            EEG.pnts=length(EEG.data);
            EEG.event=[];
            
            % Create 2 second triggers relative to the new data
            % Add new events to a loaded dataset
            nevents = 0:2*EEG.srate:length(EEG.data);
            for index = 1 : length(nevents)
                EEG.event(index).type='cue';
                EEG.event(index).latency=nevents(index);
                EEG.event(index).init_time=nevents(index)/EEG.srate;
                EEG.event(index).init_index=index;
                EEG.event(index).codes={'gidx' 1; 'cidx' 1};
                EEG.event(index).urevent=index;
                EEG.event(index).value=[];
                EEG.event(index).duration=2*EEG.srate;
            end;
            EEG = eeg_checkset(EEG, 'eventconsistency'); % Check all events for consistency
            
            EEG = pop_epoch(EEG,{'cue'},[0 2]);  
            
            pop_eegplot(EEG, 1, 1, 1);
            input('Press any key to continue...');
            
            pop_saveset(EEG,'filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) char(name(iType)) '_pir.set']);
        end
        
        
        
        
        
    end
end



