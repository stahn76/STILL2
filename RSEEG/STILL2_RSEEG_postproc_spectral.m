%
%
% STILL2 RSEEG postproc spectral
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

% outside of the scalp
rmv_ch=[114 121 1 8 14 21 25 32 38 44 ...
    57 64 69 74 82 89 95 100 ...
    120 113 107 99 94 88 81 73 68 63 56 49 43 ...
    48 119 125 128 17 126 127];



ind=[0 2 4 6];

psd.delta=[];
psd.theta=[];
psd.alpha=[];
psd.beta=[];
psd.gamma=[];

sub_psd=[];

close all

%% RSEEG (EO+EC)
% POOL = parpool('local',4); % open Parallel, use matlabpool in Matlab2013a (Killdevil)

for iSub = 13
    subId = subStruct(iSub).name;
    session_subStruct = dir([myPath '/' subId]);
    session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
    numSessions = length(session_subStruct);
    
    figure;
    for iSession = 1  : 4
%         numSessions
        
        sessionId = session_subStruct(iSession).name;
        fileStructDir = dir([ subId '/' sessionId '/*fil.mff']);
        fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'EC')));
        if isempty(fileStruct)
            fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'EO')));
        end
        fileId = fileStruct.name;
        
        for iType = 1 : 2
%             : length(name)
            
            EEG = pop_loadset('filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_' char(name(iType)) '_pir.set']);
            ch=setdiff(EEG.chaninfo.icachansind,rmv_ch);
            %             figure;
            
            [spec freq]=spectopo(EEG.data,0,EEG.srate,'plot','off');
            
%             delta = find(freq>=1 & freq<=4);
%             theta = find(freq>=4 & freq<=8);
            alpha = find(freq>=8 & freq<=12);
%             beta  = find(freq>=13 & freq<=30);
%             gamma = find(freq>=30 & freq<=50);
%             all   = find(freq>=1 & freq<=50);
            
            sub_psd=[sub_psd mean(spec(ch,alpha),2)];
            
            subplot(numSessions,length(name),ind(iSession)+iType);
            
%             psd.delta=[psd.delta mean(spec(ch,delta),2)];
%             psd.theta=[psd.theta mean(spec(ch,theta),2)];
%             psd.alpha=[psd.alpha mean(spec(ch,alpha),2)];
%             psd.beta =[psd.beta mean(spec(ch,beta),2)];
%             psd.gamma=[psd.gamma mean(spec(ch,gamma),2)];
            
                        title([char(name(iType)) ' - S' num2str(iSession)]);
                        topoplot(mean(spec(ch,alpha),2),EEG.chanlocs(ch));
%                         topoplot(sum(spec(ch,alpha),2)./sum(spec(ch,all),2),EEG.chanlocs(ch));
                        colorbar;
            
            
        end
        
    end
end



