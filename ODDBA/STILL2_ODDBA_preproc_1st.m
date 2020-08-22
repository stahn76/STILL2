%
%
% STILL2 ODDBA Preprocessing for Parallel Computing
%
% Sangtae Ahn (sangtae_ahn@med.unc.edu)
% Frohlich Lab.
%
% first written by 1/12/2017
%
%


close all
clear
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



% initial parameters
desiredFs = 250;
lowCut = 1;
highCut = 50;
nSessions=4;
close all;



%% ODDBA

POOL = parpool('local',4); % open Parallel, use matlabpool in Matlab2013a (Killdevil)

for iSub = 1
    
    % initialization for Sessions Cell
    sessionId=cell(nSessions,1);
    fileStructDir=cell(nSessions,1);
    fileStruct=cell(nSessions,1);
    fileId=cell(nSessions,1);
    EEG=cell(nSessions,1);
    PNS=cell(nSessions,1);
    temp=cell(nSessions,1);
    
    subId = subStruct(iSub).name;
    session_subStruct = dir([myPath '/' subId]);
    session_subStruct = session_subStruct(cellfun(@any,regexp({session_subStruct.name},'S[1-4]')));
    numSessions = length(session_subStruct);
    
    parfor iSession = 1 : 4
        
        sessionId{iSession} = session_subStruct(iSession).name;
        fileStructDir{iSession} = dir([ subId '/' sessionId{iSession} '/*fil.mff']);
        fileStruct{iSession} = fileStructDir{iSession}(cellfun(@any,strfind({fileStructDir{iSession}.name},'ODDBA')));
        fileId{iSession} = fileStruct{iSession}.name;
        
        
        %% Part 1 : load and save dataset
        disp(['Loading the data ' fileId{iSession} '...']);
        EEG{iSession} = pop_readegimff([myPath subId '\' sessionId{iSession} '\' fileId{iSession}],'datatype','PIB');
        EEG{iSession}.event = EEG{iSession}.event(2:end); % remove trash trigger
        
        disp(['The number of events : ' num2str(length(EEG{iSession}.event))]);
    
        
        PNS{iSession} = pop_readegimff([myPath subId '\' sessionId{iSession} '\' fileId{iSession}],'datatype','PIB');
        PNS{iSession}.event=PNS{iSession}.event(2:end);

        
%         thresh=10000;
%         width=100;
%         [pks locs] = findpeaks(PNS{2}.data(1:end),'MinPeakHeight',thresh,'MinPeakDistance',width);
%         locs=locs';
%         figure;
%         plot(PNS{2}.data);
  
%          parfor iSession = 1:4

        % Extract the ODDBA EEG data 
        
        % based on sind and eind variables
%         EEG{iSession}.etc.sind=EEG{iSession}.event(1).latency;
%         EEG{iSession}.etc.new_sind=EEG{iSession}.etc.sind-EEG{iSession}.srate*0.4; % extract 400ms data before 1st trigger
%         
%         EEG{iSession}.etc.eind=EEG{iSession}.event(length(EEG{iSession}.event)).latency;
%         EEG{iSession}.etc.new_eind=EEG{iSession}.etc.eind+EEG{iSession}.srate*0.6; % extract 600ms data after last trigger
%         
%         temp{iSession}=EEG{iSession};
%         EEG{iSession}.data=[];
%         EEG{iSession}.data=temp{iSession}.data(:,temp{iSession}.etc.new_sind:temp{iSession}.etc.new_eind);
%         
%         % Check the compatibility of the current EEG dataset
%         EEG{iSession}.xmin=0;
%         EEG{iSession}.xmax=length(EEG{iSession}.data)/EEG{iSession}.srate;
%         EEG{iSession}.times=1:length(EEG{iSession}.data);
%         EEG{iSession}.pnts=length(EEG{iSession}.data);
%         
%         EEG{iSession}.orig_event=EEG{iSession}.event; % save original events
%         EEG{iSession}.event=[];
%         EEG{iSession}.event=EEG{iSession}.orig_event; % extract event for EO + EC
%         EEG{iSession}.cutevent=temp{iSession}.etc.new_sind-1; % makes start flag 1
%         
%         % adjusting event structure
%         for index = 1 : length(temp{iSession}.event);
%             EEG{iSession}.event(index).latency=EEG{iSession}.orig_event(index).latency-EEG{iSession}.cutevent;
%             EEG{iSession}.event(index).init_time=EEG{iSession}.orig_event(index).latency/EEG{iSession}.srate;
%             EEG{iSession}.event(index).init_index=index;
%             EEG{iSession}.event(index).codes={'gidx' 1; 'cidx' 1};
%             EEG{iSession}.event(index).urevent=index;
%         end
%         
        EEG{iSession} = eeg_checkset(EEG{iSession}, 'eventconsistency'); % Check all events for consistency
        pop_saveset(EEG{iSession},'filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '_NoEpoch.set']);
        
        %% Part 2 : Preprocessing
        
        % resampling
        EEG{iSession} = STILL2_ODDBA_preproc_func(EEG{iSession},desiredFs,lowCut,highCut);
        pop_saveset(EEG{iSession},'filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '__NoEpoch_p.set']);
                
        
        
%         EEG{iSession} = pop_loadset('filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '_p.set']);
    
        
        %% infomax ICA
%         EEG{iSession}.rank=rank(double(EEG{iSession}.data));
%     
%         EEG{iSession} = pop_runica(EEG{iSession},'extended',1,'pca',EEG{iSession}.rank);
%         EEG{iSession} = eeg_checkset(EEG{iSession});
%         
%         pop_saveset(EEG{iSession},'filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '_pi.set']);
        
        
        
        %% dipole fit
        
%         %only label electrodes that are less than 1cm deviant from true 10-10 location
%         EEG{iSession}=pop_chanedit(EEG{iSession}, 'changefield',{11 'labels' 'Fz'},'changefield',{9 'labels' 'Fp2'},...
%             'changefield',{22 'labels' 'Fp1'},'changefield',{24 'labels' 'F3'},...
%             'changefield',{36 'labels' 'C3'},'changefield',{62 'labels' 'Pz'},'changefield',{92 'labels' 'P4'},...
%             'changefield',{104 'labels' 'C4'},...
%             'changefield',{124 'labels' 'F4'},'changefield',{3 'labels' 'AF4'},'changefield',{2 'labels' 'AF8'},...
%             'changefield',{16 'labels' 'AFz'},'changefield',{30 'labels' 'C1'},'changefield',{105 'labels' 'C2'},...
%             'changefield',{87 'labels' 'CP2'},'changefield',{42 'labels' 'CP3'},'changefield',{93 'labels' 'CP4'},...
%             'changefield',{47 'labels' 'CP5'},'changefield',{98 'labels' 'CP6'},'changefield',{55 'labels' 'CPz'},...
%             'changefield',{1 'labels' 'F10'},'changefield',{4 'labels' 'F2'},'changefield',{27 'labels' 'F5'},...
%             'changefield',{123 'labels' 'F6'},'changefield',{32 'labels' 'F9'},'changefield',{13 'labels' 'FC1'},...
%             'changefield',{29 'labels' 'FC3'},'changefield',{28 'labels' 'FC5'},'changefield',{6 'labels' 'FCz'},...
%             'changefield',{60 'labels' 'P1'},'changefield',{97 'labels' 'P6'},'changefield',{72 'labels' 'POz'});
%         
%         EEG{iSession}=pop_chanedit(EEG{iSession}, 'eval','chans = pop_chancenter( chans, [],[]);');
%         
%         % Fit single dipole
%         EEG{iSession} = pop_dipfit_settings( EEG{iSession}, 'hdmfile','D:\Sangtae\MATLAB\toolbox\eeglab13_6_5b\plugins\dipfit2.3\standard_BEM\standard_vol.mat',...
%             'coordformat','MNI','mrifile','D:\Sangtae\MATLAB\toolbox\eeglab13_6_5b\plugins\dipfit2.3\standard_BEM\standard_mri.mat',...
%             'chanfile','D:\Sangtae\MATLAB\toolbox\eeglab13_6_5b\plugins\dipfit2.3\standard_BEM\elec\standard_1005.elc',...
%             'coord_transform',[0.49027 -18.8194 -16.8183 0.022217 -0.017001 -1.5661 11.5973 11.3797 13.4455] ,'chansel',[1:128] );
%         EEG{iSession} = pop_multifit(EEG{iSession}, 1:EEG{iSession}.rank ,'threshold',100,'plotopt',{'normlen' 'on'});
%         
%         % Fit bilateral dipole
%         EEG{iSession} = fitTwoDipoles(EEG{iSession}, 'LRR', 35);
%         
%         pop_saveset(EEG{iSession},'filepath',[subId '/' sessionId{iSession}],'filename',[fileId{iSession}(1:end-4) '_pid.set']);
        
    end
end

msgbox('done');
delete(POOL);




