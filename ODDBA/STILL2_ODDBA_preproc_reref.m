

clear all;
close all;
clc;


%%
eeglab;

myPath='D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\Data\STILL2\EGI\';
cd(myPath);
subStruct = dir;
subStruct = subStruct(cellfun(@any,strfind({subStruct.name},'P_0')));
numSubs = length(subStruct);
session_subStruct = dir([myPath '/' subStruct(1).name]);
session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
nSessions = length(session_subStruct);

% initial parameters
desiredFs = 250;
lowPass = 1;
highPass = 50;

%%

for iSub = 1;
    subId = subStruct(iSub).name;
    session_subStruct = dir([myPath '/' subId]);
    session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
    numSessions = length(session_subStruct);
    
    for iSession = 2
        
        sessionId = session_subStruct(iSession).name;
        fileStructDir = dir([ subId '/' sessionId '/*fil.mff']);
        fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'ODDBA')));
        fileId = fileStruct.name;
        
        % load data
        disp(['Loading the data ' fileId '...']);
        EEG = pop_readegimff([myPath subId '\' sessionId '\' fileId],'datatype','EEG');
        EEG.event = EEG.event(2:end); % remove trash trigger
        
        % resampling
        EEG = pop_resample (EEG, desiredFs);
        
        % filtering
        EEG = pop_eegfiltnew(EEG, 1, 50, 1650, 0, [], 0);
        
        % clean line
        %         EEG = pop_cleanline(EEG, 'bandwidth', 2,'chanlist', [1:EEG.nbchan], 'computepower', 0, 'linefreqs', [60],...
        %             'normSpectrum', 0, 'p', 0.01, 'pad', 2, 'plotfigures', 0, 'scanforlines', 1, 'sigtype', 'Channels', 'tau', 100,...
        %             'verb', 1, 'winsize', 4, 'winstep', 4);
        
        % ASR
        EEG.originalEEG=EEG; % keep origianl EEG before ASR
        EEG = clean_rawdata(EEG,5,-1,0.8,4,20,-1);
        EEG = eeg_checkset( EEG );
        
        % Interpolate bad channel
        
        EEG = pop_interp(EEG, EEG.originalEEG.chanlocs, 'spherical');
        
        % re-referencing
        EEG.nbchan = EEG.nbchan+1;
        EEG.data(end+1,:) = zeros(1, EEG.pnts);
        EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
        EEG = pop_reref(EEG, []);
        EEG = pop_select( EEG,'nochannel',{'initialReference'});
        
        
    end
end


%%

% outside of the scalp
rmv_ch=[114 121 1 8 14 21 25 32 38 44 ...
    57 64 69 74 82 89 95 100 ...
    120 113 107 99 94 88 81 73 68 63 56 49 43 ...
    48 119 125 128 17 126 127];
ch=setdiff([1:128],rmv_ch);

% eeglab;
EEG = pop_loadset('filename','test_ICAr.set','filepath','D:\\Dropbox (Frohlich Lab)\\Sangtae\\MATLAB\\Data\\STILL2\\EGI\\P_001\\S1\\');

temp=EEG;

EEG = pop_epoch( temp, {'DIN1'}, [-0.5 1], 'epochinfo', 'yes');
% eeglab redraw


EEG = pop_rmbase( EEG, [-200    0]);
% figure;
% pop_timtopo(EEG, [-500  996], [NaN]);
figure;
pop_plottopo(EEG, ch , 'ERP on scalp channels', 0, 'ydir',1);
figure; 
pop_erpimage(EEG,1, [55],[[]],'E55',10,1,{},[],'' ,'yerplabel','\muV','erp','on','cbar','on','topo', { [55] EEG.chanlocs EEG.chaninfo } );


% figure;
% imagesc(squeeze(EEG.data(55,:,:))');
% colorbar;


% infomax ICA
if isfield(EEG.etc, 'clean_channel_mask')
    EEG.dataRank = min([rank(double(EEG.data')) sum(EEG.etc.clean_channel_mask)]);
else
    EEG.dataRank = min(rank(double(EEG.data')));
end

EEG = pop_runica(EEG,'extended',1,'pca',EEG.dataRank);


EEG=pop_eegplot(EEG, 1, 1, 1);







