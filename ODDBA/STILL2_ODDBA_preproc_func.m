function [ outEEG ] = STILL2_ODDBA_preproc_func( inEEG , desiredFs, lowCut, highCut)

%
% function outEEG = func_preproc( inEEG , desiredFs, lowCut, highCut)
%
%  1. downsampling based on desiredFs
%  2. band-pass filtering based on lowCut & highCut
%  3. artifact subspace reconstruction (ASR)
%  4. channel interpolation
%  5. re-referencing
%
% *** INPUT ***
%   EEG : eeglab structure
%   desiredFs [scalar] : sampling rate (Hz) to be downsampled
%   lowCut [scalar] : high-pass frequency (Hz)
%   highCut [scalar] : low-pass frequency (Hz)
%
% *** OUTPUT ***
%   outEEG : eeglab structure after preprocessing
%
% *** USEAGE ***
%
% outEEG = func_preproc( inEEG, desiredFs, lowPass, highPass )
%
% Frohlich Lab. Sangtae Ahn (sangtae_ahn@med.unc.edu)
%
% first written by 11/30/2016
% revised by 12/13/2016 : change parameters as GUI default setting in clean_rawdata
%

%% Resampling and band-pass filtering
EEG = pop_resample( inEEG, desiredFs);
disp(['resampled with : ' num2str(EEG.srate) 'Hz']);

EEG = pop_eegfiltnew(EEG, lowCut, highCut);
disp(['band-pass filtering from '  num2str(lowCut) ' to ' num2str(highCut)  ' Hz']);

EEG = eeg_checkset( EEG );


%% Save chanloc structure for future use (Interpolatation)
EEG.historychanlocs=EEG.chanlocs; % Save channel locs of 128 because ASR will remove bad channel loc information
EEG.historychaninfo=EEG.chaninfo;

%% Run Artifact subspace reconstruction (removes bad epoch data (PCA), bad channels)
EEG = clean_rawdata(EEG,5,-1,0.8,4,20,-1); % default setting
% EEG.badchan=find(EEG.etc.clean_channel_mask==0); %Bad chananel information from ASR
EEG = eeg_checkset( EEG );

%% Interpolate bad channels
EEG.originalEEG=EEG; % keep origianl EEG before interpolation
EEG = pop_interp(EEG, EEG.historychanlocs, 'spherical');




%% Rereference to AVREF (should be with interpolated channels included!)
EEG = eeg_checkset( EEG );

EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:) = zeros(1, EEG.pnts);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, []);
EEG = pop_select( EEG,'nochannel',{'initialReference'});


outEEG=EEG;


end

