%
%
% STILL2 CLICK postproc tf
%
% Sangtae Ahn (sangtae_ahn@med.unc.edu)
% Frohlich Lab.
%
% first written by 03/27/2017
%
%


close all
clear
clc

%% Load dataset
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\toolbox\eeglab13_6_5b');
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\toolbox\hline_vline');
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\toolbox\DrosteEffect-BrewerMap-54c4241');
addpath('D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\STILL2\CLICK');


eeglab;

pop_editoptions( 'option_savetwofiles', 1,'option_single', 0);

myPath='D:\Dropbox (Frohlich Lab)\Sangtae\MATLAB\Data\STILL2\EGI\';
% myPath='D:\Sangtae\MATLAB\Data\STILL2\EGI';
cd(myPath);
subStruct = dir;
subStruct = subStruct(cellfun(@any,strfind({subStruct.name},'P_0')));
nSub = length(subStruct);
session_subStruct = dir([myPath '/' subStruct(1).name]);
session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
nSes = length(session_subStruct);


% outside of the scalp
rmv_ch=sort([114 121 1 8 14 21 25 32 38 44 ...
    57 64 69 74 82 89 95 100 ...
    120 113 107 99 94 88 81 73 68 63 56 49 43 ...
    48 119 125 128 17 126 127]);

% inside of the scalp
ch=sort(setdiff(1:128,rmv_ch));

% based on regions
occipitalScalp=sort([59,60,61,65,66,67,70,71,62,72,75,76,77,78,83,84,85,90,91]);
LtemporalScalp=sort([28,33,34,35,40,41,45,46,47,50,51,52,58]);
RtemporalScalp=sort([92,96,97,98,101,102,103,108,109,110,115,116,117,122]);
frontalScalp=sort([2,3,4,5,123,124,9,10,11,16,15,12,19,18,24,23,22,26,27]);
centralScalp=sort([20,13,7,29,30,36,31,37,42,53,54,55,6,106,112,118,80,105,111,104,87,93,79,86]);
allScalp=sort([frontalScalp centralScalp LtemporalScalp RtemporalScalp occipitalScalp]);

% define labeling
tACS=sort([4 8 9 13 14 16 17 19]);
tDCS=sort([2 3 5 7 10 12 21]);
sham=sort([1 6 11 15 18 20 22]);

Arms={tACS,tDCS,sham};
ArmsName=['tACS';'tDCS';'sham'];
nArm=3;
nTrial=200;
srate=250;

DINs={'DIN1';'DIN2';'DIN3';'DIN4';'DIN5'};
Stim={'10Hz','20Hz','30Hz','40Hz','80Hz'};
locs=readlocs('GSN-HydroCel-128_el.sfp');

maxFreq=100;
freqBin=1;

fromTime=-0.1;
toTime=0.6;

frange=1:freqBin:maxFreq;
time=linspace(fromTime*1000,toTime*1000,(toTime-fromTime)*srate);
freq=linspace(1,maxFreq,length(frange));


close all;

%% CLICK

for iDIN = 4
    
    
    for iSub = 2 : nSub
        
        subId = subStruct(iSub).name;
        session_subStruct = dir([myPath '/' subId]);
        session_subStruct = session_subStruct(cellfun(@any,strfind({session_subStruct.name},'S')));
        numSessions = length(session_subStruct);
        
        %     figure;
        
        for iSession = 1 : 4
            
            sessionId = session_subStruct(iSession).name;
            fileStructDir = dir([ subId '/' sessionId '/*fil.mff']);
            fileStruct = fileStructDir(cellfun(@any,strfind({fileStructDir.name},'CLICK')));
            fileId = fileStruct.name;
            
            EEG = pop_loadset('filepath',[subId '/' sessionId],'filename',[fileId(1:end-4) '_pir.set']);
            EEG = pop_epoch( EEG, {  'DIN4'  }, [-0.1         0.6]);
            
            for iCh = 1 : EEG.nbchan
                wav{iCh} = SA_wavelet(EEG.data(iCh,:,:),params);
            end
            for iCh = 1 : EEG.nbchan
                ITPC=abs(mean(exp(1i*wav{iCh}.phase),3));
                figure;imagesc(time,freq,ITPC);set(gca,'YDir','normal');
                saveas(gcf,[num2str(iCh) '.png']);
                close all;
            end
            
            % extract events
            for n=1:length(EEG.event)
                latency(n,1)=EEG.event(n).latency; % Get latency of urevents that are artifatc free
                type{n}=EEG.event(n).type;
            end
            latency=int64(latency); % make sure latency is integer
            
            latencyDIN=[];
            for i = 1 : length(EEG.event)
                if strcmp(EEG.event(i).type,DINs{iDIN})
                    latencyDIN=[latencyDIN;int64(EEG.event(i).latency)];
                else
                end
            end

            
            params.Fs=EEG.srate;
            params.frange=frange;
            wave_conv=[];
            
            % iteration for channels
            for iCh = 55
                
                % calculate wavelet for individual channel due to Memory
                wave_conv=wavelet_convolution_STILL2_CLICK(EEG.data(ch(iCh),:),params);
                
                wave_conv=wavelet_convolution_STILL2_CLICK(EEG.data,params);
                
%                 Wave=WaveletConvolution(EEG.data(55,:,:),params);
%                 ITPC=abs(squeeze(mean(exp(1i*Wave.phase),4)));
%                 figure;
%                 imagesc(time,freq,ITPC);
%                 colorbar;caxis([0 0.3]);
                
                
                
                wave_conv.amplitude=squeeze(wave_conv.amplitude);
                
                
                wave_conv.phase=squeeze(wave_conv.phase);
                
                
                
                wave_conv.amplitudeZ=zscore(wave_conv.amplitude')';
                
                
                
                % ITPC
                epochWavePhase=epoch(wave_conv.phase,latencyDIN,[fromTime*EEG.srate toTime*EEG.srate]);
                epochWavePhase=rmbase(epochWavePhase,[],1:EEG.srate*0.1);
                ITPC(iCh,:,:)=abs(mean(exp(1i*epochWavePhase),3));
                
                % Amplitude
                epochWaveAmp=epoch(wave_conv.amplitude,latencyDIN,[fromTime*EEG.srate toTime*EEG.srate]);
                epochWaveAmp=rmbase(epochWaveAmp,[],1:EEG.srate*0.1);
                Amp(iCh,:,:)=mean(epochWaveAmp,3);
                
                % Z-scored Amplitude
                epochWaveAmpZ=epoch(wave_conv.amplitudeZ,latencyDIN,[fromTime*EEG.srate toTime*EEG.srate]);
                epochWaveAmpZ=rmbase(epochWaveAmpZ,[],1:EEG.srate*0.1);
                AmpZ(iCh,:,:)=mean(epochWaveAmpZ,3);
                
                
            end
            
            %             allWaveAmp{iSub,iSession}=waveAmp;
            allITPC{iSub,iSession}=ITPC;
            allAmp{iSub,iSession}=Amp;
            allAmpZ{iSub,iSession}=AmpZ;
            
            clear ITPC Amp AmpZ EEG latency type ;
            
            
        end
        
    end
end



%% Individual ITPC on topographical distribution

timeRange=[26:150];
freqRange=[40];
baseSes=1;

for iArm = 1 : nArm
    targetArm = Arms{iArm};
    topofig;
    for iSub = 1 : length(targetArm)
        
        for iSes = 1 : 4
            
            topoITPC = mean(mean(allITPC{Arms{iArm}(iSub),iSes}(:,freqRange,timeRange),3),2)-mean(mean(allITPC{Arms{iArm}(iSub),baseSes}(:,freqRange,timeRange),3),2);
            subplot(length(Arms{iArm}),4,(4*iSub-4)+iSes);
            topoplot(topoITPC,locs(ch),'numcontour',0,'intrad',0.53);
            set(gca,'YDir','normal');
            colorbar;
            %         caxis([-0.1 0.1]);
            title([ArmsName(iArm,:) ' : ' 'S' num2str(iSes)]);
            xlabel('time (ms)');
            ylabel('frequency (Hz)');
            
        end
        
    end
    
    
    clear topoITPC;
    
end




