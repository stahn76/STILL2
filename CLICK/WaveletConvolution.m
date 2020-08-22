function wave_conv=WaveletConvolution(data,params)
%-------------------------------------------------------------------------%
% Function to convolve a short segment of data with Gabor wavelet at
% specified sampling rate, frequency range and resolution
%
% Input: data   - ch (1-channel) x ns x ntrial
%                   ns - number of samples
%                   ntrial - number of trials
%        params - wavelet - cell array of wavelets generate by
%                           CreateWavelet.m
%               - Fs      - sampling frequency in Hz 
%               - frange  - 1xn vector for frequency range for spectrogram 
%              
% Output: wave_conv - struct with amplitude and phase
%-------------------------------------------------------------------------
% - SA 8/6/2019 changed the function from convn to conv



%%
% Create wavelets if not specified in params
if ~isfield(params,'wavelet')
    if ~isfield(params,'Fs')
        Fs=500; % Set 500 Hz to be default sampling frequency
    else 
        Fs=params.Fs;
    end
    if ~isfield(params,'frange')
        frange=1:0.5:50; % Set default frequency range
    else 
        frange=params.frange;
    end
       
    fprintf(['No wavelet specified; Initializing wavelets with' ...
        'Fs=%d Hz, %d to %d Hz at %1.1f Hz bin \n'],...
        Fs,frange(1),frange(end),(frange(end)-frange(1))/(length(frange)-1));
    
    wavelet=cell(1,length(frange)); % Initialize cell to store wavelets
    
    for iter =1:length(frange) % Create wavelets
        wavelet_temp=CreateWaveletShort(frange(iter),Fs);
        wavelet{iter}={wavelet_temp{1},wavelet_temp{2}};
    end
else
    wavelet=params.wavelet;
    frange=params.frange;
end
%%
% get number of sample and trial
nSample=size(data,2);
nTrial = size(data,3);

% Initialize zero matrices to store convolved signal
nTrial = size(data,3);

if nTrial==1
    wave_conv_amp=zeros(num_channels,num_samples,length(frange));
    wave_conv_phase=zeros(num_channels,num_samples,length(frange));
    
    for iter=1:length(frange)
        re=convn(data,wavelet{iter}{1},'same');
        im=convn(data,wavelet{iter}{2},'same');
        conv_sig=re+1i*im;
        wave_conv_amp(:,:,iter)=abs(conv_sig);% Convolved signal amplitude
        wave_conv_phase(:,:,iter)=angle(conv_sig);%atan2(re,im);% Convolved signal phase
        
        if mod(frange(iter),10)==0
            pct=iter*100/length(frange);
            fprintf('%3d%% ... ',round(pct));
        end
    end
    
else
    wave_conv_amp=zeros(num_channels,num_samples,length(frange),nTrial);
    wave_conv_phase=zeros(num_channels,num_samples,length(frange),nTrial);
    
    for iTrial = 1 : nTrial
        for iter=1:length(frange)
            re=convn(data(:,:,iTrial),wavelet{iter}{1},'same');
            im=convn(data(:,:,iTrial),wavelet{iter}{2},'same');
            conv_sig=re+1i*im;
            wave_conv_amp(:,:,iter,iTrial)=abs(conv_sig);% Convolved signal amplitude
            wave_conv_phase(:,:,iter,iTrial)=angle(conv_sig);%atan2(re,im);% Convolved signal phase
            
%             if mod(frange(iter),10)==0
%                 pct=iter*100/length(frange);
%                 fprintf('%3d%% ... ',round(pct));
%             end
        end
        fprintf('Convolving ..... Trial # %d.... \n',iTrial);       
    end
    
end
    


%% outputs
wave_conv.amplitude=wave_conv_amp;
wave_conv.phase=wave_conv_phase;

% wave_conv=wave_conv_amp;

