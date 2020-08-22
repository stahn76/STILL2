function wave_conv=wavelet_convolution_STILL2_CLICK(data,params)
%-------------------------------------------------------------------------%
% Function to convolve a short segment of data with Gabor wavelet at
% specified sampling rate, frequency range and resolution
%
% Input: data   - nch x ns matrix
%                   nch - number of channels
%                   ns - number of samples
%        params - wavelet - cell array of wavelets generate by
%                           CreateWavelet.m
%               - Fs      - sampling frequency in Hz 
%               - frange  - 1xn vector for frequency range for spectrogram 
%              
% Output: wave_conv - struct with amplitude and phase
%-------------------------------------------------------------------------%
% - SA 05/30/2014


num_channels=size(data,1); %Or trials in SOME cases
num_samples=size(data,2);

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
       
    fprintf('No wavelet specified; Initializing wavelets with Fs = %d Hz for frequencies %d to %d Hz at %1.1f Hz resolution',Fs,frange(1),frange(end),(frange(end)-frange(1))/(length(frange)-1))
    
    wavelet=cell(1,length(frange)); % Initialize cell to store wavelets
    
    for iter =1:length(frange) % Create wavelets
        wavelet_temp=CreateWaveletShort(frange(iter),Fs);
        wavelet{iter}={wavelet_temp{1},wavelet_temp{2}};
    end
else
    wavelet=params.wavelet;
    frange=params.frange;
end
% Initialize zero matrices to store convolved signal
wave_conv_amp=zeros(length(frange),num_channels,num_samples);
wave_conv_phase=zeros(length(frange),num_channels,num_samples);

for iter=1:length(frange)
    re=convn(data,wavelet{iter}{1},'same');
    im=convn(data,wavelet{iter}{2},'same');
    conv_sig=re+1i*im;
    wave_conv_amp(iter,:,:)=abs(conv_sig);% Convolved signal amplitude
    wave_conv_phase(iter,:,:)=angle(conv_sig);%atan2(re,im);% Convolved signal phase
    if mod(frange(iter),10)==0
        pct=iter*100/length(frange);
        fprintf('%3d%% ... ',round(pct));
    end
end
fprintf('\n');

%% outputs
wave_conv.amplitude=wave_conv_amp;
wave_conv.phase=wave_conv_phase;

% wave_conv=wave_conv_amp;

