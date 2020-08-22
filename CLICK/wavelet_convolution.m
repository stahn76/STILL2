function wave_conv=wavelet_convolution(data,params)
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



num_channels = size(data,1);
num_samples = size(data,2);

% Create wavelets if not specified in params
if ~isfield(params,'wavelet')
    if ~isfield(params,'Fs')
        Fs=800; % Set 800 Hz to be default sampling frequency
    else 
        Fs=params.Fs;
    end
    if ~isfield(params,'frange')
        frange=5:0.5:15; % Set default frequency range
    else 
        frange=params.frange;
    end
       
    fprintf('No wavelet specified;\n Initializing wavelets with Fs = %d Hz for frequencies %d to %d Hz at %1.1f Hz resolution\n',Fs,frange(1),frange(end),(frange(end)-frange(1))/(length(frange)-1))
    
    wavelet=cell(1,length(frange)); % Initialize cell to store wavelets
    
    for iter =1:length(frange) % Create wavelets
        wavelet{iter} = CreateWavelet_SA(frange(iter),Fs,3.5);         
    end
else
    wavelet=params.wavelet;
    frange=params.frange;
end
% Initialize zero matrices to store convolved signal
wave_conv_amp=zeros(length(frange),num_channels,num_samples);
wave_conv_phase=zeros(length(frange),num_channels,num_samples);

for iter=1:length(frange)
    conv_sig = convn(data,wavelet{iter},'same');% Be careful about convn (    
    wave_conv_amp(iter,:,:) = abs(conv_sig);% Convolved signal amplitude
    wave_conv_phase(iter,:,:) = angle(conv_sig);% Convolved signal phase
    if mod(frange(iter),10) == 0
        pct=iter*100/length(frange);
        fprintf('%3d%% ... ',round(pct));
    end
end
fprintf('\n');

wave_conv.amplitude=wave_conv_amp;
wave_conv.phase=wave_conv_phase;