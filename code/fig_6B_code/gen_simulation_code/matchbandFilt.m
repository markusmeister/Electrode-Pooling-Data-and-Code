function [ filtNoiseMat, rms_noise ] = matchbandFilt( noiseMat )
%A function that filtered the noiseMat (chno by data-len) with a third
%order butterworth filter [300,5000]
%   Input: generated noise with size chno by data-len
%   Output: filtNoiseMat: the Band Passed Input noise, in orriginal dim
%   Output: rms_noise of filtNoiseMat

%Setting the butterworth filter
samplingrate=30e3;
filt_band = [300,5000];
filt_order = 3;
filt_passtype = 'bandpass';
[b,a] = butter(filt_order,filt_band/(samplingrate/2),filt_passtype);

%filtfilt the Noise matrix
filtNoiseMat = filtfilt(b, a, noiseMat')';

%record the RMS level of the noise
rms_noise = std(reshape(filtNoiseMat,1,[]));

end

