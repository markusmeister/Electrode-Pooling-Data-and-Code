function [ scale_ratio ] = matchAmpscaler( waveformCells, target_p2p_amp )
%Determines the scaling ratio to scale the waveforms bank all with max
%negative peak of 400 to the median target p2p
%   Input: WaveformCells: a 1 by 12 cell that contains the waveform
%   templates all with same negative peak of 200
%   Input: target_p2p_amp: The target median amp of peak2peak amplitude

tmp_lis = [];

%get all p2p amp of the templates
for i = 1:size(waveformCells,2)
tmp_lis = [tmp_lis, max(waveformCells{i})-min(waveformCells{i})];
end

%determine the scaling ratio by the 
scale_ratio = target_p2p_amp / median(tmp_lis);
end

