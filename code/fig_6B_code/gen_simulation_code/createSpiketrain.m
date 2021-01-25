function [ ground_truth_cell, legit_spikeMat ] = createSpiketrain(fr, simLength,numTrials)
%Function that returns a binary matrix of numTrials by simulation length*sampling
% for creating simulation volatage trace in main script.
% ground truth cell holds the spike time for later scoring
%   Detailed explanation goes here
% fr: unit as spikes per second
% simLength: Length of the simulation in Seconds
% numTrials: Number of Max pool in the simulation

% Output:
% Ground truth Cell: Ground truth spike time places in a cell with the
% order of simulation
% Legit_spikeMat: a binary matrix numTrials by simulation length*sampling
% rate that places spike timings for convolving with spikeWaveform in the
% main script

% fix rnd
%seednum = 101;
%rng(seednum,'twister');
% Use: 2ms abs ITI
sampling_rate = 30000; % 30000 bins per S
abs_refractory_period = 2 * 10^-3; % 2 ms (in S)
refractory_bin_num = sampling_rate * abs_refractory_period; % time bin number of 2ms in sampling bins

[spikeMat, tVec] = poissonSpikeGen(fr, simLength, numTrials);

% find the time bins 
spike_bins = find(spikeMat'); % transpose to t by cell
diff_spike_bins = diff(spike_bins);
legit_logical_idxs = (diff_spike_bins > refractory_bin_num); % the logical positions where ITI are > 50 ms
legit_logical_idxs = [1;legit_logical_idxs]; % supplement the first spike back
legit_logical_idxs =logical(legit_logical_idxs);
%% create the final outputs

% logical indexing the spikes that are > refractory period
legit_spike_bins = spike_bins(legit_logical_idxs);
legit_spikeMat = zeros(size(spikeMat))';
legit_spikeMat(legit_spike_bins) = 1;
legit_spikeMat = legit_spikeMat';

ground_truth_cell = cell(1,numTrials);
for i = 1:numTrials
    % shift by 31 bins as the peak in waveform is at the 31th bin
    ground_truth_cell{1,i} = find(legit_spikeMat(i,:))+31; 
end

% function that generates the poisson spikes
function [spikeMat, tVec] = poissonSpikeGen(fr, tSim, nTrials)
% tSim : Unit in S
% fr: spikes per S
dt = 1/ 30000; % s
nBins = floor(tSim/dt);
spikeMat = rand(nTrials, nBins) < fr*dt;
tVec = 0:dt:tSim-dt;

end

end

