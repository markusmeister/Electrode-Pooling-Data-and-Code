% -------------------------------------------------------------------------
% Code to generate Data for plotting Figure 6 in:
% 
% Kyu Hyun Lee, Yu-Li Ni, Jennifer Colonell, Bill Karsh, Jan Putzeys,
% Marius Pachitariu, Timothy D. Harris, and Markus Meister (2021)
% Electrode pooling: boosting the yield of extracellular recordings with
% switchable silicon probes.

% This Script to count the Recovered Units with an accuracy score > 0.8 from sorted units' 
% spiketimes with respect to ground truth units' spiketimes.
% Confusion matrix was calculated using the matching algo. from 
% Barnett et al (2016), 
% https://github.com/ahbarnett/validspike
% -------------------------------------------------------------------------



% change to your root folder of the simulation
simroot = 'D:\Repo\Electrode-Pooling-Data-and-Code';
cd(simroot); %

% Add helper functions 
addpath(fullfile(simroot,'code','fig_6B_code','validation_code')) 
%load meta file: Information of parameter, groundtruth, maxsort success
load(fullfile(simroot,'data','simdata','simMetafile.mat'),'simMetafile');

% run the validation steps for each Param
Yield_M = zeros(12,15); % Place holder
% For each parameter, count the recovered units from pooling 1:12 tetrodes
for i = 1:size(simMetafile,1)
    [POOL_SCORE] = validation_fn( simroot, simMetafile{i,1}, simMetafile{i,2}, simMetafile{i,3} );
    Yield_M(1:length(POOL_SCORE.recovery_counts),i) = POOL_SCORE.recovery_counts';
end    
% save to the organized shape of each criteria
% In the format of header followed by # of units recovered above criteria
header = simMetafile(:,1);
fid = fopen( fullfile(simroot,'data','simdata','Pooling_data_reparse.csv'), 'w' );
fprintf( fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', header{:,1});
fprintf(fid,'%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d\n',Yield_M');
fclose( fid );

