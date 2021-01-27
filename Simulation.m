%main for re-extracting accuracy Yield Pool data needed for plotting fig6
% (The recovery counts that goes into /data/Pooling_data.csv)

% change to your root folder of the simulation
simroot = 'D:\Repo\Electrode-Pooling-Data-and-Code';
cd(simroot); %

% helper functions 
addpath(fullfile(simroot,'code','fig_6B_code','validation_code')) 
%load meta file: Information of parameter, groundtruth, maxsort success
load(fullfile(simroot,'data','simdata','simMetafile.mat'),'simMetafile');

% run the validation steps for each simpool of each grid param
%simroot, simfolderName, grdTruthMat_Name, maxpoolednum
Yield_M = zeros(12,15);

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

