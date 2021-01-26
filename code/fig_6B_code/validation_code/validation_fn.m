function [POOL_SCORE] = validation_fn( simroot, simfolderName, grdTruthMat_Name, Max_Pool )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Update the Folder name that contains the kilosorted folders
simDatafolder = fullfile(simroot,'data' ,'simdata', 'gridsearch', simfolderName);

% update ground truth cell for matching sim vs truth
load(fullfile(simDatafolder, grdTruthMat_Name),'ground_truth_cell');

% Update Fig title and successful sorted pool num
fig_title = simfolderName; 

% creates save dir called validation
if exist(fullfile(simDatafolder, 'validation_'),'dir')~=7
    [status,msg] = mkdir(fullfile(simDatafolder, 'validation_'));
    assert(status,'Could not create directory ''%s'': %s',fullfile(simDatafolder, 'validation_'),msg);
    
end
confusion_Mat_savedir = fullfile(simDatafolder, 'validation_');

% for loop every thing
for N_pool = 1:Max_Pool % loop thru all the confusion matrixes
     
    sorted_path = fullfile(simDatafolder , sprintf('simpool%s',num2str(N_pool) ) , 'spiketrain_good.mat');
    % load sorted spike time
    load(sorted_path, 'sp')
    
    
    % border case of only one unit entry, sp != cell array but a vector, have to cast
    % it to cell so that the format is correct for the later function call
    %if N_pool == 1
    if ~isa(sp,'cell')    
        sp_test = cell(1,1);
        sp_test{1} = sp;
        sp = sp_test;
    end
    % formating savepath
    confusion_Mat_filesavepath = fullfile(confusion_Mat_savedir, sprintf('confusion_matrix%s.mat', num2str(N_pool)) );
    perm_filesavepath = fullfile(confusion_Mat_savedir, sprintf('L2perm%s.mat', num2str(N_pool)) );
    % Run thru the Barnett/Magland et al algo.
    [truth_T_list, truth_L_list, sim_T_list, sim_L_list]=helper_confusion_matrix2(N_pool, sp, ground_truth_cell);
    % times_label function only accept row >= column. Format the inputs 
    %if max(sim_L_list)> max(truth_L_list)
    if max(sim_L_list)>= max(truth_L_list)    
        [confusion_matrix, T1, T2, T2w, permL2]=times_labels_confusion_matrix(truth_T_list,truth_L_list,sim_T_list,sim_L_list);
        % default takes identify +- 3 /30,000 sec as the same spike
        % save the result
        save(confusion_Mat_filesavepath, 'confusion_matrix')
        save(perm_filesavepath, 'permL2')
    else
        
        [confusion_matrix, T1, T2, T2w, permL2]=times_labels_confusion_matrix(sim_T_list,sim_L_list,truth_T_list,truth_L_list);
        % default takes identify +- 3 /30,000 sec as the same spike
        % transpose so that ground truth axis is still in the rows
        confusion_matrixT = confusion_matrix' ;
        % save the result
        save(confusion_Mat_filesavepath, 'confusion_matrixT')
        save(perm_filesavepath, 'permL2')
        
    end
    
    clear truth_T_list truth_L_list sim_T_list sim_L_list

end

%% scoring and visualization
% construct the scoring struct
POOL_SCORE = struct;
POOL_SCORE.pools_idx = 1:Max_Pool;
POOL_SCORE.max_recovery_counts = zeros(1,Max_Pool);
POOL_SCORE.recovery_counts = zeros(1,Max_Pool);
%POOL_SCORE.recovery_purities = cell(1,Max_Pool); % place holder for storing the purities of a a pool test
POOL_SCORE.good_unit_purity_thresh = 0.8; % the units have to pass threshold to make it into a recovered unit

for i = 1:Max_Pool
    
    test_mat = load(fullfile(confusion_Mat_savedir , sprintf('confusion_matrix%d.mat', i))); 
    try 
        test_mat = test_mat.confusion_matrixT;
    catch 
        % when there are matrix with sim units > truth units
        %warning( sprintf('confusion_matrix%d has more sim units then truth units', i) );
        test_mat = test_mat.confusion_matrix;
    end
    % save the confusion matrix visualization
    helper_save_confusionmat(confusion_Mat_savedir, i, test_mat)
    
    % accuracy score
    [max_unit_num, unit_magland_accuracies ] = calc_pool_magland_accuracy( test_mat );
    POOL_SCORE.max_recovery_counts(i) = max_unit_num;
    POOL_SCORE.unit_magland_accuracies{1,i} = unit_magland_accuracies;
    % return yield
    POOL_SCORE.recovery_counts(1,i)=sum(POOL_SCORE.unit_magland_accuracies{1,i}(1:POOL_SCORE.max_recovery_counts(i))> POOL_SCORE.good_unit_purity_thresh);

end

%% counting and plotting using magland accuracy
% yield plot

fig= figure;
pool_t = 1:Max_Pool;
% extract the qualified units
qualified_units_counts = zeros(1,Max_Pool);

for i = 1:Max_Pool
    
    qualified_units_counts(1,i)=sum(POOL_SCORE.unit_magland_accuracies{1,i}(1:POOL_SCORE.max_recovery_counts(i))> POOL_SCORE.good_unit_purity_thresh);

end
hold on
scatter(pool_t, qualified_units_counts, 80, 'k', 'LineWidth',3)
plot(pool_t, pool_t, '-.', 'linewidth', 4)
legend('Actual Recovery', 'Best Case Recovery' , 'location', 'northwest')
title(fig_title)
xlabel('# of Pools')
ylabel('# of good unit Yield')
xlim([0 12]);ylim([0 12])

x0=10;
y0=10;
width=800;
height=600;
set(gcf,'position',[x0,y0,width,height])
saveas(fig, fullfile(confusion_Mat_savedir, sprintf('UnitsYield_%s.png',fig_title) ))
%cla(ax)
close(fig)
%% plot distribution of accuracy
fig= figure;
%pool_t = 1:Max_Pool;
% extract the qualified units
%qualified_units_counts = zeros(1,Max_Pool);

for i = 1:Max_Pool
    hold on
    scat=scatter(ones(1,POOL_SCORE.max_recovery_counts(i))*i ,...
        POOL_SCORE.unit_magland_accuracies{1,i}(1:POOL_SCORE.max_recovery_counts(i)),80, 'LineWidth',3);

end

lin=plot(0:13, ones(1,14)*POOL_SCORE.good_unit_purity_thresh,'-.', 'color', 'r');
xlim([0,13])
%ylim([-0.5,1.1])
ylabel('Accuracy score')
xlabel('# of Pools')
set(gcf,'position',[x0,y0,width,height])
set(gca, 'TickDir', 'out')
legend([scat lin],{'sorted unit accuracy score', '0.8 threshold'},'Location','southwest')

title('Scoring: Accuracy Distribution')
saveas(fig, fullfile(confusion_Mat_savedir, sprintf('Accuracy Distribution_%s.png',fig_title) ))

close(fig)
% Save the Scoring Struct
save(fullfile(confusion_Mat_savedir, 'POOL_SCORE.mat') , 'POOL_SCORE')


%%
function [truth_T_list, truth_L_list, sim_T_list, sim_L_list]=helper_confusion_matrix2(N_pool, sp_cell, ground_truth_cell)

% Input N pool: int of how many banks one is using
 % sp cell: a one by n cell array containing the kilosorted units spike
 % times
% ground_truth_cell: cell array containing the original ground
% truth of sim units

% Output: truth_T_list : an 1 d array of spike times of all units
% truth L list : an 1 d array of corresponding unit number (1:n) that
% corresponds to the spike times
% sim_T_list: an 1 d array of spike times of all units of the simulation
% units 
% sim_L_list: an 1d array of corresponding unit number (1: available sorted unit)
% that corresponds to the spiking array

% construct truth_T_list, truth_L_list
truth_T_list = [];
truth_L_list = [];
for tru_i = 1:N_pool
    truth_T_list = [truth_T_list, ground_truth_cell{tru_i} ]; % load the ground truth spiketimes
    truth_L_list = [truth_L_list, ones(1, length(ground_truth_cell{tru_i}))*tru_i ]; % create the corresponding idx
end

% construct the sim_T_list sim_L_list
sim_T_list =[];
sim_L_list = [];
for sim_i = 1:size(sp_cell,2)
    
    %sim_T_list = [sim_T_list; sp_cell{sim_i}]; % access and append the spike times of the sorted units
    % using column convention
    sim_T_list = [sim_T_list; reshape(sp_cell{sim_i},[],1)]; % access and append the spike times of the sorted units
    sim_L_list = [sim_L_list, ones(1, length(sp_cell{sim_i})) *sim_i] ; % create labels
end
    sim_T_list = sim_T_list';

end


end

