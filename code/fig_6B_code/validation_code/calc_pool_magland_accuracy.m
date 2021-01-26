function [max_units_num,unit_magland_accuracies ] = calc_pool_magland_accuracy( conf_mat )
%function that calculates the purity, defined as 
% (matches - wrongly merged spike - wrongly splitted - FP - FN) / ground
% truth
%   Detailed explanation goes here
% input: the confusion matrix from the magland script output
% output: 
% max_unit_num: max possible recovery unit count independent of quality
% unit_magland_accuracies: a vector of purities (float between 0 to 1)

max_units_num = min(size(conf_mat)) -1;

unit_magland_accuracies=zeros(1,size(conf_mat,1)-1);% place holder

% forloop for calculating each recovered unit's purity
for unit_i = 1:max_units_num
    
    %ground_truth_spike_counts = sum(conf_mat(unit_i,:)); %row sum
    recovered_spike_counts = conf_mat(unit_i,unit_i); % diagonal count that match truth
    %unit_accuracy = (recovered_spike_counts*3 - sum(conf_mat(unit_i,:)) - sum(conf_mat(:, unit_i)) ) / ground_truth_spike_counts;
    unit_accuracy = (recovered_spike_counts) / ...
    ( sum(conf_mat(unit_i,:)) + sum(conf_mat(:, unit_i)) - recovered_spike_counts ); % count everything once
    unit_magland_accuracies(unit_i) = unit_accuracy;
end

% add-in the unrecovered units as purity = 0 

end

