function [] = helper_save_confusionmat(savepath, pool_num, confusion_matrix)
% function that saves both the visualization of the confusion matrix and
% the numerical value of the matrix

% input: savepath (str) of the save path directory
% input: pool_num: (int), int of number of pools used in this
% save numerical val
    

    % save the figure
    % get time diff of pk trough
    fig = figure;
    ax =axes;
    imagesc(confusion_matrix);
    colorbar;
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 4])
    fig_name = (sprintf('pooling confusion matrix%d', pool_num));
    title(fig_name)
    xlabel('sorted unit num')
    ylabel('ground truth unit num')
    
    saveas(fig, fullfile(savepath, sprintf('pool%d.png',pool_num) ))
    cla(ax)
    close(fig)

end
