% -------------------------------------------------------------------------
% Code to generate Data for plotting Figure 9 in:
% 
% Kyu Hyun Lee, Yu-Li Ni, Jennifer Colonell, Bill Karsh, Jan Putzeys,
% Marius Pachitariu, Timothy D. Harris, and Markus Meister (2021)
% Electrode pooling: boosting the yield of extracellular recordings with
% switchable silicon probes.

% This Script visualize accuracy scores from one example simulation (Fig9)
% 
% Accuracy score was defined using criteria from 
% Barnett et al (2016), 
% https://github.com/ahbarnett/validspike
% plotSpread function from
% Jonas (2021). plot spread points (beeswarm plot) 
% (https://www.mathworks.com/matlabcentral/fileexchange/37105-plot-spread-points-beeswarm-plot)
%  MATLAB Central File Exchange. Retrieved April 27, 2021.
% -------------------------------------------------------------------------

% change to your root folder of the simulation
simroot = 'D:\Repo\Electrode-Pooling-Data-and-Code';
cd(simroot); %

% Add helper functions 
addpath(fullfile(simroot,'code','plotSpread')) 

%load meta file: Information of accuracy scores
load(fullfile(simroot,'data','simdata','gridsearch',...
    '380_uV_10hz_noise_5_9_2_exp1','validation_','POOL_SCORE.mat'),...
    'POOL_SCORE');

% Bee Swarm plot (use plotSpread function from matlab central)

f = figure('Position', [0   0   1200   600]);
A = subplot('Position', [0.1, 0.12, 0.85, 0.8]);
set(A, 'FontSize', 14)
A.TickDir = 'out';

ax = gca;
ax.FontSize=14;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;
%ax.linewidth=1.5;
ax.TickDir='out';

hold on

lin=plot(0:13, ones(1,14)*POOL_SCORE.good_unit_purity_thresh,'-.', 'color', 'k');
xlim([0,13])
ylim([-0.01,1.05])

plotSpread(POOL_SCORE.unit_magland_accuracies,'distributionColors','k',...
    'distributionMarkers','o')
xlabel('Number of Pools' ,'FontSize',18);
ylabel('Accuracy Score','FontSize',16);
legend([lin],{'threshold'},'Location','southwest')


fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print(fig,fullfile(simroot,'figs',... 
'Accuracy_Swarm.pdf'),'-dpdf')