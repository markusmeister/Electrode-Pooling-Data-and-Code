% -------------------------------------------------------------------------
% Code to generate Figure 5 in:
% 
% Kyu Hyun Lee, Yu-Li Ni, Jennifer Colonell, Bill Karsh, Jan Putzeys,
% Marius Pachitariu, Timothy D. Harris, and Markus Meister (2021)
% Electrode pooling: boosting the yield of extracellular recordings with
% switchable silicon probes
% -------------------------------------------------------------------------
clear all

figure('position',[100, 100, 800, 650]);

% -------------------------------------------------------------------------
% A) Sample neurons
% -------------------------------------------------------------------------

% Panel label
ax_AA = subplot('position', [0.02, 0.97, 0.03, 0.03]);
text(ax_AA, 0, 0.5, 'A', 'FontSize', 13,  'FontName','Arial');
axis(ax_AA,'off');

% Positions of the subplots (waveforms and ISI)
ax1 = subplot('Position',[0.05 0.65 0.2 0.4]);
ax2 = subplot('Position',[0.28 0.65 0.2 0.4]);
ax3 = subplot('Position',[0.52 0.65 0.2 0.4]);
ax4 = subplot('Position',[0.76 0.65 0.2 0.4]);
ax12 = subplot('Position',[0.06 0.55 0.2 0.08]);
ax22 = subplot('Position',[0.29 0.55 0.2 0.08]);
ax32 = subplot('Position',[0.53 0.55 0.2 0.08]);
ax42 = subplot('Position',[0.77 0.55 0.2 0.08]);

% Indices of sample cells (one from bank0 and one from bank1 that were
% detected by the same set of electrodes)
bank0_cellidx = 315;
bank1_cellidx = 192;

% Plot waveforms and ISI
plot_wfs([ax1, ax2, ax3, ax4, ax12, ax22, ax32, ax42], bank0_cellidx, bank1_cellidx);

% Set YLim for the ISI plots
ylim(ax12,[0,260]);
xlim(ax12,[0,0.05]*1000);
ylabel(ax12, 'Count');
ylim(ax22,[0,260]);
xlim(ax22,[0,0.05]*1000);
ylim(ax32,[0,30]);
xlim(ax32,[0,0.05]*1000);
ylim(ax42,[0,30]);
xlim(ax42,[0,0.05]*1000);

% Axis labels
ax_l = subplot('Position',[0.48, 0.49, 0.04, 0.04]);
text(ax_l, 0.5,0.5,'Inter-spike interval (ms)','HorizontalAlignment','center',...
    'FontName','Arial','FontSize',8);
axis(ax_l, 'off');

% -------------------------------------------------------------------------
% B) Similarity matrix
% -------------------------------------------------------------------------

% Panel label
ax_A = subplot('position', [0.02, 0.46, 0.05, 0.05]);
text(ax_A, 0, 0.5, 'B', 'FontSize', 13, 'FontName','Arial');
axis(ax_A,'off');

% Load waveforms and cell IDs of manually sorted split- and pooled-mode cells
load(fullfile('.','data','manual.mat'),'cellid_s', 'cellid_bank01', 'wfs_s', 'wf_bank01');

% Load probe geometry
load(fullfile('.','data','chmap_neuropix3B.mat'),'ycoords');

% Organize the split-mode cells waveforms and IDs by depth
max_amp_loc_s = zeros(1,size(wfs_s,2));
for i=1:size(wfs_s,2)
    wf = wfs_s{1,i};
    amp = max(wf,[],2)-min(wf,[],2);
    [~,maxind] = max(amp);
    max_amp_loc_s(i) = ycoords(maxind);
end
clear i wf amp maxind
[~, inds] = sort(max_amp_loc_s,'descend');
cellid_s = cellid_s(inds);
wfs_s = wfs_s(inds);
clear inds

% Organize the pooled-mode cells waveforms and IDs by depth
max_amp_loc_p = zeros(1,size(wf_bank01,2));
for i=1:size(wf_bank01,2)
    wf = wf_bank01{1,i};
    amp = max(wf,[],2)-min(wf,[],2);
    [~,maxind] = max(amp);
    max_amp_loc_p(i) = ycoords(maxind);
end
clear i wf amp maxind
[~,inds]=sort(max_amp_loc_p,'descend');
cellid_bank01 = cellid_bank01(inds);
wf_bank01 = wf_bank01(inds);
clear inds

% Threshold for cosine similarity score; if below this then consider not
% good match
cosine_sim_score_th = 0.9;

% Compute similarity matrix (based on the 20 channels around the peak
% channel)
nchans = 20;
s = compute_sim2(cellid_s, cellid_bank01, wfs_s, wf_bank01, nchans);

% Red-blue colormap (nonlinear to emphasize the middle range)
colsize=71;
col = ones(colsize,3); 
col(1:ceil(colsize/2),2) = (linspace(0,1,ceil(colsize/2))).^(0.3);
col(1:ceil(colsize/2),3) = (linspace(0,1,ceil(colsize/2))).^(0.3);
col(ceil(colsize/2):colsize,1) = (linspace(1,0,ceil(colsize/2))).^(0.3);
col(ceil(colsize/2):colsize,2) = (linspace(1,0,ceil(colsize/2))).^(0.3);
col = flipud(col);

% Plot similarity matrix
ax1 = subplot('position', [0.07, 0.1, 0.35,0.35]);
hold(ax1, 'on');
imagesc(ax1, 'Cdata', s, [-1,1]);
ax1.FontSize = 8;
ax1.FontName = 'Arial';
xlabel(ax1, 'Pooled-mode cells', 'FontSize', 9, 'FontName', 'Arial'); 
ylabel(ax1, 'Split-mode cells', 'FontSize', 9, 'FontName', 'Arial'); 
axis(ax1, 'tight'); 
colormap(ax1, col);
colorbar(ax1,'FontSize',8);
ax1.TickDir = 'out';
ax1.YDir = 'reverse';

[rowInd, colInd] = find(s>cosine_sim_score_th);
plot(ax1, colInd, rowInd, 'k.','MarkerSize',1);
hold(ax1, 'off');

% Plot a distribution of similarity scores
ax2 = subplot('position', [0.48, 0.1, 0.17, 0.35]);
hold(ax2, 'on');
h = histogram(ax2, s(:),-1:0.05:1,'DisplayStyle','stairs','LineWidth',1.5,...
    'EdgeColor',[0.3 0.3 0.3]);
plot(ax2, [cosine_sim_score_th,cosine_sim_score_th],[0,max(h.Values)*1.05],'k--');
hold(ax2, 'off');
xlim(ax2, [-1,1]);
ylim(ax2, [0,max(h.Values)*1.05]);
xlabel(ax2, 'Cosine similarity','FontSize', 9, 'FontName', 'Arial'); 
ylabel(ax2, 'Count','FontSize', 9, 'FontName', 'Arial'); 
ax2.Box = 'off';
ax2.TickDir='out';
ax2.FontSize = 8;
ax2.FontName = 'Arial';
clear h

ax2_inset = axes('Position',[0.61 0.37 0.05 0.07]);
box(ax2_inset,'on');
h = histogram(ax2_inset, s(:),-1:0.05:1,'DisplayStyle','stairs','LineWidth',1.5,...
    'EdgeColor',[0.3 0.3 0.3]);
xlim(ax2_inset, [0.7,1]);
ylim(ax2_inset, [0,200]);
ax2_inset.TickDir='out';

clear h
% -------------------------------------------------------------------------
% C) Hot sorting
% -------------------------------------------------------------------------

% Panel label
ax_hotsort_label = subplot('position', [0.67, 0.46, 0.05, 0.05]);
text(ax_hotsort_label, 0, 0.5, 'C', 'FontSize', 13,  'FontName','Arial');
axis(ax_hotsort_label,'off');

% Axis where the plot will go
ax_hotsort = subplot('position',[0.73, 0.1, 0.22, 0.35]);

bin_edges = 0:25:600;
p_match_manual = match_amp('manual', bin_edges);
p_match_coldsort = match_amp('coldsort', bin_edges);
p_match_hotsort = match_amp('hotsort', bin_edges);

p_match_manual(isnan(p_match_manual))=0;
p_match_coldsort(isnan(p_match_coldsort))=0;
p_match_hotsort(isnan(p_match_hotsort))=0;

color_manual = [0,158,119]/255;
color_hotsort = [217,95,0]/255;
color_coldsort = [110,110,170]/255;

hold(ax_hotsort,'on');
histogram(ax_hotsort,'BinEdges',bin_edges,'BinCounts',p_match_manual,...
    'DisplayStyle','stairs','LineWidth',1.5,'EdgeColor',color_manual);
histogram(ax_hotsort,'BinEdges',bin_edges,'BinCounts',p_match_coldsort,...
    'DisplayStyle','stairs','LineWidth',1.5,'EdgeColor',color_coldsort); 
histogram(ax_hotsort,'BinEdges',bin_edges,'BinCounts',p_match_hotsort,...
    'DisplayStyle','stairs','LineWidth',1.5,'EdgeColor',color_hotsort); 
text(160, 0.3, 'Manual','FontSize', 8, 'FontName', 'Arial','Color',color_manual); 
text(160, 0.25, 'Cold sort','FontSize', 8, 'FontName', 'Arial','Color',color_coldsort); 
text(160, 0.2, 'Hot sort','FontSize', 8, 'FontName', 'Arial','Color',color_hotsort); 
plot(ax_hotsort,[0,270],[0.5,0.5],'k--');
hold(ax_hotsort,'off');
xlim(ax_hotsort,[0,270])
ax_hotsort.TickDir='out';
xlabel(ax_hotsort,'Spike amplitude (\muV)');
ylabel(ax_hotsort,'Prob(match)');

% -------------------------------------------------------------------------
% Save output
% -------------------------------------------------------------------------
% Uncomment to save figure as svg
% saveas(gcf,fullfile('.','figs','fig5-InVivo.svg'))
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];
print('-bestfit',fig,fullfile('.','figs','fig5-InVivo.pdf'),'-dpdf')

%%
% -------------------------------------------------------------------------
% Functions called
% -------------------------------------------------------------------------
function x = id_to_ind(ids, cellids)
x = zeros(1,length(ids));
for i=1:length(ids)
    x(i) = find(cellids==ids(i));
end
clear i
end

function amp = find_p2p_amp_at_peak_channel(wfs)
% -------------------------------------------------------------------------
% Given a cell array of waveforms (each of which is channels x samples)
% returns the peak-to-peak amplitude at the channel with highest amplitude
% -------------------------------------------------------------------------
amp = zeros(1,length(wfs));
for i = 1:length(wfs)    
    S = wfs{1,i};
    p2p = max(S,[],2)-min(S,[],2); % peak-to-peak amplitude
    amp(i) = max(p2p);
end
clear i

end

function p_match = match_amp(type, bin_edges)
if strcmp(type, 'manual')
% Load waveforms and cell IDs from manual sorting
    load(fullfile('.','data','manual.mat'),'cellid_s', 'cellid_bank01', 'wfs_s', 'wf_bank01');
elseif strcmp(type, 'coldsort')
    load(fullfile('.','data','coldsort.mat'),'cellid_s', 'cellid_bank01', 'wfs_s', 'wf_bank01');
elseif strcmp(type, 'hotsort')
    load(fullfile('.','data','hotsort.mat'),'cellid_s', 'cellid_bank01', 'wfs_s', 'wf_bank01');
else
    error('recording type undefined');
end

% Threshold for cosine similarity score; if below this then consider not
% good match
cosine_sim_score_th = 0.9;

% Compute similarity matrix (based on the 20 channels around the peak
% channel)
nchans = 20;
s = compute_sim2(cellid_s, cellid_bank01, wfs_s, wf_bank01, nchans);

% Get good matches for manual sort
m = gen_match(s, cellid_s, cellid_bank01);

% Get IDs of split-mode cells that found a match
split_match_ids = m(1, m(3,:)>cosine_sim_score_th);

amps_s = find_p2p_amp_at_peak_channel(wfs_s);

[total,~] = histcounts(amps_s, bin_edges);
[matched,~] = histcounts(amps_s(id_to_ind(split_match_ids, cellid_s)), bin_edges);
p_match = matched./total;
end

function plot_wfs(ax, cell_bank0, cell_bank1)
% -------------------------------------------------------------------------
% Plots the waveforms and ISI of specified cells
% 
% Parameters
% ----------
% ax : array of 8 Axis objects
%   [ax1, ax2, ax3, ax4, ax12, ax22, ax32, ax42]
% cell_bank0 : int
%   index of cell from bank0
% cell_bank1 : int
%   index of cell from bank1
% -------------------------------------------------------------------------

% Set color for the two banks
col_bank0 = [228,26,28]/255;
col_bank1 = [55,126,184]/255;

% Load probe geometry
load(fullfile('.','data','chmap_neuropix3B.mat'), 'xcoords', 'ycoords');

% Load waveforms and cell IDs from manual sorting
load(fullfile('.','data','manual.mat'),'cellid_s', 'cellid_bank01', 'wfs_s', 'wf_bank01');

% Compute cosine similarity between split and pooled cells based on the
% waveforms that exceed a threshold of 25 microV
th = 25;
s = compute_sim(cellid_s, cellid_bank01, wfs_s, wf_bank01, th);
m = gen_match(s, cellid_s, cellid_bank01);

cellidx_bank0_pooled_match = m(2,m(1,:)==cell_bank0);
cellidx_bank1_pooled_match = m(2,m(1,:)==cell_bank1+10000);

% Set waveform (choose the best matching pooled cell)
S0 = wfs_s{1,cellid_s==cell_bank0};
P0 = wf_bank01{1,cellid_bank01==cellidx_bank0_pooled_match};

S1 = wfs_s{1,cellid_s==(cell_bank1+10000)};
P1 = wf_bank01{1,cellid_bank01==cellidx_bank1_pooled_match};

p2p_amp = max(S0,[],2)-min(S0,[],2);
chanstoplot = find(p2p_amp>th)';

% Fit split cell waveform to pooled cell waveform
pc0 = get_pc(S0, P0);
pc1 = get_pc(S1, P1);

% Scale factor for plotting waveform
yfactor = 6.5;
xfactor = 4;

% Define time vector
t = (1:size(S0,2))/xfactor;

% Scale wfs
S0 = S0/yfactor;
P0 = P0/yfactor;
S1 = S1/yfactor;
P1 = P1/yfactor;

for i = ax(1:4)
    hold(i, 'on');
end

for i=chanstoplot
    % Waveform, split neuron 1
    plot(ax(1), t+xcoords(i), S0(i,:)+ycoords(i), 'color', col_bank0, 'LineWidth', 1);
    % Waveform, pooled neuron 1
    plot(ax(2), t+xcoords(i), P0(i,:)+ycoords(i), 'color', [0,0,0], 'LineWidth', 1);
    % Waveform, pooled neuron 2
    plot(ax(3), t+xcoords(i), P1(i,:)+ycoords(i), 'color', [0,0,0], 'LineWidth', 1);
    % Waveform, split neuron 2
    plot(ax(4), t+xcoords(i), S1(i,:)+ycoords(i), 'color', col_bank1, 'LineWidth', 1);
    % Pooling coefficient
    text(ax(2), t(end)+xcoords(i), P0(i,end)+ycoords(i), num2str(round(pc0(i),2)),...
        'HorizontalAlignment', 'left', 'FontName','Arial', 'FontSize', 8);
    text(ax(3), t(end)+xcoords(i), P1(i,end)+ycoords(i), num2str(round(pc1(i),2)),...
        'HorizontalAlignment', 'left', 'FontName','Arial', 'FontSize', 8);
end
clear i

% Vertical scale bar
plot(ax(4), [max(xcoords)+length(t)/xfactor+10,  max(xcoords)+length(t)/xfactor+10],...
    [min(ycoords(chanstoplot)),  min(ycoords(chanstoplot))+50/yfactor], 'k', 'LineWidth', 2);
text(ax(4), max(xcoords)+length(t)/xfactor+13,  min(ycoords(chanstoplot))+30/yfactor,...
    {'50'; '\muV'},'HorizontalAlignment','left', 'FontName', 'Arial', 'FontSize', 8);
% Horizontal scale bar
plot(ax(4), [max(xcoords)+length(t)/xfactor+10, max(xcoords)+0.5*length(t)/xfactor+10],...
    [min(ycoords(chanstoplot)),  min(ycoords(chanstoplot))],'k','LineWidth',2);
text(ax(4), max(xcoords)+0.75*length(t)/xfactor+10,  min(ycoords(chanstoplot))-15/yfactor,...
    '1 ms', 'HorizontalAlignment','center', 'FontName', 'Arial', 'FontSize', 8);

% Set X and YLims
for i=ax(1:4)
    xlim(i, [min(xcoords)-3, max(xcoords)+30]);
    ylim(i, [min(ycoords(chanstoplot))-10, max(ycoords(chanstoplot))+30]);
    axis(i, 'off');
    hold(i, 'off');
end

% Plot ISI
plot_isi(ax(5),'bank0',cell_bank0)
plot_isi(ax(6),'bank01',cellidx_bank0_pooled_match)
plot_isi(ax(7),'bank01',cellidx_bank1_pooled_match)
plot_isi(ax(8),'bank1',cell_bank1)

end

function plot_isi(ax, bank, cell)
% -------------------------------------------------------------------------
% Plots inter-spike interval distribution for a given cell
% 
% Parameters
% ----------
% ax : Axis object
%   where plot output will go
% bank : str
%   specifies which bank the cell is from (0, 1, or pooled)
% cell : int
%   cell ID
% -------------------------------------------------------------------------

% Sampling rate
samplingrate = 30e3;
% Load spike trains
if strcmp(bank, 'bank0')
    load(fullfile('.','data','manual_bank0_sp.mat'),'sp','cellid');
elseif strcmp(bank, 'bank1')
    load(fullfile('.','data','manual_bank1_sp.mat'),'sp','cellid');
else
    load(fullfile('.','data','manual_bank01_sp.mat'),'sp','cellid');
end
% Find spike train
spikes = double(sp{1,cellid==cell})/samplingrate;
% Plot ISI for 100 ms
bins = 0:100;
h = histogram(ax, diff(spikes)*1000, bins,'FaceColor',[0,0,0],'EdgeColor','none');
ax.TickDir = 'out';
ax.Box = 'off';
ax.FontSize = 7.5;
ax.FontName = 'Arial';
title(ax, ['Mean firing rate: ', num2str(round(length(spikes)/spikes(end),2)), ' sp/s'],...
    'FontSize',9,'FontName','Arial','FontWeight','normal', 'FontSize', 8);
ylim(ax, [0, max(h.Values)*1.1]);
end

function m = gen_match(s, cellid_s, cellid_bank01)
% -------------------------------------------------------------------------
% Matches split to pooled based on similarity score
% 
% Parameters
% ----------
% s : array, (a, b)
%   similarity matrix
% cellid_s : array, (1,a)
%   IDs of split cells
% cellid_bank01 : array, (1,b)
%   IDs of pooled cells
% 
% Output
% ------
% m : array, (3, k) where k = min(# split cells, # pooled cells)
%   each column is a match
%   first row is split-mode cell ID
%   second row is pooled-mode cell ID
%   third row is the cosine similarity score
% -------------------------------------------------------------------------
tmp_s = s;
tmp_cellid_s = double(cellid_s);
tmp_cellid_bank01 = double(cellid_bank01);
m = [];
while size(tmp_s, 1)>1 && size(tmp_s, 2)>1
    [i,j] = find(tmp_s==max(tmp_s(:)));
    m = [m, [tmp_cellid_s(i);tmp_cellid_bank01(j);tmp_s(i,j)]];
    tmp_s = tmp_s(setdiff(1:size(tmp_s,1),i)',setdiff(1:size(tmp_s,2),j));
    tmp_cellid_s = tmp_cellid_s(setdiff(1:length(tmp_cellid_s), i));
    tmp_cellid_bank01 = tmp_cellid_bank01(setdiff(1:length(tmp_cellid_bank01), j));
    
end
clear i j tmp_s S P pc
end

function s = compute_sim(cellid_s, cellid_bank01, wfs_s, wf_bank01, th)
% -------------------------------------------------------------------------
% Computes cosine similarity between given sets of waveforms; uses only 
% waveforms that exceed a threshold (microV)
% 
% Parameters
% ----------
% cellid_s : array, (1,a)
%   IDs of split cells
% cellid_bank01 : array, (1,b)
%   IDs of pooled cells
% wfs_s : array
%   average waveforms of split cells
% wf_bank01 : array
%   average waveforms of pooled cells
% th : float
%   waveforms whose peak-to-peak amplitude exceeds this threshold (unit:
%   microV) are included in similarity score computation
% 
% Output
% ------
% s : array, (a, b)
%   similarity matrix
% -------------------------------------------------------------------------

s = zeros(length(cellid_s), length(cellid_bank01));
for i = 1:length(cellid_s)
    
    S = wfs_s{1,i};

    p2p = max(S,[],2)-min(S,[],2); % peak-to-peak amplitude
    
    Ss = S(p2p>th,:);
    
    norm_S = norm(Ss(:));
    
    for j = 1:length(cellid_bank01)
        P = wf_bank01{1,j};
        Ps = P(p2p>th,:);
        
        norm_P = norm(Ps(:));
        S_dot_P = dot(Ss(:), Ps(:));
        
        s(i,j) = S_dot_P / norm_S / norm_P;
        
    end
    clear j P norm_P S_dot_P
end
clear i p2p S norm_S Ss Ps
end

function s = compute_sim2(cellid_s, cellid_bank01, wfs_s, wf_bank01, nchan)
% -------------------------------------------------------------------------
% Computes cosine similarity between given sets of waveforms; uses only 
% waveforms at channels near the peak channel (set by nchan)
% 
% Parameters
% ----------
% cellid_s : array, (1,a)
%   IDs of split cells
% cellid_bank01 : array, (1,b)
%   IDs of pooled cells
% wfs_s : array
%   average waveforms of split cells
% wf_bank01 : array
%   average waveforms of pooled cells
% nchan : int
%   number of channels around the peak channel whose waveforms will be used
%   to compute the similairty scores
% 
% Output
% ------
% s : array, (a, b)
%   similarity matrix
% -------------------------------------------------------------------------

s = zeros(length(cellid_s), length(cellid_bank01));
for i = 1:length(cellid_s)
    
    S = wfs_s{1,i};
    
    p2p = max(S,[],2)-min(S,[],2); % peak-to-peak amplitude
    % Channel around peak channel (but flanked by 1 and 384)
    chans = max([1,find(p2p==max(p2p))-nchan]):min([384,find(p2p==max(p2p))+nchan]);
    
    Ss = S(chans,:);
    
    for j = 1:length(cellid_bank01)
        P = wf_bank01{1,j};
        Ps = P(chans,:);
        
        s(i,j) = dot(Ss(:), Ps(:)) / norm(Ss(:)) / norm(Ps(:));
        
    end
    clear j P norm_P S_dot_P
end
clear i p2p S norm_S Ss Ps
end

function pc = get_pc(S, P)
% computes pooling coefficient given set of waveforms S and P
pc = zeros(1,size(S,1));
for i= 1:size(S,1)        
    wf_is = S(i,:);
    wf_ip = P(i,:);
    pc(i) = wf_is(:)\wf_ip(:); 
end
clear i
end
