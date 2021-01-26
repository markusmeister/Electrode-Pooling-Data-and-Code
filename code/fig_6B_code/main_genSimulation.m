% Example code of generating simulated traces
% with amp, noise, fr as free parameters

% change to your root folder of the simulation
% Sim with adding noise, and filter the final signal+noise traces
clear all;
simroot = 'D:\Repo\Electrode-Pooling-Data-and-Code';
cd(simroot); %
% helper functions 
addpath(fullfile(simroot,'code','fig_6B_code','gen_simulation_code'))
% update seed num to generate simulations of this noise level
seednum = 1001; 
rng(seednum,'twister');
s = rng;

%% Create ground truth for scoring and poisson spike train of the simulation
firing_rate = 10; % 10 hz
simLength = 60;  % default=600 sec; 60sec for github example
numTrials = 12; % simulate 12 unit's spike time
% generate spike timing in binaries, record ground truth of that run
[ground_truth_cell, legit_spikeMat ] = createSpiketrain(firing_rate, simLength,numTrials);

% Noise settings
% 5.7, 9, 1.6 is close to filtered saline data, dial in pre filtered
% level; 
prefilt_ratio = 0.52; % scale the noise amp up by 1/ratio for post butterfilter match
N_com_std = 5.7 / prefilt_ratio; % Npx 5.7
N_bio_std = 9 / prefilt_ratio;  % 9/15
N_ele_std = 1.6 / prefilt_ratio; % 

% Spike amp setting
% match the p2p amp percentile target from Allen data vpm
target_p2p_amp = 380; 



%% load the template waveforms and foot print ratios
% load footprints template
load(fullfile(simroot,'data','simdata', 'sim_materials',...
'final_footprint.mat'),'final_footprint');

% load spike waveforms templates
load(fullfile(simroot, 'data','simdata','sim_materials',...
    'final_12_template.mat'),'final_12_template');

%% update dir_name for every exp if change parameter !
savedir = fullfile(simroot,'data', 'simdata', 'sim_example') ;
dir_name = sprintf('%d_uV_10hz_noise_5_9_2_seed_%d_sim_example',...
    target_p2p_amp , seednum);

if exist(fullfile(savedir, dir_name),'dir')~=7
    [status,msg] = mkdir(fullfile(savedir, dir_name));
    assert(status,'Could not create directory ''%s'': %s',fullfile(savedir, dir_name),msg);
end
%% scale the templates to match the input amp

scale_ratio = matchAmpscaler( final_12_template, target_p2p_amp );
scale_arr = ones(1, length(final_12_template));
scale_arr(1:end)= scale_ratio; 
final_12_template = helper_scale_template_bank_func(final_12_template,scale_arr);
%% create voltage traces
chno_num = 4; % # of electrodes
M = 12; % max banks # to pool, do not exceed dimension of template cell & fr time template

%Filter to match recording system bandpass
%set the butterworth filter
samplingrate=30e3;
filt_band = [300,10000];
filt_order = 3;
filt_passtype = 'bandpass';
[b,a] = butter(filt_order,filt_band/(samplingrate/2),filt_passtype);

for n = 1:M % outer loop for providing place holder and saving the result of pooling n sets of tetrodes, with Max of M
    
    % place holder for saving the simulated traces
    spikeWaveform_by_t = zeros(chno_num, size(legit_spikeMat,2)+ size(final_12_template{1,1},2) -1); 
    rng(s); % make noise deterministic within the same seed run
    
for m = 1:n % indexting banks of the tetrodes banks to pool
    tmp_max_wave_form = final_12_template{1,m};% load scaled template
    
    %convolve wave form with spike with scaled foot print
    Conv_waveform_by_t_ = conv2(tmp_max_wave_form,  final_footprint(:, m) * legit_spikeMat(m, :));
    
    % create private noise
    % electrode noise are independent for the 4 electrodes
    ele_noise = normrnd(0, N_ele_std, size(Conv_waveform_by_t_ ));
    % bio noise correlated among the 4 electrode
    bio_noise = repmat(normrnd(0, N_bio_std, [1, size(Conv_waveform_by_t_ ,2 )]) , 4,1);
    
    % add private noises
    Conv_waveform_by_t_ = Conv_waveform_by_t_ + ele_noise + bio_noise;
    
    % update the place holder with this pool of tetrode, scaled by the
    % fraction of # of pool
    spikeWaveform_by_t = spikeWaveform_by_t + (Conv_waveform_by_t_ / n); 
end
    
    rng(seednum+1); 
    % add common noise
    common_noise = normrnd(0, N_com_std, size(spikeWaveform_by_t ));
    spikeWaveform_by_t = spikeWaveform_by_t + common_noise;
    
    %Filtfilt the voltagetrace
    spikeWaveform_by_t = filtfilt(b, a, spikeWaveform_by_t')';
    
    % Bandpass filtered RMS measurement
    [~, ele_RMS] = matchbandFilt(filtfilt(b, a, ele_noise')' );
    [~, bio_RMS] = matchbandFilt(filtfilt(b, a, bio_noise')');
    [~, common_RMS] = matchbandFilt(filtfilt(b, a, common_noise')');     
        
    %save binary .dat file in its specific folder for Kilosorting
    simPoolDir = sprintf('simpool%s', num2str(n));
    if exist(fullfile(savedir, dir_name, simPoolDir),'dir')~=7
        [status,msg] = mkdir(fullfile(savedir, dir_name, simPoolDir));
        assert(status,'Could not create directory ''%s'': %s',fullfile(savedir, dir_name),msg);
    end
    fidW     = fopen(fullfile(savedir ,dir_name, simPoolDir, ['testmergepool_' num2str(n) '.dat']), 'w');
    fwrite(fidW, int16(spikeWaveform_by_t) , 'int16'); % save matrix to .dat
    fclose(fidW); % all done
    
    % save a visualization of the sim sec
    fig = figure;
    ax =axes;
    hold on
    % visualize the first .5 sec (15000 data points)
    for plot_i = 1:chno_num
        plot(1:15000, spikeWaveform_by_t(plot_i, 1:15000)+ plot_i*200 , 'LineWidth', 2, 'color', 'k')
    end
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8 6])
    fig_name = (sprintf('pool total %d', m));
    title(fig_name)
    print('-dpng', fullfile(savedir, dir_name, sprintf('%s.png',fig_name)), '-r100');
    
    cla(ax)
    close(fig)
    
end

%% save construct the param 
Param_pool_sim = struct;
Param_pool_sim.seednum = seednum;
Param_pool_sim.N_com_std = common_RMS; % the post band passed RMS 
Param_pool_sim.N_bio_std = bio_RMS;
Param_pool_sim.N_ele_std = ele_RMS;
Param_pool_sim.chno_num = chno_num; % sim tetrode num
Param_pool_sim.scale_ratio = scale_ratio;
Param_pool_sim.targetp2pAmp = target_p2p_amp;
%Param_pool_sim.rev_opt = rev_opt;
% save the params & Ground Truth spiketimes
save(fullfile(savedir, dir_name, 'Param_pool_sim.mat'), 'Param_pool_sim')
save(fullfile(savedir, dir_name, sprintf('ground_truth_cell_%d_seed_%d_fr', seednum, firing_rate)), 'ground_truth_cell')

