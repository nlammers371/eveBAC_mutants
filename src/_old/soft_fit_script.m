% script to perform soft fits to traces using existing HMM parameters
clear
close all
% ID variables
% project = '20190514_eveBAC2spot'; 
project = '20190613_eveGtMut_eS1';
dataPath = ['../dat/' project '/'];
figPath = ['../fig/' project '/'];
mkdir(figPath);
% load trace data
load([dataPath 'nucleus_struct.mat'])
% load HMM results
hmm_path = "D:\Data\Nick\projects\all_the_stripes_cd\dat\eve7stripes_inf_2018_04_28\w7_t20_alpha14_f1_cl1_no_ends1_tbins1\K3_summary_stats\hmm_results_t_window30_t_inf40.mat";
load(hmm_path)

% define fit params
Tres_inf = 20;
Tres_dat = nucleus_struct(1).TresInterp;
alpha = 1.4;
w = 7;
K = 3;
use_ind = 1;
A_log = log(reshape(hmm_results(use_ind).A_mean,K,K));
v = hmm_results(use_ind).initiation_mean;
sigma = sqrt(hmm_results(use_ind).noise_mean);
pi0_log = log([1 1 1]/3); 
eps = 1e-4;
% generate cell arrays of times and traces
minDP = 15;
particle_index = [nucleus_struct.ParticleID];
qc_indices = find([nucleus_struct.inference_flag]);
fluo_values = cell(numel(qc_indices),1);
time_values = cell(numel(qc_indices),1);
ap_values = cell(numel(qc_indices),1);
rm_indices = [];
%%
for i = 1:numel(qc_indices)
    fluo_raw = nucleus_struct(qc_indices(i)).fluo_interp;
    time_raw = nucleus_struct(qc_indices(i)).time_interp;
    time_raw_raw = nucleus_struct(qc_indices(i)).time;
    ap_raw = nucleus_struct(qc_indices(i)).apPos;
    start_i = find(~isnan(fluo_raw),1);
    stop_i = find(~isnan(fluo_raw),1,'last');
    fluo_raw = fluo_raw(start_i:stop_i);
    time_raw = time_raw(start_i:stop_i);
    time_new = time_raw(1):Tres_inf:time_raw(end);
    fluo_new = interp1(time_raw,fluo_raw,time_new);   
    ap_new = interp1(time_raw_raw,ap_raw,time_new);   
    if numel(fluo_new) < minDP
        rm_indices = [rm_indices i];
    end
    fluo_values{i} = fluo_new;
    time_values{i} = time_new;
    ap_values{i} = ap_new;
end 
ft = ~ismember(1:numel(fluo_values),rm_indices);
fluo_values = fluo_values(ft);
ap_values = ap_values(ft);
qc_indices = qc_indices(ft);
time_values = time_values(ft);
% fit
disp('conducting single trace fits...')
tic 
local_em_outputs = local_em_MS2_reduced_memory (fluo_values, ...
              v', sigma, pi0_log, A_log, K, w, alpha, 1, eps);
toc
soft_fit_struct = local_em_outputs.soft_struct;
soft_fit_struct.time_cell = time_values;
soft_fit_struct.fluo_cell = fluo_values;
soft_fit_struct.ap_cell = ap_values;
soft_fit_struct.particle_index = particle_index(qc_indices);
save([dataPath 'soft_fit_struct.mat'],'soft_fit_struct')