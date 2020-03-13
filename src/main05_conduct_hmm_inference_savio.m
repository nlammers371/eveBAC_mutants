% main05_conduct_hmm_inference(project,varargin)
%
% DESCRIPTION
% Script to conduct HMM inference
%
% ARGUMENTS
% project: master ID variable 
%
% modelPath: file path to folder containing hmmm scripts
%
% w: Integer corresponding number of time steps for Pol II to transcribe
% gene
%
%
% OPTIONS
% dropboxFolder: Path to data folder where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% savio: if 1, indicates we are running inference on savio cluster
% K: number of states
% minDp: min data points needed to be included in inferece
%
% OUTPUT: nucleus_struct_protein: compiled data set with protein samples

% function output = main05_conduct_hmm_inference_savio%(project, DropboxFolder, varargin)
clear
close all
warning('off','all') %Shut off Warnings

% basic inputs
project = 'eveGtMut_WTS1_normAP';
% project = 'eveGtS2-WT';
DataPath = ['../dat/' project '/'];
% default path to model scripts
modelPath = './utilities';
savioFlag = 1;
awsFlag = 0;
K = 3;
w = 7;
ec_flag = false;
dpBootstrap = 1;
nBoots = 2;
stripe_bin_flag = true;
SampleSize = 3000;
maxWorkers = 16;
t_start = 20;
minDP = 2*w;
if contains(project,'NullS1')
    stripe_set = 2:7;
else
    stripe_set = 1:7;
end

%%%%% These options generally remain fixed 
n_localEM = 24; % set num local runs
n_steps_max = 500; % set max steps per inference
eps = 1e-4; % set convergence criteria
min_dp_per_inf = 1000; % inference will be aborted if fewer present 
%%%%%%%%%%%%%%

addpath(modelPath); % Route to utilities folder
% check that we have proper fields
if ~dpBootstrap
    warning('Bootstrap option not selected. Setting nBoots to 1')
    nBoots = 1;
end

%-------------------------------System Vars-------------------------------%
alphaFrac = 1302 / 6444;
alpha = alphaFrac*w;
%----------------------------Set Write Paths------------------------------%
d_type = '';
if dpBootstrap
    d_type = '_dp';
end
load([DataPath '/trace_struct.mat'],'trace_struct') % load data
% Set write path (inference results are now written to external directory)

% set write path
if savioFlag
    out_prefix = '/global/scratch/nlammers/'; %hmmm_data/inference_out/';
else    
    out_prefix = '../out/';
end
outDir = [out_prefix project '/K' num2str(K) 'w' num2str(w) '_ec' num2str(ec_flag)];
mkdir(outDir);

Tres = trace_struct(1).TresInterp; % Time Resolution
% filter for quality traces of sufficient length
trace_struct_filtered = [];
for i = 1:length(trace_struct)
    temp = struct;
    time = trace_struct(i).time_interp;
    time_raw = trace_struct(i).time;
    fluo = trace_struct(i).fluo_interp;  
    time_ft = time / 60 >= t_start;
    time_ft_raw = time_raw / 60 >= t_start;
    if sum(time_ft) >= minDP
        temp.fluo = fluo(time_ft);
        temp.time = time(time_ft);     
        temp.qc_flag = trace_struct(i).qc_flag;
        stripe_id = mode(trace_struct(i).Stripe(time_ft_raw));
        ap_pos = nanmean(trace_struct(i).APPosParticleNorm);
        if ap_pos < 37 && stripe_id == 2
            stripe_id = 1.5;
        end
        temp.Stripe = stripe_id;   
        temp.ec_flag = ~ismember(temp.Stripe,stripe_set);
        temp.MeanFluo = nanmean(fluo(time_ft));
        temp.ParticleID = trace_struct(i).ParticleID;    
        temp.N = sum(time_ft);
        trace_struct_filtered = [trace_struct_filtered temp];    
    end
end
trace_struct_filtered = trace_struct_filtered([trace_struct_filtered.qc_flag]==1);

% generate fluorescence bins
if ~ec_flag
    stripe_id_vec = [trace_struct_filtered.Stripe];
else
    stripe_id_vec = [trace_struct_filtered.ec_flag];
end
stripe_index = unique(stripe_id_vec);
fluo_bin_cell = cell(1,numel(stripe_index));
for s = 1:numel(stripe_index)
    % filter for stripe
    stripe_ft = stripe_index(s) == stripe_id_vec;
    stripe_indices = find(stripe_ft);
    fluo_vec = [trace_struct_filtered(stripe_ft).MeanFluo];
    nTotal = sum([trace_struct_filtered(stripe_ft).N]);
    nBins = ceil(nTotal/SampleSize)+1;
    prctile_vec = linspace(0.025,.975,nBins);
    % get quantile bins
    fluo_quantiles = quantile(fluo_vec,prctile_vec);
    fluo_bin_cell{s} = fluo_quantiles;
    fluo_id_vec_temp = discretize(fluo_vec,fluo_quantiles);
    iter = 1;
    for ind = stripe_indices
        trace_struct_filtered(ind).FluoBin = fluo_id_vec_temp(iter);
        trace_struct_filtered(ind).FluoBins = fluo_quantiles;
        iter = iter + 1;
    end
end

fluo_id_vec = [trace_struct_filtered.FluoBin];
%%% Conduct Inference
% iterate through designated groups
rng('shuffle')
stripe_index = 3;
for s = 1:numel(stripe_index)
    fluo_bins = fluo_bin_cell{s};
    for t = 1:length(fluo_bins)-1
        iter_filter = fluo_id_vec == t & stripe_id_vec == stripe_index(s);
        for b = 1:nBoots
            iter_start = now;
            local_struct = struct;    
            output = struct;

            % Use current time as unique inference identifier 
            inference_id = num2str(round(10e5*now));

            % Generate filenames            
            fName_sub = ['hmm_results_t' inference_id '.mat'];                
            out_file = [outDir '/' fName_sub];  

            % Extract fluo_data        
            inference_set = trace_struct_filtered(iter_filter);        
            skip_flag = 0;
            set_size = length([inference_set.fluo]);                 
            if isempty(inference_set)
                skip_flag = 1;
            elseif set_size < min_dp_per_inf                    
                skip_flag = 1;                    
            end
            if skip_flag
                warning('Too few data points. Skipping')                
            else 
                sample_index = 1:length(inference_set);
                if dpBootstrap                        
                    ndp = 0;    
                    sample_ids = [];                    
                    %Reset bootstrap size to be on order of set size for small bins
                    if set_size < SampleSize
                        SampleSize = ceil(set_size/100)*100;
                    end
                    while ndp < SampleSize
                        tr_id = randsample(sample_index,1);
                        sample_ids = [sample_ids tr_id];
                        ndp = ndp + length(inference_set(tr_id).time);
                    end
                    fluo_data = cell([length(sample_ids), 1]);    
                    time_data = cell([length(sample_ids), 1]);    
                    sample_particles = [inference_set(sample_ids).ParticleID];
                    for tr = 1:length(sample_ids)
                        fluo_data{tr} = inference_set(sample_ids(tr)).fluo;                    
                        time_data{tr} = inference_set(sample_ids(tr)).time;                    
                    end            
                else % Take all relevant traces if not bootstrapping
                    fluo_data = cell([length(inference_set), 1]);    
                    time_data = cell([length(inference_set), 1]);    
                    for tr = 1:length(inference_set)
                        fluo_data{tr} = inference_set(tr).fluo;
                        time_data{tr} = inference_set(tr).time;                    
                    end
                    sample_particles = [inference_set.ParticleID];
                end
                % Random initialization of model parameters
                param_init = initialize_random (K, w, fluo_data);
                % Approximate inference assuming iid data for param initialization                
                local_iid_out = local_em_iid_reduced_memory(fluo_data, param_init.v, ...
                                    param_init.noise, K, w, alpha, 500, 1e-4);
                noise_iid = 1/sqrt(exp(local_iid_out.lambda_log));
                v_iid = exp(local_iid_out.v_logs);            
                p = gcp('nocreate');
                if isempty(p)
                    parpool(maxWorkers); %6 is the number of cores the Garcia lab server can reasonably handle per user.
                elseif p.NumWorkers > maxWorkers
                    delete(gcp('nocreate')); % if pool with too many workers, delete and restart
                    parpool(maxWorkers);
                end
                parfor i_local = 1:n_localEM % Parallel Local EM                
                    % Random initialization of model parameters
                    param_init = initialize_random_with_priors(K, noise_iid, v_iid);
                    % Get Intial Values
                    pi0_log_init = log(param_init.pi0);
                    A_log_init = log(param_init.A);
                    v_init = param_init.v;                        
                    noise_init = param_init.noise;
                    %--------------------LocalEM Call-------------------------%
                    local_out = local_em_MS2_reduced_memory_truncated(fluo_data, ...
                        v_init, noise_init, pi0_log_init', A_log_init, K, w, ...
                        alpha, n_steps_max, eps);                    
                    %---------------------------------------------------------%                
                    % Save Results 
                    local_struct(i_local).inference_id = inference_id;
                    local_struct(i_local).subset_id = i_local;
                    local_struct(i_local).logL = local_out.logL;                
                    local_struct(i_local).A = exp(local_out.A_log);
                    local_struct(i_local).v = exp(local_out.v_logs).*local_out.v_signs;
                    local_struct(i_local).r = exp(local_out.v_logs).*local_out.v_signs / Tres;                                
                    local_struct(i_local).noise = 1/exp(local_out.lambda_log);
                    local_struct(i_local).pi0 = exp(local_out.pi0_log);
                    local_struct(i_local).total_steps = local_out.n_iter;               
                    local_struct(i_local).soft_struct = local_out.soft_struct;               
                end
                [~, max_index] = max([local_struct.logL]); % Get index of best result                    
                % Save parameters from most likely local run
                output.pi0 =local_struct(max_index).pi0;                        
                output.r = local_struct(max_index).r(:);          
                output.noise = local_struct(max_index).noise;
                output.A = local_struct(max_index).A(:);
                output.A_mat = local_struct(max_index).A;            
                % get soft-decoded structure
                output.soft_struct = local_struct(max_index).soft_struct;
                % Info about run time
                output.total_steps = local_struct(max_index).total_steps;                                  
                output.total_time = 100000*(now - iter_start);            
                % other inference characteristics                                
                output.dp_bootstrap_flag = dpBootstrap;   
                output.iter_id = b;                     
                output.particle_ids = sample_particles;
                output.FluoBin = t;
                output.Stripe = stripe_index(s);                
                if dpBootstrap                                    
                    output.N = ndp;
                end
                output.w = w;
                output.alpha = alpha;
                output.deltaT = Tres; 
                output.sampleSize = SampleSize; 
                % save inference data used
                output.fluo_bins = fluo_bins;
                output.fluo_data = fluo_data;
                output.time_data = time_data;
                output.ec_flag = ec_flag;
            end
            output.skip_flag = skip_flag;
            save(out_file, 'output');           
        end  
    end
end
