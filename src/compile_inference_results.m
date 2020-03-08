clear 
close all

% set path to inference results
project = 'eveGtS2-NullS1';
% project = 'eveGtS2-WT';
w = 7;
K = 3;
ec_flag = true;
ResultsPath = ['../out/' project '/K' num2str(K) 'w' num2str(w) '_ec' num2str(ec_flag) '/'];
FileList = dir([ResultsPath '/*.mat']);
% make figure path
FigPath = ['../fig/' project '/w' num2str(w) '_K' num2str(K) '_ec' num2str(ec_flag) '/'];
mkdir(FigPath);


% read files
inference_results = struct;
iter = 1;
for f = 1:numel(FileList)
    load([ResultsPath FileList(f).name])
    fnames = fieldnames(output);
    if numel(fnames) > 5
        for fn = 1:numel(fnames) 
            inference_results(iter).(fnames{fn}) = output.(fnames{fn});
        end
        iter = iter + 1;
    end
end

% Extract average parameter trends
Tres = inference_results(1).deltaT;
stripe_bin_vec = [inference_results.Stripe];

stripe_index = unique(stripe_bin_vec);

% initialize arrays
fluo_axis_cell = cell(1,numel(stripe_index));

amp_mean_cell = cell(1,numel(stripe_index));
amp_ste_cell = cell(1,numel(stripe_index));

dur_mean_cell = cell(1,numel(stripe_index));
dur_ste_cell = cell(1,numel(stripe_index));

freq_mean_cell = cell(1,numel(stripe_index));
freq_ste_cell = cell(1,numel(stripe_index));

for s = 1:numel(stripe_index)
    stripe_results = inference_results(stripe_index(s)==stripe_bin_vec);
    fluo_bin_vec = [stripe_results.FluoBin];
    fluo_index = unique(fluo_bin_vec);
    fluo_axis = (stripe_results(1).fluo_bins(1:end-1)+stripe_results(1).fluo_bins(2:end))/2;
    fluo_axis_cell{s} = fluo_axis;

    for f = 1:numel(fluo_index)
        iter_indices = find(fluo_bin_vec==fluo_index(f));
        amp_vec = NaN(1,numel(iter_indices));
        dur_vec = NaN(1,numel(iter_indices));
        freq_vec = NaN(1,numel(iter_indices));
        for i = 1:numel(iter_indices)
            [r,ri] = sort(stripe_results(iter_indices(i)).r);
            A = stripe_results(iter_indices(i)).A_mat(ri,ri);
%             [V,D] = eig(A);
%             [~,di] = max(diag(D));
%             ss = V(:,di) / sum(V(:,di));
%             % mean initiation rate
%             amp_vec(i) = (ss(2)*r(2) + ss(3)*r(3)) / (ss(2) + ss(3));
            % calculate effective 2 state transition matrix
            R = logm(A) / Tres;
            if ~isreal(R) || sum(R(:)<0) > K
                out = prob_to_rate_fit_sym(A, Tres, 'gen', .005, 1);            
                R = out.R_out;     
            end
            [V,D] = eig(R);
            [~,di] = max(diag(D));
            ss = V(:,di) / sum(V(:,di));
            % mean initiation rate
            amp_vec(i) = (ss(2)*r(2) + ss(3)*r(3)) / (ss(2) + ss(3));        
            % convert o rates
%             R_eff = logm(A_eff) / Tres;
            dur_vec(i) = -1/R(1,1) *(1/ss(1) - 1);
            freq_vec(i) = -R(1,1);
        end
        % take average and ste
        amp_mean = amp_mean_cell{s};
        amp_mean_cell{s} = [amp_mean nanmean(amp_vec)];
        amp_ste = amp_ste_cell{s};
        amp_ste_cell{s} = [amp_ste nanstd(amp_vec)];
        
        dur_mean = dur_mean_cell{s};
        dur_mean_cell{s} = [dur_mean nanmean(dur_vec)];
        dur_ste = dur_ste_cell{s};
        dur_ste_cell{s} = [dur_ste nanstd(dur_vec)];
        
        freq_mean = freq_mean_cell{s};
        freq_mean_cell{s} = [freq_mean nanmean(freq_vec)];
        freq_ste = freq_ste_cell{s};
        freq_ste_cell{s} = [freq_ste nanstd(freq_vec)];
    end
end
        
        
%%% Plot
close all


MarkerSize = 60;

if ec_flag
    ectopic_ids = [2];
    stripe_lgd = {'ectopic regions','stripe regions'};  
    cm1 = flipud(brewermap(9,'Blues'));
    cm1 = repmat(cm1(5,:),9,1);
    cm2 = brewermap(9,'Reds');
    cm2 = repmat(cm2(5,:),9,1);
else
    ectopic_ids = [1 3];
    stripe_lgd = {'"stripe 0"','stripe 1','1-2 interstripe','stripe 2', 'stripe 3', 'stripe 4'};  
    cm1 = flipud(brewermap(9,'Blues'));
    cm2 = brewermap(9,'Reds');
end
% plot initiation rate 
init_fig = figure;
hold on
sc = [];
for s = 1:numel(stripe_index)
    e = errorbar(fluo_axis_cell{s},amp_mean_cell{s}*60,amp_ste_cell{s}*60,'o','Color','black');
    e.CapSize = 0;
    cmap = cm2;
    if ismember(s,ectopic_ids)
        cmap = cm1(3:end,:);
    end
    sc = [sc scatter(fluo_axis_cell{s},amp_mean_cell{s}*60,MarkerSize,'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black')];
end

xlabel('fluorescence (au)')
ylabel('burst amplitude (au/min)')
legend(sc,stripe_lgd{:},'Location','southeast','Fontsize',12)
grid on
box on
set(gca,'FontSize',14)
saveas(init_fig,[FigPath 'init_rate_trends_ec' num2str(ec_flag) '.png'])   
saveas(init_fig,[FigPath 'init_rate_trends_ec' num2str(ec_flag) '.pdf'])   

% plot burst frequency
freq_fig = figure;
hold on
sc = [];
for s = 1:numel(stripe_index)
%     x_axis = FluoBins+rand(size(FluoBins))*1000-2000;
    e = errorbar(fluo_axis_cell{s},freq_mean_cell{s}*60,freq_ste_cell{s}*60,'o','Color','black');
    e.CapSize = 0;
    cmap = cm2;
    if ismember(s,ectopic_ids)
        cmap = cm1(3:end,:);
    end
    sc = [sc scatter(fluo_axis_cell{s},freq_mean_cell{s}*60,MarkerSize,'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black','MarkerFaceAlpha',1)];
end

xlabel('fluorescence (au)')
ylabel('burst frequency (1/min)')
% legend(sc,stripe_lgd{:},'Location','northwest','Fontsize',12)
grid on
box on
set(gca,'FontSize',14)
saveas(freq_fig,[FigPath 'burst_freq_trends_ec' num2str(ec_flag) '.png'])
saveas(freq_fig,[FigPath 'burst_freq_trends_ec' num2str(ec_flag) '.pdf'])

% plot burst duration
dur_fig = figure;
hold on
sc = [];
for s = 1:numel(stripe_index)
    e = errorbar(fluo_axis_cell{s},dur_mean_cell{s}*60,dur_ste_cell{s}*60,'o','Color','black');
    e.CapSize = 0;
    cmap = cm2;
    if ismember(s,ectopic_ids)
        cmap = cm1(3:end,:);
    end
    sc = [sc scatter(fluo_axis_cell{s},dur_mean_cell{s}*60,MarkerSize,'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black')];
end

xlabel('fluorescence (au)')
ylabel('burst duration (min)')
grid on
box on
% legend(sc,stripe_lgd{:},'Location','northwest','Fontsize',12)
set(gca,'FontSize',14)
saveas(dur_fig,[FigPath 'burst_dur_trends_ec' num2str(ec_flag) '.png'])
saveas(dur_fig,[FigPath 'burst_dur_trends_ec' num2str(ec_flag) '.pdf'])


results_struct = struct;

for s = 1:numel(stripe_index)
    results_struct(s).fluo_axis = fluo_axis_cell{s};    
    results_struct(s).region = stripe_index(s);
    results_struct(s).ectopic_flag = ismember(s,ectopic_ids);
    % dur
    results_struct(s).burst_duration_mean = dur_mean_cell{s};
    results_struct(s).burst_duration_ste = dur_ste_cell{s};
    % amp
    results_struct(s).burst_amplitude_mean = amp_mean_cell{s};
    results_struct(s).burst_amplitude_ste = amp_ste_cell{s};
    % freq
    results_struct(s).burst_frequency_mean = freq_mean_cell{s};
    results_struct(s).burst_frequency_ste = freq_ste_cell{s};
end

save(['../out/' project '/K' num2str(K) 'w' num2str(w) '_ec' num2str(ec_flag) '_results.mat'],'results_struct')
    
    

