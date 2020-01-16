clear 
close all

% set path to inference results
project = 'eveGtS2-NullS1';
% project = 'eveGtS2-WT';
w = 6;
K = 3;
ResultsPath = ['../out/' project '/w' num2str(w) '_K' num2str(K) '/'];
FileList = dir([ResultsPath '/*.mat']);
% make figure path
FigPath = ['../fig/' project '/w' num2str(w) '_K' num2str(K) '/'];
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

%% Extract average parameter trends
Tres = inference_results(1).deltaT;
FluoBinVec = [inference_results.FluoBin];
FluoIndex = unique(FluoBinVec);
StripeBinVec = [inference_results.Stripe];
StripeIndex = unique(StripeBinVec);
FluoBins = inference_results(1).FluoBins(1:end-1) + diff(inference_results(1).FluoBins)/2;
% initialize arrays
amp_mean_array = NaN(numel(FluoIndex),numel(StripeIndex));
amp_ste_array = NaN(numel(FluoIndex),numel(StripeIndex));

dur_mean_array = NaN(numel(FluoIndex),numel(StripeIndex));
dur_ste_array = NaN(numel(FluoIndex),numel(StripeIndex));

freq_mean_array = NaN(numel(FluoIndex),numel(StripeIndex));
freq_ste_array = NaN(numel(FluoIndex),numel(StripeIndex));

for s = 1:numel(StripeIndex)
    for f = 1:numel(FluoIndex)
        iter_indices = find(FluoBinVec==FluoIndex(f) & StripeBinVec==StripeIndex(s));
        amp_vec = NaN(1,numel(iter_indices));
        dur_vec = NaN(1,numel(iter_indices));
        freq_vec = NaN(1,numel(iter_indices));
        for i = 1:numel(iter_indices)
            [r,ri] = sort(inference_results(iter_indices(i)).r);
            A = inference_results(iter_indices(i)).A_mat(ri,ri);
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
        amp_mean_array(f,s) = nanmean(amp_vec);
        amp_ste_array(f,s) = nanstd(amp_vec);
        
        dur_mean_array(f,s) = nanmean(dur_vec);
        dur_ste_array(f,s) = nanstd(dur_vec);
        
        freq_mean_array(f,s) = nanmean(freq_vec);
        freq_ste_array(f,s) = nanstd(freq_vec);
    end
end
        
        
%% Plot
close all
cm1 = flipud(brewermap(9,'Blues'));
cm2 = brewermap(9,'Reds');
ectopic_ids = [1 3];
MarkerSize = 60;

stripe_lgd = {'"stripe 0"','stripe 1','1-2 interstripe','stripe 2', 'stripe 3', 'stripe 4'};  
% plot initiation rate 
init_fig = figure;
hold on
sc = [];
for s = 1:numel(StripeIndex)
    e = errorbar(FluoBins,amp_mean_array(:,s)*60,amp_ste_array(:,s)*60,'o','Color','black');
    e.CapSize = 0;
    cmap = cm2;
    if ismember(s,ectopic_ids)
        cmap = cm1(3:end,:);
    end
    sc = [sc scatter(FluoBins,amp_mean_array(:,s)*60,MarkerSize,'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black')];
end

xlabel('fluorescence (au)')
ylabel('burst amplitude (au/min)')
legend(sc,stripe_lgd{:},'Location','southeast','Fontsize',12)
grid on
box on
set(gca,'FontSize',14)
saveas(init_fig,[FigPath 'init_rate_trends.png'])   
saveas(init_fig,[FigPath 'init_rate_trends.pdf'])   

%% plot burst frequency
freq_fig = figure;
hold on
sc = [];
for s = 1:numel(StripeIndex)
    x_axis = FluoBins+rand(size(FluoBins))*1000-2000;
    e = errorbar(x_axis,freq_mean_array(:,s)*60,freq_ste_array(:,s)*60,'o','Color','black');
    e.CapSize = 0;
    cmap = cm2;
    if ismember(s,ectopic_ids)
        cmap = cm1(3:end,:);
    end
    sc = [sc scatter(x_axis,freq_mean_array(:,s)*60,MarkerSize,'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black','MarkerFaceAlpha',1)];
end

xlabel('fluorescence (au)')
ylabel('burst frequency (1/min)')
% legend(sc,stripe_lgd{:},'Location','northwest','Fontsize',12)
grid on
box on
set(gca,'FontSize',14)
saveas(freq_fig,[FigPath 'burst_freq_trends.png'])

%% plot burst duration
dur_fig = figure;
hold on
sc = [];
for s = 1:numel(StripeIndex)
    e = errorbar(FluoBins,dur_mean_array(:,s)/60,dur_ste_array(:,s)/60,'o','Color','black');
    e.CapSize = 0;
    cmap = cm2;
    if ismember(s,ectopic_ids)
        cmap = cm1(3:end,:);
    end
    sc = [sc scatter(FluoBins,dur_mean_array(:,s)/60,MarkerSize,'MarkerFaceColor',cmap(s,:),'MarkerEdgeColor','black')];
end

xlabel('fluorescence (au)')
ylabel('burst duration (min)')
grid on
box on
% legend(sc,stripe_lgd{:},'Location','northwest','Fontsize',12)
set(gca,'FontSize',14)
saveas(dur_fig,[FigPath 'burst_dur_trends.png'])


