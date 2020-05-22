clear
close all

% specify paths and projects
DataPath = '../dat/';
FigPath = '../fig/act_comparisons/';
mkdir(FigPath);
projectGtMutES1 = 'eveGtMut_NullS1_normAP';
projectGtMutWT = 'eveGtMut_WTS1_normAP';
% projectWT = 'eveWT';

% load data
load([DataPath projectGtMutES1 '\trace_struct.mat'])
traces_gt_es1 = trace_struct;
load([DataPath projectGtMutWT '\trace_struct.mat'])
traces_gt_wt = trace_struct;
% load([DataPath projectWT '\trace_struct.mat'])
% traces_wt = trace_struct;
clear trace_struct

% compare raw fluorescence values
fluo_bins = linspace(0,6e5,100);
fluo_hist_fig = figure;
hold on
% histogram([traces_wt.fluo],fluo_bins,'Normalization','probability')
histogram([traces_gt_wt.fluo],fluo_bins,'Normalization','probability')
histogram([traces_gt_es1.fluo],fluo_bins,'Normalization','probability')
legend('wts1-gts2','es1-gts2')
xlabel('fluorescence')
ylabel('share')
set(gca,'Fontsize',14)
saveas(fluo_hist_fig,[FigPath 'fluo_hist.png'])
%%
% Generat heatmaps of activity over time
ap_index = 15:60;
time_index = 0:50;
fluo_hm_gt_wt = NaN(numel(time_index),numel(ap_index));
fluo_hm_gt_es1 = NaN(numel(time_index),numel(ap_index));
frac_hm_gt_wt = zeros(numel(time_index),numel(ap_index));
frac_hm_gt_es1 = zeros(numel(time_index),numel(ap_index));

time_vec_gt_wt = [traces_gt_wt.time]/60;
time_vec_gt_es1 = [traces_gt_es1.time]/60;

ap_vec_gt_wt = [traces_gt_wt.APPosParticle];
ap_vec_gt_es1 = [traces_gt_es1.APPosParticle];

fluo_vec_gt_wt = [traces_gt_wt.fluo];
fluo_vec_gt_es1 = [traces_gt_es1.fluo];
% 
% [N,Xedges,Yedges,binX_wt,binY_wt] = histcounts2(time_vec_gt_wt,ap_vec_gt_wt,ap_index,time_index);
% [N,Xedges,Yedges,binX_es1,binY_es1] = histcounts2(time_vec_gt_wt,ap_vec_gt_wt,ap_index,time_index);

for a = 1:numel(ap_index)
    for t = 1:numel(time_index)
        at_filter_wt = round(time_vec_gt_wt) == time_index(t) & round(ap_vec_gt_wt) == ap_index(a);
        at_filter_es1 = round(time_vec_gt_es1) == time_index(t) & round(ap_vec_gt_es1) == ap_index(a);
        % fluo
        fluo_hm_gt_wt(t,a) = nanmean(fluo_vec_gt_wt(at_filter_wt)/(1e5));
        fluo_hm_gt_es1(t,a) = nanmean(fluo_vec_gt_es1(at_filter_es1)/(1e5));
        % fraction active
        frac_hm_gt_wt(t,a) = sum(at_filter_wt);
        frac_hm_gt_es1(t,a) = sum(at_filter_es1);
    end
end
        
        
%% fluorescence
cmap1 = flipud(brewermap(128,'Spectral'));
cmap2 = flipud(brewermap(128,'RdBu'));
% close all

wt_fig = figure;
colormap(cmap1);
p = pcolor(flipud(fluo_hm_gt_wt));
p.EdgeAlpha = 0;
cb = colorbar;
ylabel(cb,'spot intensity (au)')
xlabel('ap position (%)')
ylabel('minutes into nc14')
set(gca,'Fontsize',14)
set(gca,'YTick',0:10:50,'Yticklabels',fliplr(0:10:50))
set(gca,'XTick',(15:10:60) - 14,'Xticklabels',15:10:60)
caxis([0 3])
saveas(wt_fig,[FigPath 'fluo_s1wt_s2gt.png'])

es1_fig = figure;
colormap(cmap1);
p = pcolor(flipud(fluo_hm_gt_es1));
p.EdgeAlpha = 0;
cb = colorbar;
ylabel(cb,'spot intensity (au)')
xlabel('ap position (%)')
ylabel('minutes into nc14')
set(gca,'Fontsize',14)
set(gca,'YTick',0:10:50,'Yticklabels',fliplr(0:10:50))
set(gca,'XTick',(15:10:60) - 14,'Xticklabels',15:10:60)
caxis([0 3])
saveas(es1_fig,[FigPath 'fluo_s1wt_s2gt.png'])

diff_fig = figure;
colormap(cmap2);
imagesc(fluo_hm_gt_es1-fluo_hm_gt_wt);
% p.EdgeAlpha = 0;
cb = colorbar;
ylabel(cb,'spot intensity (au)')
xlabel('ap position (%)')
ylabel('minutes into nc14')
set(gca,'Fontsize',14)
set(gca,'YTick',0:10:50,'Yticklabels',fliplr(0:10:50))
set(gca,'XTick',(15:10:60) - 14,'Xticklabels',15:10:60)
caxis([-1.5 1.5])
saveas(diff_fig,[FigPath 'fluo_s1wt_s2gt.png'])

%% approximate fractional activity
cmap1 = flipud(brewermap(128,'Spectral'));
close all
wt_norm_factor = prctile(frac_hm_gt_wt(:),99);
es1_norm_factor = prctile(frac_hm_gt_es1(:),99);


wt_fig = figure;
colormap(cmap1);
p = pcolor(flipud(frac_hm_gt_wt / wt_norm_factor));
p.EdgeAlpha = 0;
cb = colorbar;
xlabel('ap position (%)')
ylabel('minutes into nc14')
set(gca,'Fontsize',14)
ylabel(cb,'fraction active')
set(gca,'YTick',0:10:50,'Yticklabels',fliplr(0:10:50))
set(gca,'XTick',(15:10:60) - 14,'Xticklabels',15:10:60)
caxis([0 1])
saveas(wt_fig,[FigPath 'fraction_active_s1wt_s2gt.png'])

es1_fig = figure;
colormap(cmap1);
p = pcolor(flipud(frac_hm_gt_es1 / es1_norm_factor));
p.EdgeAlpha = 0;
cb = colorbar;
xlabel('ap position (%)')
ylabel('minutes into nc14')
set(gca,'Fontsize',14)
ylabel(cb,'fraction active')
set(gca,'YTick',0:10:50,'Yticklabels',fliplr(0:10:50))
set(gca,'XTick',(15:10:60) - 14,'Xticklabels',15:10:60)
caxis([0 1])
saveas(es1_fig,[FigPath 'fraction_active_s1null_s2gt.png'])

diff_fig = figure;
colormap(cmap2);
imagesc(frac_hm_gt_es1 / es1_norm_factor - frac_hm_gt_wt / wt_norm_factor);
% p.EdgeAlpha = 0;
cb = colorbar;
xlabel('ap position (%)')
ylabel('minutes into nc14')
ylabel(cb,'fraction active')
set(gca,'Fontsize',14)
set(gca,'YTick',0:10:50,'Yticklabels',fliplr(0:10:50))
set(gca,'XTick',(15:10:60) - 14,'Xticklabels',15:10:60)
caxis([-0.8 0.8])
saveas(diff_fig,[FigPath 'fraction_active_s1null_diff_s1wt.png'])