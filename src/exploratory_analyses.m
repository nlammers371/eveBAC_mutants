clear 
close all
% set ID variables and file paths
addpath('utilities')
project = '20190613_eve1spot';
% project = '20190613_eveGtMut_eS1';
dataPath = ['../dat/' project '/'];
figPath = ['../fig/' project '/'];
mkdir(figPath);
% load 
load([dataPath 'nucleus_struct.mat'])

% extract useful vectors
ap_vec = [nucleus_struct.apPos]*100;
time_vec = [nucleus_struct.time]/60;
fluo_vec = [nucleus_struct.fluo]*31.5/30;
if size(fluo_vec,1) == 2
    fluo_vec = [fluo_vec(1,:)];% fluo_vec(2,:)];
%     time_vec = [time_vec time_vec];
%     ap_vec = [ap_vec ap_vec];
end
fluo_vec(isnan(fluo_vec)) = 0;
% take a look at some basic spatiotemporal trends
ap_grid = 1:100;
time_grid = 1:60;
mf_grid = NaN(numel(time_grid),numel(ap_grid));
% populate grid
for a = 1:numel(ap_grid)
    for t = 1:numel(time_grid)
        ap = ap_grid(a);
        time = time_grid(t);
        tr_ft = round(ap_vec)==ap & round(time_vec)==time;
        mf_grid(t,a) = nanmean(fluo_vec(tr_ft));
    end
end

% make figure
hm_fig = figure;
p = pcolor(flipud(mf_grid));
p.EdgeAlpha = .2;
ylim([11,60])
xlim([15,60])
set(gca,'ytick',0:5:60,'yticklabels',fliplr(0:5:60))
set(gca,'xtick',0:5:60,'xticklabels',0:5:60)
ylabel('minutes into nc14')
xlabel('AP position')
h = colorbar;
ylabel(h,'transcription rate (AU)')
colormap(flipud(brewermap(128,'RdYlBu'))) 
set(gca,'Fontsize',12);
caxis([0 2e5])
saveas(hm_fig,[figPath 'activity_heatmap.tif'])