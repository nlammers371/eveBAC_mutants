clear 
close all
project = '20190613_eveGtMut_eS1';
% project = '20190613_eveGtMut_eS1';
dataPath = ['../dat/' project '/'];
figPath = ['../fig/' project '/'];
mkdir(figPath);
% load 
load([dataPath 'nucleus_struct.mat'])
load([dataPath 'soft_fit_struct.mat'])


% extract useful vectors
ap_vec = [soft_fit_struct.ap_cell{:}]*100;
time_vec = [soft_fit_struct.time_cell{:}]/60;
fluo_vec = [soft_fit_struct.fluo_cell{:}];
hmm_index = NaN(size(fluo_vec));
hmm_sub_index = NaN(size(fluo_vec));
iter = 1;
for i = 1:numel(soft_fit_struct.fluo_cell)
    fluo = soft_fit_struct.fluo_cell{i};
    hmm_index(iter:iter+numel(fluo)-1) = i;
    hmm_sub_index(iter:iter+numel(fluo)-1) = 1:numel(fluo);
    iter = iter + numel(fluo);
end
% ref grids
ap_grid = 0:1:100;
t_inc = 1;
time_grid = 0:t_inc:50;
hmm_particle_index = soft_fit_struct.particle_index;
% tr_particle_index = [nucleus_struct.ParticleID];
% initialize arrays
occ_grid = NaN(numel(time_grid),numel(ap_grid));
kon_grid = NaN(numel(time_grid),numel(ap_grid));
koff_grid = NaN(numel(time_grid),numel(ap_grid));
% iterate through spatiotemporal regions
for a = 1:numel(ap_grid)
    for t = 1:numel(time_grid)
        ap = ap_grid(a);
        time = time_grid(t);
        tr_filter = ceil(time_vec/t_inc)*t_inc == time & ceil(ap_vec) == ap;
        index_list = hmm_index(tr_filter);
        index_u = unique(index_list);
        sub_index_list = hmm_sub_index(tr_filter);
        % iterate
        iter = 1;
        active_steps = 0;
        total_steps = 0;
        a_mat = zeros(3,3);
        if numel(index_list)>5
            for i = index_list            
                sub_indices = sub_index_list(index_list==i);
                p_zz = exp(soft_fit_struct.p_zz_log_soft{i}(:,:,sub_indices(1:end-1)));
                p_z = exp(soft_fit_struct.p_z_log_soft{i}(:,sub_indices));
               
                % occupancy first                
                active_steps = active_steps + sum(sum(p_z(2:3,:)));
                total_steps = total_steps + sum(p_z(:));
                % rates
                a_raw = sum(p_zz,3);
                a_mat = a_mat + a_raw;
            end     
            occ_grid(t,a) = active_steps / total_steps;
            a2 = zeros(2,2);
            a2(1,1) = a_mat(1,1);
            a2(2,2) = sum(sum(a_mat(2:3,2:3)));
            a2(1,2) = sum(a_mat(1,2:3));
            a2(2,1) = sum(a_mat(2:3,1));
            a2 = a2 ./ sum(a2);
            if true%~any(isnan(a2(:)))
%                 k2 = logm(a2) / 20;
                kon_grid(t,a) = a2(2,1);
                koff_grid(t,a) = a2(1,2);
            end
        end
    end
end


%%
% make figure
occ_fig = figure;
p = pcolor(flipud(imgaussfilt(occ_grid,1)));
p.EdgeAlpha = .2;
ylim([11,50])
xlim([23,58])
set(gca,'ytick',0:5:60,'yticklabels',fliplr(0:5:60))
set(gca,'xtick',0:5:60,'xticklabels',0:5:60)
ylabel('minutes into nc14')
xlabel('AP position')
h = colorbar;
ylabel(h,'fractional occupancy')
colormap(flipud(brewermap(128,'RdYlBu'))) 
set(gca,'Fontsize',12);
% caxis([0 2e5])
saveas(occ_fig,[figPath 'occupancy_heatmap.tif'])

% make figure
kon_fig = figure;
p = pcolor(flipud(imgaussfilt(kon_grid,1)));
p.EdgeAlpha = .2;
ylim([11,50])
xlim([23,58])
set(gca,'ytick',0:5:60,'yticklabels',fliplr(0:5:60))
set(gca,'xtick',0:5:60,'xticklabels',0:5:60)
ylabel('minutes into nc14')
xlabel('AP position')
h = colorbar;
ylabel(h,'activation probability')
colormap(flipud(brewermap(128,'RdYlBu'))) 
set(gca,'Fontsize',12);
% caxis([0 2e5])
saveas(kon_fig,[figPath 'kon_heatmap.tif'])

% make figure
koff_fig = figure;
p = pcolor(flipud(imgaussfilt(koff_grid,1)));
p.EdgeAlpha = .2;
ylim([11,50])
xlim([23,58])
set(gca,'ytick',0:5:60,'yticklabels',fliplr(0:5:60))
set(gca,'xtick',0:5:60,'xticklabels',0:5:60)
ylabel('minutes into nc14')
xlabel('AP position')
h = colorbar;
ylabel(h,'de-activation probability')
colormap(flipud(brewermap(128,'RdYlBu'))) 
set(gca,'Fontsize',12);
% caxis([0 2e5])
saveas(koff_fig,[figPath 'koff_heatmap.tif'])