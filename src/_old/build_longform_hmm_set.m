clear 
close all
project = '20190613_eveGtMut_eS1';
% project = '20190613_eveGtMut_eS1';
dataPath = ['../dat/' project '/'];
writePath = ['../out/' project '/'];
mkdir(writePath)
figPath = ['../fig/' project '/'];
mkdir(figPath);
% load 
load([dataPath 'nucleus_struct.mat'])
load([dataPath 'soft_fit_struct.mat'])
%%
% pull x and y info into soft structr
ss_particle_index = [soft_fit_struct.particle_index];
nc_particle_index = [nucleus_struct.ParticleID];
xPos_cell = cell(size(ss_particle_index));
yPos_cell = cell(size(ss_particle_index));
for i = 1:numel(ss_particle_index)
    ind = find(ss_particle_index(i) == nc_particle_index);
    nc_time =  round(nucleus_struct(ind).time_interp);
    ss_time =  round(soft_fit_struct.time_cell{i});
    xPos_cell{i} = interp1(nc_time,nucleus_struct(ind).xPos_interp,ss_time);
    yPos_cell{i} = interp1(nc_time,nucleus_struct(ind).yPos_interp,ss_time);
end
soft_fit_struct.xPos_cell = xPos_cell;
soft_fit_struct.yPos_cell = yPos_cell;

% extract vectors
ap_vec = [soft_fit_struct.ap_cell{:}]*100;
time_vec = [soft_fit_struct.time_cell{:}]/60;
x_vec = [soft_fit_struct.xPos_cell{:}];
y_vec = [soft_fit_struct.yPos_cell{:}];
fluo_vec = [soft_fit_struct.fluo_cell{:}];
particle_index = soft_fit_struct.particle_index;
%% Iterate through structure and extract promoter states
p_z_log_cell = soft_fit_struct.p_z_log_soft;
promoter_state_vec = NaN(size(ap_vec));
particle_id_vec = NaN(size(ap_vec));
iter = 1;
for i = 1:numel(p_z_log_cell)
    [~, z_vec_bin] = max(exp(p_z_log_cell{i}));
    promoter_state_vec(iter:iter+numel(z_vec_bin)-1) = z_vec_bin>1;
    particle_id_vec(iter:iter+numel(z_vec_bin)-1) = repelem(particle_index(i),numel(z_vec_bin));
    iter = iter+numel(z_vec_bin);
end

results_array = [double(particle_id_vec') double(time_vec') double(ap_vec') double(x_vec')...
    double(y_vec') double(fluo_vec') double(promoter_state_vec')];

results_table = array2table(results_array, 'VariableNames', {'particle_id',...
    'time', 'ap', 'x_pos', 'y_pos', 'fluo', 'promoter_state'});

save([writePath 'results_table.mat'],'results_table')
