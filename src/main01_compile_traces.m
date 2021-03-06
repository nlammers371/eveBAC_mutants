% main01_compile_data(project,DropboxFolder,varargin)
%
% DESCRIPTION
% Funcion to compile relevant outputs from image analysis pipeline across
% multiple experiments

%
% ARGUMENTS
% project: master ID variable (should match a tab name in the Data Status
% sheet)
% DropboxFolder: full file path to folder containing compiled imaging
% results

% OPTIONS
% dropboxFolder: Pass this option, followed by the path to data folder 
%                where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% first_nc: script defaults to taking only nc14 traces. If you want
%           earlier traces, pass 'first_nc', followed by desired nuclear cycle
%           number
%
% OUTPUT: nucleus_struct: compiled data set contain key nucleus and
% particle attributes

function trace_struct = main01_compile_traces(project,FolderPath,varargin)
addpath('./utilities')
% set defaults
firstNC = 14;
minDP = 15;
pctSparsity = 50;
% two_spot_flag = contains(project, '2spot');
min_time = 0*60; % take no fluorescence data prior to this point
TresInterp = 20; 
for i = 1:numel(varargin)
    if ischar(varargin{i}) && i < numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

DataPath = ['../dat/' project '/'];
mkdir(DataPath);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Set Path Specs, ID Vars %%%%%%%%%%%%%%%%%%%%%%%%

    
% make filepath
mkdir(DataPath);
% assign save names
trace_name = [DataPath 'trace_struct.mat']; % names for compid elipse struct
file_list = dir([FolderPath '*.mat']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Obtain Relevant Filepaths %%%%%%%%%%%%%%%%%%%%%%%

% Generate master structure with info on all nuclei and traces in
% constituent sets

disp('compiling data...')
trace_struct = [];
% Loop through filenames    
for i = 1:length(file_list)      
    % extract data structures
    cp_name = [FolderPath file_list(i).name];
    processed_data = load(cp_name); % processed particles       
    % Extract compiled particles structure         
    cp_particles = processed_data.CompiledParticles;    
    if iscell(cp_particles)
        cp_particles = cp_particles{1};
    end   
    
    % check for 3D fit data    
    threeD_flag = isfield(cp_particles,'Fluo3DRaw');        
    % determine whether there is AP info
    ap_flag = isfield(cp_particles, 'APpos');
    ap_norm_flag = isfield(cp_particles,'APPos_Normalized');
    % set identifier
    setID = i;            
    % pull trace, time, and frame variables
    time_raw = processed_data.ElapsedTime*60; % time vector   
    traces_raw = processed_data.AllTracesVector; % array with a column for each trace 
    % check to see if traces are stored in cell array
    if iscell(traces_raw)
        traces_raw = traces_raw{1};
    end
    %%%%% Basic data characteristics %%%%%%
    frames_raw = 1:length(time_raw); % Frame list       
    % find index for first frame
    first_frame = max([1, processed_data.(['nc' num2str(firstNC)])]); 
    % filter trace mat and time
    traces_clean = traces_raw(first_frame:end,:);
    time_clean = time_raw(first_frame:end);    
    time_clean = time_clean - min(time_clean); % Normalize to start of nc    
    frames_clean = frames_raw(first_frame:end);                 
   
    % Iterate through traces
    temp_struct = struct;
    i_pass = 1;
    for j = 1:size(traces_clean,2)  
        % Raw fluo trace
        raw_trace = traces_clean(:,j); 
        % Get nucleus ID
        temp_struct(i).Nucleus = cp_particles(j).schnitz;    
        % Find first and last expression frames, requiring that spots
        % cannot apear earlier than some fixed time (min_time)
        trace_start = find(time_clean>=min_time&~isnan(raw_trace'),1);
        trace_stop = find(~isnan(raw_trace),1,'last');        
        %Creat versions with all intervening frames present (missing frames
        %appear as NaNs)
        trace_full = raw_trace(trace_start:trace_stop)';                                   
        trace_frames_full = frames_clean(trace_start:trace_stop); 
        time_full = time_clean(trace_start:trace_stop);
        % skip particles not in specified nc range
        if sum(~isnan(trace_full)) == 0
            continue
        end   
        
        % Find intersection btw full frame range and CP frames        
        raw_pt_frames = cp_particles(j).Frame;
        cp_filter1 = ismember(trace_frames_full,raw_pt_frames);
        cp_filter2 = ismember(raw_pt_frames,trace_frames_full);
%         raw_pt_frames = raw_pt_frames(ismember(raw_pt_frames,trace_frames_full));
        % Perform qc tests    
        nDP = sum(~isnan(trace_full));
        sparsity = prctile(diff(find(~isnan(trace_full))),pctSparsity);
        qc_flag = nDP >= minDP && sparsity == 1;             
        % Identifier variable                        
        particle = cp_particles(j).OriginalParticle;                        
        ParticleID = eval([num2str(setID) '.' sprintf('%04d',particle)]);

        % Record particle identifier
        temp_struct(i_pass).ParticleID = ParticleID; 
        % record fluorescence info             
        temp_struct(i_pass).fluo = trace_full; % note that
        temp_struct(i_pass).time = time_full; % note that
        % add offset info        
        temp_struct(i_pass).fluoOffset = NaN(size(cp_filter1));
        temp_struct(i_pass).fluoOffset(cp_filter1) = cp_particles(j).Off(cp_filter2);
        % x, y, and z info                                
        temp_struct(i_pass).xPosParticle = NaN(size(cp_filter1));
        temp_struct(i_pass).yPosParticle = NaN(size(cp_filter1));
        temp_struct(i_pass).zPosParticle = NaN(size(cp_filter1)); 
        temp_struct(i_pass).Stripe = NaN(size(cp_filter1)); 
        temp_struct(i_pass).xPosParticle(cp_filter1) = cp_particles(j).xPos(cp_filter2);
        temp_struct(i_pass).yPosParticle(cp_filter1) = cp_particles(j).yPos(cp_filter2);   
%         temp_struct(i_pass).zPosParticle(cp_filter1) = cp_particles(j).zPos(cp_filter2);
        % extract bin info'
%         temp_struct(i_pass).Stripe(cp_filter1) = cp_particles(j).Stripe(cp_filter2);
%         temp_struct(i_pass).FluoBin = cp_particles(j).particleFluoBin;
%         temp_struct(i_pass).MeanFluo = cp_particles(j).AverageFluor;
%         temp_struct(i_pass).FluoBinEdges = cp_particles(j).BinEdges;
        if ap_flag
            temp_struct(i_pass).APPosParticle = NaN(size(cp_filter1));
            temp_struct(i_pass).APPosParticle(cp_filter1) = cp_particles(j).APpos(cp_filter2)*100;            
        end
        if ap_norm_flag
            temp_struct(i_pass).APPosParticleNorm = NaN(size(cp_filter1));
            temp_struct(i_pass).APPosParticleNorm(cp_filter1) = cp_particles(j).APPos_Normalized(cp_filter2)*100;
        end
        % 3D info
        if threeD_flag
            temp_struct(i_pass).xPosParticle3D = cp_particles(j).xPosGauss3D;            
            temp_struct(i_pass).yPosParticle3D = cp_particles(j).yPosGauss3D;              
            temp_struct(i_pass).zPosParticle3D = cp_particles(j).zPosGauss3D;            
            temp_struct(i_pass).fluo3D = cp_particles(j).Fluo3DRaw;
        end
        % add qc info
        temp_struct(i_pass).N = nDP;
        temp_struct(i_pass).sparsity = sparsity;        
        temp_struct(i_pass).qc_flag = qc_flag; 
        % save file path
        temp_struct(i_pass).source = cp_name;
        i_pass = i_pass + 1;
    end      
    trace_struct = [trace_struct  temp_struct];        
end

disp('interpolating data...')
% generate interpolation fields
if threeD_flag
    interp_fields = {'fluo','fluo3D'};
else
    interp_fields = {'fluo'};
end
interpGrid = 0:TresInterp:60*60;
for i = 1:numel(trace_struct)
    time_vec = trace_struct(i).time;
    fluo_vec = trace_struct(i).fluo;
    start_i = find(~isnan(fluo_vec),1);
    stop_i = find(~isnan(fluo_vec),1,'last');    
    time_vec = time_vec(start_i:stop_i);       
    if numel(time_vec)>1%nucleus_struct(i).qc_flag == 1
        time_interp = interpGrid(interpGrid>=time_vec(1)&interpGrid<=time_vec(end));
        for  j = 1:numel(interp_fields)
            vec = trace_struct(i).(interp_fields{j})(start_i:stop_i);
            %Look for clusters of 6 or more NaNs
            kernel = [1,1,1,1,1];
            vec_nans = isnan(vec);
            vec_conv = conv(kernel,vec_nans);
            vec_conv = vec_conv(3:end-2);
            z_ids = find(vec_conv==numel(kernel));
            z_ids = unique([z_ids-1 z_ids z_ids+1]); % get set of z_ids    
            vec(z_ids) = 0; % set clusters to zeros    
            vec(vec<0) = 0; % deal with negative values    
            % find single dp "blips". These will be replaced via interpolation
            % interpolate remaining NaNs    
            query_points = time_vec(isnan(vec));
            interp_t = time_vec(~isnan(vec));
            interp_f = vec(~isnan(vec)); 
            new_f = interp1(interp_t,interp_f,query_points);      
            vec(ismember(time_vec,query_points)) = new_f;        
            % Interpolate to standardize spacing               
            trace_struct(i).([interp_fields{j} '_interp']) = interp1(time_vec,vec,time_interp);
        end
    else
        time_interp = NaN;
        for  j = 1:numel(interp_fields)
            trace_struct(i).([interp_fields{j} '_interp']) = NaN;
        end
    end
    trace_struct(i).time_interp = time_interp;
    trace_struct(i).TresInterp = TresInterp;
end
% call function to calculate average psf difs
% save
save(trace_name ,'trace_struct') 

disp('done.')
