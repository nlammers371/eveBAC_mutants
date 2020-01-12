clear 
close all

% set path to inference results
project = 'eveGtS2-NullS1';
ResultsPath = ['../out/' project '/'];
FileList = dir([ResultsPath '/*.mat']);
% make figure path
FigPath = ['../fig/' project '/'];
mkdir(FigPath);

% read files
inference_results = struct;
for f = 1:numel(FileList)
    load([ResultsPath FileList(f).name])
    fnames = fieldnames(output);
    if numel(fnames) > 5
        for fn = 1:numel(fnames) 
            inference_results(f).(fnames{fn}) = output.(fnames{fn});
        end
    end
end

