function all_dir = add_path_create_dir(homepath)
% Add code path and creat data path
%% Add code path to the environment
folder = fileparts(which(mfilename));
addpath(genpath(folder));

%% Mac or Windows
if ispc
    mydash = '\';
else
    mydash = '/';
end

%% Creat path
if nargin == 0
    all_dir.save_dir = [folder, mydash, 'Learn_Kernel_output',mydash];
else
    all_dir.save_dir = [homepath, mydash, 'Learn_Kernel_output',mydash];
end

all_dir.data_dir	= [all_dir.save_dir,'data',mydash]; 
all_dir.infer_dir	= [all_dir.save_dir,'inference',mydash]; 
all_dir.fig_dir     = [all_dir.save_dir,'figure',mydash]; 


if ~exist(all_dir.save_dir,'dir');           mkdir(all_dir.save_dir);   end
if ~exist(all_dir.data_dir,'dir');           mkdir(all_dir.data_dir);   end
if ~exist(all_dir.infer_dir,'dir');          mkdir(all_dir.infer_dir);  end
if ~exist(all_dir.fig_dir,'dir');            mkdir(all_dir.fig_dir);    end


addpath(all_dir.save_dir);
addpath(all_dir.data_dir);
addpath(all_dir.infer_dir);
addpath(all_dir.fig_dir);

% %% Save data to a local DIR, OUTSIDE of the repo: (in your user home DIR)
% SAVE_DIR = [SAVE_DIR, mydash,'QuanjunFei',mydash,'MF',mydash];
% if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end
% addpath(SAVE_DIR)
%
%
% dataDIR                 = [SAVE_DIR,'data_output',      mydash];     % save data here; will be read by code
% datafigDIR              = [dataDIR, 'figures',  mydash];           % data fig DIR, prevously figDIR_solu
% inferenceDIR            = [SAVE_DIR,'inference_output', mydash];
% inference_PDE_DIR       = [inferenceDIR, 'PDE_error',  mydash];
% inference_PS_DIR        = [inferenceDIR, 'PS_error',  mydash];
%

%
% addpath(dataDIR);
% addpath(inferenceDIR);

end