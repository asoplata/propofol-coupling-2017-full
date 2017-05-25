run_name = 'p1d1w2f7to8r2_20170109_triple_switch';
%{
# Run Simulation file
- Project 1:
    - Propofol PAC investigation
- Direction 1:
    - Thalamus-only propofol PAC modeling
- Writing item 2:
    - Journal article on results of p1d1
- Figures 7-8:
    - Illustrating and showing example of phase-amplitude coupling across
    corticothalamic UP/DOWN SWO states in action
- Subfigures:
    - N/A
- Run 2:
    - Running long simulation AND passing `CalcPAC` to run it in batch, trying
    fewer bins for smoother picture
- Date Created:
    - 20170109
- Inherits from:
    - 'p1d1q36c5i3_20161106_triple_switch_sfn2016poster.m'
%}

% Define equations of cell model (same for all populations)
eqns={
  'dV/dt=Iapp+@current';
};

first_onset =       0; % in ms

% time_end = 3000; % in milliseconds
% second_onset =   1000;
% third_onset =    2000;
% numcells = 4;

time_end = 11000; % in milliseconds
second_onset =   1000;
third_onset =    6000;
numcells = 50;


square_waves_freq = 0.6; % in Hz
% Note: the uses of 'onset' etc. overlap between mechanisms, and while in
%   general this is poor programming practice, it can act as a discrepancy-check
%   here

g_PYsyn = 0.05;
prob_cxn = 0.5;

% Create DynaSim specification structure
s=[];
s.populations(1).name='TC';
s.populations(1).size=numcells;
s.populations(1).equations=eqns;
% this is where the magic happens
s.populations(1).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                 'iLeakChing2010TC','iKLeakChing2010TC',...
                                 'CaBufferChing2010TC','iTChing2010TC',...
                                 'iHChing2010TCSwitchSwa'};
s.populations(1).parameters={   'gH1',0.01,          'gH2',0.003,       'gH3',0.0018,...
                             'onset1',first_onset,'onset2',second_onset,'onset3',third_onset,...
                             'Iapp',0};
s.populations(2).name='RE';
s.populations(2).size=numcells;
s.populations(2).equations=eqns;
s.populations(2).mechanism_list={'iNaChing2010RE','iKChing2010RE',...
                                 'iTChing2010RE','iLeakChing2010RE'};
s.populations(2).parameters={'Iapp',0};

s.connections(1).direction='TC->RE';
s.connections(1).mechanism_list={'iAMPAChing2010'};
s.connections(1).parameters={'gAMPA',0.08};
s.connections(2).direction='RE->TC';
s.connections(2).mechanism_list={'iGABAAChing2010Switch','iGABABChing2010'};
s.connections(2).parameters={  'spm1',1,            'spm2',2,             'spm3',3,...
                             'onset1',first_onset,'onset2',second_onset,'onset3',third_onset,...
                             'gGABAA_base',0.069,'tauGABAA_base',5,...
                             'gGABAB',0.001};
s.connections(3).direction='TC->TC';
s.connections(3).mechanism_list={'iPSUSW','iSquare','iTonicSwitch'};
s.connections(3).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,...
                             'tau_i',10,'prob_cxn',prob_cxn,...
                             'poisson_square_freq',square_waves_freq,...
                             'square_amp',0.5,'square_freq',square_waves_freq,...
                             'stim1',0.7,         'stim2',0.0,          'stim3',-0.3,...
                             'onset1',first_onset,'onset2',second_onset,'onset3',third_onset,...
                             };
s.connections(4).direction='RE->RE';
s.connections(4).mechanism_list={'iPSUSW','iSquare','iGABAAChing2010Switch','iTonicSwitch'};
s.connections(4).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,...
                             'tau_i',10,'prob_cxn',prob_cxn,...
                             'poisson_square_freq',square_waves_freq,...
                             'square_amp',0.5,'square_freq',square_waves_freq,...
                             'spm1',1,            'spm2',2,             'spm3',3,...
                             'stim1',0.7,         'stim2',0.2,          'stim3',-0.3,...
                             'onset1',first_onset,'onset2',second_onset,'onset3',third_onset,...
                             'gGABAA_base',0.069,'tauGABAA_base',5};


%   '(RE->RE,RE->TC)',   'gGABAA_base',      [3*0.069];
%   '(RE->RE,RE->TC)',   'tauGABAA_base',    [3*5];
%   'TC',                'gH',       [0.0032];
%   '(RE->TC,RE->RE)',   'spm',      [1,2,3,4,5,6,7,8];
vary={
  '(TC,RE)',           'Iapp',     [0.0];
};

%% Set simulation parameters
% How much RAM, options: 8G?, 24, 48, 96, 128
% are segfaults from too long of a sim?
% memlimit = '8G';
% memlimit = '16G';
% memlimit = '48G';
% memlimit = '96G';
memlimit = '254G';

% Save data/results to this directory. If just a single name, will
%   save to that directory name in the current directory from which it's run.
%   Will create directory if it does not exist.
data_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                  run_name);

% Flags
cluster_flag =      1;
overwrite_flag =    1;
save_data_flag =    1;
% Even if `save_data_flag` is 0, if running on cluster this must be off too in
%   order to not save data?
save_results_flag = 1;
verbose_flag =      1;
compile_flag =      0;
disk_flag =         0;
downsample_factor = 1;

% local run of the simulation,
%   i.e. in the interactive session you're running this same script in
SimulateModel(s,'save_data_flag',save_data_flag,'study_dir',data_dir,...
              'cluster_flag',cluster_flag,'verbose_flag',verbose_flag,...
              'overwrite_flag',overwrite_flag,'tspan',[0 time_end],...
              'save_results_flag',save_results_flag,'solver','euler',...
              'memlimit',memlimit,'compile_flag',compile_flag,...
              'disk_flag',disk_flag,'downsample_factor',downsample_factor,...
              'vary',vary,...
              'plot_functions',{@PlotData,@PlotData,@CalcPAC},...
              'plot_options',{{'plot_type','waveform','format','png'},...
                              {'plot_type','power','format','png','xlim',[0 40]},...
                              {'n_bins',10}});
exit
