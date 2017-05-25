run_name = 'p1d1w2f6s1r1_20170108_tau_t';
%{
# Run Simulation file
- Project 1:
    - Propofol PAC investigation
- Direction 1:
    - Thalamus-only propofol PAC modeling
- Writing item 2:
    - Journal article on results of p1d1
- Figure 6:
    - Extreme Tau_T and propofol multiplier (spm) comparisons to look at
    network frequency.
- Subfigure 1:
    - Plot of network frequency across tau_t changes
- Run 1:
    - There should probably only be a single run!
- Date Created:
    - 20170108
- Inherits from:
    - 'p1d1q43c1i1_20160919_cOFF_f3_TcurrTau.m'
%}

% Define equations of cell model (same for all populations)
eqns={
  'dV/dt=Iapp+@current';
};

time_end = 8000; % in milliseconds
numcells = 50;

g_PYsyn = 0.05;

% Create DynaSim specification structure
s=[];
s.populations(1).name='TC';
s.populations(1).size=numcells;
s.populations(1).equations=eqns;
% this is where the magic happens
s.populations(1).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                 'iLeakChing2010TC','iKLeakChing2010TC',...
                                 'CaBufferChing2010TC','iTChing2010TCTauParam',...
                                 'iHChing2010TCSwapped'};
s.populations(1).parameters={'Iapp',0};
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
s.connections(2).mechanism_list={'iGABAAChing2010','iGABABChing2010'};
s.connections(2).parameters={'gGABAA_base',0.069,'spm',3,'tauGABAA_base',5,...
                             'gGABAB',0.001};
% s.connections(3).direction='TC->TC';
% s.connections(3).mechanism_list={'iPoissonSpiketrainUncorr'};
% s.connections(3).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,...
%                              'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500};
s.connections(3).direction='RE->RE';
s.connections(3).mechanism_list={'iGABAAChing2010'};
s.connections(3).parameters={'gGABAA_base',0.069,'spm',3,'tauGABAA_base',5};
% s.connections(4).direction='RE->RE';
% s.connections(4).mechanism_list={'iPoissonSpiketrainUncorr','iGABAAChing2010'};
% s.connections(4).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,...
%                              'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500,...
%                              'gGABAA_base',0.069,'spm',1,'tauGABAA_base',5};

%   '(RE->RE,RE->TC)',   'gGABAA_base',      [3*0.069];
%   '(RE->RE,RE->TC)',   'tauGABAA_base',    [3*5];
vary={
  'TC',                'gH',       [0.0032];
  '(TC,RE)',           'Iapp',     [0.3];
  'TC',                'tauHproportion',[0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 5, 10];
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
cluster_flag = 1;
overwrite_flag = 1;
save_data_flag = 0;
% Even if `save_data_flag` is 0, if running on cluster this must be off too in
%   order to not save data?
save_results_flag = 0;
verbose_flag = 1;
compile_flag = 0;
disk_flag = 0;
downsample_factor = 6;

% local run of the simulation,
%   i.e. in the interactive session you're running this same script in
SimulateModel(s,'save_data_flag',save_data_flag,'study_dir',data_dir,...
              'cluster_flag',cluster_flag,'verbose_flag',verbose_flag,...
              'overwrite_flag',overwrite_flag,'tspan',[0 time_end],...
              'save_results_flag',save_results_flag,'solver','euler',...
              'memlimit',memlimit,'compile_flag',compile_flag,...
              'disk_flag',disk_flag,'downsample_factor',downsample_factor,...
              'vary',vary,...
              'plot_functions',{@PlotData,@PlotData},...
              'plot_options',{{'plot_type','waveform','format','png'},...
                              {'plot_type','power','format','png',...
                               'xlim',[0 40]}});
exit
