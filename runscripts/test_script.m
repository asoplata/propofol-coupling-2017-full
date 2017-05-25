%{
This is a starter script to show both how to use the DynaSim mechanisms I've
made for simulating the (Ching et al. 2010) model in Dynasim, available here
(https://github.com/asoplata/ching2010_tcre_dynasim_mechanisms) while also
giving you a good skeleton to start experimenting yourself.

The easiest way to get started with this is
1. `git clone` DynaSim to where you want to put the code
2. Add `addpath(genpath('/your/path/to/dynasim'))` to your
     `~/Documents/MATLAB/startup.m` file. Create the necessary folders and files
     if need be.
3. `git clone` this repo into '/your/path/to/dynasim/models', i.e. the 'models'
     subdirectory of your copy of the DynaSim repo.
4. Believe it or not...that should be it! You should be able to start MATLAB in
     any directory and run this script successfully! Let me know if there are
     problems, at austin.soplata 'at symbol' gmail 'dot' com

This is descended from DynaSim's `tutorial.m` script, at
https://github.com/DynaSim/DynaSim/blob/master/demos/tutorial.m
%}

% Define equations of cell model (same for all populations)
eqns={
  'dV/dt=Iapp+@current';
};

time_end = 200; % in milliseconds
numcells = 10;

% You can program parameters just like you can anything else in MATLAB code
g_PYsyn = 0.7;

% Create DynaSim specification structure
s=[];
s.populations(1).name='TC';
s.populations(1).size=numcells;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                 'CaBufferChing2010TC','iTChing2010TC',...
                                 'iHChing2010TC','iLeakChing2010TC',...
                                 'iKLeakChing2010TC'};
s.populations(1).parameters={'Iapp',0,'gH',0.001};
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
s.connections(2).parameters={'gGABAA_base',0.069,'spm',1,'tauGABAA_base',5,'gGABAB',0.001};
s.connections(3).direction='TC->TC';
s.connections(3).mechanism_list={'iPoissonSpiketrainCorr'};
s.connections(3).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500};
s.connections(4).direction='RE->RE';
s.connections(4).mechanism_list={'iPoissonSpiketrainUncorr','iGABAAChing2010'};
s.connections(4).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500,'gGABAA_base',0.069,'spm',1,'tauGABAA_base',5};

% "Vary" parameters, akw parameters to be varied -- run a simulation for all combinations of values
vary={
  'TC',              'gH',        [0.001];
  '(TC,RE)',         'Iapp',      [0, 0.2];
  '(RE->RE,RE->TC)', 'spm',       [1];
  '(TC->TC,RE->RE)', 'g_esyn',    [0.05];
  '(TC->TC,RE->RE)', 'prob_cxn',  [0.25];
  };


% Set simulation parameters
memlimit = '8G'; % unnecessary if running locally
% Obviously, set your own data directory
data_dir = '/projectnb/crc-nak/asoplata/dynasim_data/base_dynasim_run';

% % Local run of the simulation, i.e. in the interactive session you're running
% % this same script in
% SimulateModel(s,'save_data_flag',1,'study_dir',data_dir,...
%               'vary',vary,'cluster_flag',0,'verbose_flag',1,...
%               'overwrite_flag',1,'tspan',[0 time_end],'solver','euler')

% Uncomment here to instead send this simulation to the cluster, and get some
%   automated plotting on the results.
SimulateModel(s,'save_data_flag',1,'study_dir',data_dir,...
              'vary',vary,'cluster_flag',0,'verbose_flag',1,...
              'tspan',[0 time_end],'solver','euler','memory_limit',memlimit,...
              'plot_functions',{@PlotData,@PlotData,@PlotData},...
              'plot_options',{{'plot_type','waveform','format','png'},...
                              {'plot_type','rastergram','format','png'},...
                              {'plot_type','power','format','png','xlim',[0 40]}});

% I usually do stuff from the command line instead of a MATLAB session, so I use
% this:

% exit
