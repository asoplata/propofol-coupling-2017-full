%{
# p1d1q40c2i1_20160727_htypo_cortOFF_large
- Question 40:
  redoing main gh-pol planes for p1d1w2 (a la p1d1q35c4) , but with H-current
  typo fixed
- Computer model run 2:
  - using p1d1q35c4i1_fig4_cortOFF as template
  - doing sizes of sims of what I THOUGHT I was doing: 50 of each cell type for
    8 seconds
- Iteration 1:
  - cort OFF
%}

% Save data/results/whatever to which directory. If just a single name, will
%   save to that directory name in the current directory from which it's run.
data_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                  'p1d1q40c2i1_20160727_htypo_cortOFF_large');

% Define equations of cell model (same for all populations)
eqns={
  'dV/dt=Iapp+@current';
};

time_end = 6000; % in milliseconds
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
                                 'CaBufferChing2010TC','iTChing2010TC',...
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
s.connections(2).parameters={'gGABAA_base',0.069,'spm',1,'tauGABAA_base',5,...
                             'gGABAB',0.001};
% s.connections(3).direction='TC->TC';
% s.connections(3).mechanism_list={'iPoissonSpiketrainUncorr'};
% s.connections(3).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,...
%                              'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500};
s.connections(3).direction='RE->RE';
s.connections(3).mechanism_list={'iGABAAChing2010'};
s.connections(3).parameters={'gGABAA_base',0.069,'spm',1,'tauGABAA_base',5};
% s.connections(4).direction='RE->RE';
% s.connections(4).mechanism_list={'iPoissonSpiketrainUncorr','iGABAAChing2010'};
% s.connections(4).parameters={'g_esyn',g_PYsyn,'rate',12,'T',time_end,...
%                              'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500,...
%                              'gGABAA_base',0.069,'spm',1,'tauGABAA_base',5};

vary={
  'TC',                'gH',       [0.0003, 0.0004, 0.0006, 0.0007, 0.0010, 0.0013, 0.0018, 0.0024, 0.0032, 0.0042, 0.0056, 0.0075, 0.01, 0.0133, 0.0178, 0.0237, 0.0316, 0.0422, 0.0562, 0.0750, 0.1, 0.1334, 0.1778, 0.2371, 0.3162];
  '(TC,RE)',           'Iapp',     [-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5];
  '(RE->RE,RE->TC)',   'spm',      [1,2,3];
};

% % Example with reliable spikes/oscs
% vary={
%   'TC',                'gH',       [0.002, 0.004];
%   '(TC,RE)',           'Iapp',     [0.0,0.1]
%   '(RE->RE,RE->TC)',   'spm',      [3];
% };

%% Set simulation parameters

% How much RAM, options: 8G?, 24, 48, 96, 128
% are segfaults from too long of a sim?
% memlimit = '8G';
% memlimit = '16G';
% memlimit = '48G';
% memlimit = '96G';
memlimit = '254G';

% Flags
cluster_flag = 1;
overwrite_flag = 1;
save_data_flag = 0;
% Even if `save_data_flag` is 0, if running on cluster this must be off too in
%   order to not save data?
save_results_flag = 1;
verbose_flag = 1;
compile_flag = 0;
disk_flag = 0;

% local run of the simulation,
%   i.e. in the interactive session you're running this same script in
SimulateModel(s,'save_data_flag',save_data_flag,'study_dir',data_dir,...
              'cluster_flag',cluster_flag,'verbose_flag',verbose_flag,...
              'overwrite_flag',overwrite_flag,'tspan',[0 time_end],...
              'save_results_flag',save_results_flag,...
              'solver','euler','memlimit',memlimit,'compile_flag',compile_flag,...
              'disk_flag',disk_flag,...
              'vary',vary,...
              'plot_functions',{@PlotData,@PlotData},...
              'plot_options',{{'plot_type','waveform','format','png'},...
                              {'plot_type','power','format','png',...
                               'xlim',[0 40]}});
exit

% for ii=1:2
%   figure(ii)
%   subplot 811
%   plot(data(ii).time, data(ii).TC_V)
%   title('V')
%   subplot 812
%   plot(data(ii).time, data(ii).TC_iTChing2010TC_ITChing2010TC)
%   title('I T')
%   subplot 813
%   plot(data(ii).time, data(ii).TC_CaBufferChing2010TC_CaBuffer)
%   title('Ca')
%   subplot 814
%   plot(data(ii).time, data(ii).TC_iNaChing2010TC_INaChing2010TC)
%   title('I Na')
%   subplot 815
%   plot(data(ii).time, data(ii).TC_iNaChing2010TC_mNa)
%   hold on
%   plot(data(ii).time, data(ii).TC_iNaChing2010TC_mNa.^3,'r')
%   hold off
%   title('mNa and mNa^3')
%   subplot 816
%   plot(data(ii).time, data(ii).TC_iKChing2010TC_IKChing2010TC)
%   title('I K')
%   subplot 817
%   plot(data(ii).time, data(ii).TC_iLeakChing2010TC_ILeakChing2010TC)
%   title('I Leak')
%   subplot 818
%   plot(data(ii).time, data(ii).TC_iKLeakChing2010TC_IKLeakChing2010TC)
%   title('I K Leak')
% end
% 
% 
% 
% figure(3)
% X = -100:0.01:60;
% v_shift = 35;
% plot(X,((.32.*(13-(X+v_shift)))./(exp((13-(X+v_shift))./4)-1)))
% title('mNa activation')
% 
% for ii=1:2
%   figure(3+ii)
% 
%   subplot 211
%   plot(data(ii).time, data(ii).TC_V)
%   title('V')
% 
%   subplot 212
%   plot(data(ii).time, data(ii).TC_iNaChing2010TC_mNa)
%   hold on
%   plot(data(ii).time, data(ii).TC_iNaChing2010TC_mNa.^3,'r')
%   hold off
%   title('mNa and mNa^3')
% end
% 

% figure(2)
% PlotData(data,'variable','V','plot_type','waveform')
% 
% figure(3)
% power_data = CalcPower(data)
% PlotData(power_data,'variable','V','plot_type','power')
