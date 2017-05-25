run_name = 'p1d1w2f10_20170104_gen_py_spiketrain';
%{
This script should do all the MATLAB work required to generate the cortical
component of Figure 2 for paper "p1d1j1-pac-paper". Note that this script is not
sufficient to generate ALL the figures or even all parts of some figures, as
some parts are done via Inkscape, etc.
%}

% Initialize workspace & clear memory/everything
clear
clf
close all

% Set all figure backgrounds to white
set(gcf,'color','w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 2, model diagram (cortical component)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Figure 2, let's use some actual simulated cortical spiking data for the
%   its representation within the model diagram, and let's also generate a plot
%   for the simple step function of tonic charge.

%% Simulated cortical spiking data
% Let's simulate a tiny 2 cell model to get two exemplar spike trains

% Define equations of cell model (same for all populations)
eqns={
  'dV/dt=Iapp+@current';
};

time_end = 2000; % in milliseconds
numcells = 2;

% Create DynaSim specification structure
s=[];
s.populations(1).name='TC';
s.populations(1).size=numcells;
s.populations(1).equations=eqns;
s.populations(1).mechanism_list={'iNaChing2010TC','iKChing2010TC',...
                                 'iLeakChing2010TC'};
s.populations(1).parameters={'Iapp',0,'gH',0.001};
s.connections(1).direction='TC->TC';
% Note that the Poisson spiketrain needs the total time input too, or else
%   you'll get a matrix-size compatibility error
s.connections(1).mechanism_list={'iPoissonSpiketrainUncorr'};
s.connections(1).parameters={'g_esyn',0.05,'rate',12,'T',time_end,'tau_i',10,'prob_cxn',0.5,'jitter_stddev',500};

% Save data.
data_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                  run_name);

% Local run of the simulation,
%   i.e. in the interactive session you're running this same script in
data.fig10 = SimulateModel(s,'save_data_flag',1,'study_dir',data_dir,...
                             'cluster_flag',0,'verbose_flag',1,...
                             'overwrite_flag',1,'tspan',[0 time_end],'solver','euler');

% Set our detection threshold for what counts as a "spike" in the units of
%   unitless quantities of Poisson Spiketrain 'activation'
threshold = 2;

% Note: this is from my personal library, not DynaSim
spike_times = timeSeriesToSpikeTimes(data.fig10.model.fixed_variables.TC_TC_iPoissonSpiketrainUncorr_Ge',...
                                     data.fig10.time, threshold);

figure(1)
% Plot the rasters for the spiketrain
subplot(1,2,1)
hold on
% Note: this is third-party plotting software downloaded from MathWorks Central
[raster_x, raster_y] = plotSpikeRaster(spike_times,'PlotType','vertline',...
                       'VertSpikeHeight',0.67);
text(0.29*time_end,0.5,'Artificial PY spikes','fontsize',14)
axis off
ylim([0,5])
hold off

% Plot a step function to illustrate the addition of tonic
% depolarization/hyperpolarization, aka applied current
subplot(1,2,2)
step_function = zeros(1000,1);
step_function(500:1000) = 1;
plot(step_function,'k')
ylim([-1,1.5])
axis off
text(75,1.25,'Applied Current (in $$\pm \frac{\mu A}{cm^2})$$','interpreter','latex','fontsize',14,'fontweight','bold')

saveas(gcf,'fig10_cortex_input_generated.png')
