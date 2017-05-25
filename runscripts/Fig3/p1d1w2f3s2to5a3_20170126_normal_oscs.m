run_name = 'p1d1w2f3s2to5a3_20170126_normal_oscs';
%{
# Analyse Simulation file
- DEPENDENCIES:
    - This requires `p1d1q42c1i1...` through `p1d1q42c1i4...` simulation files
    to have run, and have their data accessible
    - Run this TWICE in a live session to get consistent plotting, because
    MATLAB is dumb.
- Project 1:
    - Propofol PAC investigation
- Direction 1:
    - Thalamus-only propofol PAC modeling
- Writing item 2:
    - Journal article on results of p1d1
- Figure 3:
    - Baseline gH-background excitation plane along with regular and zoomed
    voltage traces of each of regular 4 behavior patterns of the hyperpolarized
    thalamus, in addition to the T-current window in/activation curves. The 4
    stereotypical identity patterns come from sim # 595 silent hyperpolarization, 607
    sub-alpha oscillation, 1153 spindling/transients, and 1177 silent from
    p1d1q40c2i1
- Subfigures 2-5:
    - Everything except the gH-background excitation plane
- Analysis 3:
    - Plotting just the reproduced TC #1 voltage traces, drawing a box over the
    TC T-current activation window, and showing a zoom around the window for each
    simulation. Now removing all ticklabels. Now restoring ticklabels to make
    FIG files for michelle
- Date Created:
    - 20170106
- Inherits from:
    - `p1d1q42c1i1to4a1_20160919_cOFF_f3_repro_all_ans.m`
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Figure specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all figure backgrounds to white
set(gcf,'color','w');

le_fontsize = 22;
le_linewidth = 1.5;
cornflowerblue = [0.4 0.6 0.92]; % because fight club
tc_neuron = 7;

outgoing_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                      run_name);
if exist(outgoing_dir,'dir') == 0
    mkdir(outgoing_dir)
end
outgoing_plots_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                            run_name, '/plots/');
if exist(outgoing_plots_dir,'dir') == 0
    mkdir(outgoing_plots_dir)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Plot TC T-current window
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Let's plot the T-current's activation curves

figno = 1;
voltage = [-95:0.01:0];

% Judging from visual inspection, we can safely say the EFFECTIVE (as in, very
% high probability of leading to burst) T-current window of
%   activation is limited to
tcurr_window.ceiling = -74;
tcurr_window.floor = -80;

tcurr_Minf = 1./(1+exp((-((voltage+2)+57))./6.2));
% Must be squared because that's what the T-current Current equation sees
p1 = plot((tcurr_Minf.^2), voltage, 'b',...
     'LineWidth',le_linewidth)

% Now let's add the T-current in-activation curve! (better name: lack of
%   inactivation curve...grrr....)
tcurr_Hinf = 1./(1+exp((((voltage+2)+81))./4));
hold on
p2 = plot((tcurr_Hinf), voltage, 'r',...
     'LineWidth',le_linewidth)
hold off

legend([p1 p2],{'I_T activation m_T','I_T inactivation h_T'},'Location','northeast')
% title('I_T activation and inactivation curves')
% xlabel('Activation (unitless)')
% ylabel('Voltage (mV)')
xlim([0 0.5])
ylim([-85 -55])

hold on
u = [0;          0;            0.5;            0.5];
v = [tcurr_window.floor; tcurr_window.ceiling; tcurr_window.ceiling; tcurr_window.floor];
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
% Oldest version that does this transparency is 2014b...really. ugh MATLAB
set(p,'FaceAlpha',0.3)
hold off

set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
                   '_Tcurr_activ.fig'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Our 'vary' use previously was
%   - sim #595  silent hyperpolarization
%       vary={
%         'TC',                'gH',       [0.0018];
%         '(TC,RE)',           'Iapp',     [-0.3];
%   - sim #607  sub-alpha oscillation
%     'TC',                'gH',       [0.0018];
%     '(TC,RE)',           'Iapp',     [0.1];
%   - sim #1153 spindling/transients
%     'TC',                'gH',       [0.01];
%     '(TC,RE)',           'Iapp',     [-0.3];
%   - sim #1177 silent hyperpolarization
%     'TC',                'gH',       [0.01];
%     '(TC,RE)',           'Iapp',     [0.4];

if ~exist('shpol','var')
  shpol_dir =    'p1d1q42c1i1_20160919_cOFF_f3_repro_silenthpol';
  subalpha_dir = 'p1d1q42c1i2_20160919_cOFF_f3_repro_subalpha';
  spindling_dir ='p1d1q42c1i3_20160919_cOFF_f3_repro_spindling';
  sdpol_dir =    'p1d1q42c1i4_20160919_cOFF_f3_repro_silentdpol';

  shpol = load(strcat(strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                             shpol_dir),...
                      '/data/', 'study_sim1_data.mat'));
  subalpha = load(strcat(strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                             subalpha_dir),...
                      '/data/', 'study_sim1_data.mat'));
  spindling = load(strcat(strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                             spindling_dir),...
                      '/data/', 'study_sim1_data.mat'));
  sdpol = load(strcat(strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                             sdpol_dir),...
                      '/data/', 'study_sim1_data.mat'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Plot relevant/interesting raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01;
downsample_factor = 6;
% Beginning of timespan we're interested in showing, in ms (after the initial
%   condition transient)
timebegin = 1;
% End of timespan we're interested in showing, in ms
timelimit = 8000;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];

number_bins = 1000;

figno = figno + 1;
figure(figno)
plot(shpol(1).time(timespan_indices),...
     shpol(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
% Add a transparent plane to illustrate the zoom
u = [timebegin;          timebegin;            timelimit;            timelimit];
v = [tcurr_window.floor; tcurr_window.ceiling; tcurr_window.ceiling; tcurr_window.floor];
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
% Oldest version that does this transparency is 2014b...really. ugh MATLAB
set(p,'FaceAlpha',0.3)
% The only way to be REALLY SURE that it will actually increase the fontsize is
%   to do it separately for EVERY axes object, AFTER everything's already been
%   plotted
%   I.e. Must change FontSize AFTER `plot()` has been called, since `plot()`
%   "resets" the axes object...seriously I HATE MATLAB
%   https://www.mathworks.com/matlabcentral/answers/131236-changing-fontsize-in-a-figure
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 60])
xlim([0 8000])
% title('TC#1 Voltage, silent hyperpolarization')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_vtrace_shpol.fig'))


figno = figno + 1;
figure(figno)
plot(shpol(1).time(timespan_indices),...
     shpol(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
set(p,'FaceAlpha',0.3)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 -60])
xlim([0 8000])
% title('Zoomed TC#1 Voltage, silent hyperpolarization')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoom_vtrace_shpol.fig'))

figno = figno + 1;
figure(figno)
plot(subalpha(1).time(timespan_indices),...
     subalpha(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
% Add a transparent plane to illustrate the zoom
u = [timebegin;          timebegin;            timelimit;            timelimit];
v = [tcurr_window.floor; tcurr_window.ceiling; tcurr_window.ceiling; tcurr_window.floor];
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
% Oldest version that does this transparency is 2014b...really. ugh MATLAB
set(p,'FaceAlpha',0.3)
% The only way to be REALLY SURE that it will actually increase the fontsize is
%   to do it separately for EVERY axes object, AFTER everything's already been
%   plotted
%   I.e. Must change FontSize AFTER `plot()` has been called, since `plot()`
%   "resets" the axes object...seriously I HATE MATLAB
%   https://www.mathworks.com/matlabcentral/answers/131236-changing-fontsize-in-a-figure

set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 60])
xlim([0 8000])
% title('TC#1 Voltage, sub-alpha oscillation')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_vtrace_subalpha.fig'))


figno = figno + 1;
figure(figno)
plot(subalpha(1).time(timespan_indices),...
     subalpha(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
set(p,'FaceAlpha',0.3)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 -60])
xlim([0 8000])
% title('Zoomed TC#1 Voltage, sub-alpha oscillation')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoom_vtrace_subalpha.fig'))


figno = figno + 1;
figure(figno)
plot(spindling(1).time(timespan_indices),...
     spindling(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
% Add a transparent plane to illustrate the zoom
u = [timebegin;          timebegin;            timelimit;            timelimit];
v = [tcurr_window.floor; tcurr_window.ceiling; tcurr_window.ceiling; tcurr_window.floor];
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
% Oldest version that does this transparency is 2014b...really. ugh MATLAB
set(p,'FaceAlpha',0.3)
% The only way to be REALLY SURE that it will actually increase the fontsize is
%   to do it separately for EVERY axes object, AFTER everything's already been
%   plotted
%   I.e. Must change FontSize AFTER `plot()` has been called, since `plot()`
%   "resets" the axes object...seriously I HATE MATLAB
%   https://www.mathworks.com/matlabcentral/answers/131236-changing-fontsize-in-a-figure

set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 60])
xlim([0 8000])
% title('TC#1 Voltage, spindling')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_vtrace_spindling.fig'))


figno = figno + 1;
figure(figno)
plot(spindling(1).time(timespan_indices),...
     spindling(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
set(p,'FaceAlpha',0.3)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 -60])
xlim([0 8000])
% title('Zoomed TC#1 Voltage, spindling')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoom_vtrace_spindling.fig'))


figno = figno + 1;
figure(figno)
plot(sdpol(1).time(timespan_indices),...
     sdpol(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
% Add a transparent plane to illustrate the zoom
u = [timebegin;          timebegin;            timelimit;            timelimit];
v = [tcurr_window.floor; tcurr_window.ceiling; tcurr_window.ceiling; tcurr_window.floor];
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
% Oldest version that does this transparency is 2014b...really. ugh MATLAB
set(p,'FaceAlpha',0.3)
% The only way to be REALLY SURE that it will actually increase the fontsize is
%   to do it separately for EVERY axes object, AFTER everything's already been
%   plotted
%   I.e. Must change FontSize AFTER `plot()` has been called, since `plot()`
%   "resets" the axes object...seriously I HATE MATLAB
%   https://www.mathworks.com/matlabcentral/answers/131236-changing-fontsize-in-a-figure

set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 60])
xlim([0 8000])
% title('TC#1 Voltage, silent depolarization')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_vtrace_sdpol.fig'))


figno = figno + 1;
figure(figno)
plot(sdpol(1).time(timespan_indices),...
     sdpol(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
hold on
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
set(p,'FaceAlpha',0.3)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-95 -60])
xlim([0 8000])
% title('Zoomed TC#1 Voltage, silent depolarization')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoom_vtrace_sdpol.fig'))
