run_name = 'p1d1w2f5s1to6a5_20170126_spindling_comparison_ans';
%{
# Analyse Simulation file
- DEPENDENCIES:
    - Requires 'p1d1w2f5s1to6r1_20170107_spindling_comparison.m' to have run
- Project 1:
    - Propofol PAC investigation
- Direction 1:
    - Thalamus-only propofol PAC modeling
- Writing item 2:
    - Journal article on results of p1d1
- Figure 5:
    - Explanations of propofol GABA-A sustained alpha mechanisms, using
    comparisons of baseline-silent depolarization & highdose-alpha, and
    baseline-spindling & highdose-alpha. Includes extreme Tau_T and propofol
    multiplier (spm) comparisons.
- Subfigures 1-6:
    - Just the baseline-spindling and highdose-alpha cases for comparison. Using
    slightly depolarized state 0.1, higher gh 0.0032
- Analysis 5:
    - Plotting spindle traces in new spindle dark green, on cluster NOT local
    Arch, because stupid plots axis tick labels are different sizes on local
    Arch vs cluster. I hate matlab again. Now removing all ticklabels. Now
    restoring all ticklabels to give Michelle FIG files
- Date Created:
    - 20170126
- Inherits from:
    - `p1d1q41c4i5a1_20160905_cOFF_propomech_f5_spindling_ans.m`
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Figure specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all figure backgrounds to white
set(gcf,'color','w');

le_fontsize = 22;
le_linewidth = 1.5;
cornflowerblue = [0.4 0.6 0.92]; % because fight club
spindling_green = [0 0.75 0];
tc_neuron = 20;

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
%% 1. Load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_name = 'p1d1w2f5s1to6r1_20170107_spindling_comparison';
data_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                  data_name);

if ~exist('baseline','var')
  baseline = load(strcat(data_dir, '/data/', 'study_sim1_data.mat'));
  highdose = load(strcat(data_dir, '/data/', 'study_sim2_data.mat'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Plot relevant/interesting raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01;
downsample_factor = 6;
% Beginning of timespan we're interested in showing, in ms (after the initial
%   condition transient)
timebegin = 30;
% End of timespan we're interested in showing, in ms
timelimit = 2000;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];

number_bins = 1000;
figno = 1;

% From 'p1d1q42c1i1to4a1_20160919_cOFF_f3_repro_all_ans':
%   Judging from visual inspection, we can safely say the EFFECTIVE (as in, very
%     high probability of leading to burst) T-current window of
%     activation is limited to
tcurr_window.ceiling = -74;
tcurr_window.floor = -80;

% figno = figno + 1;
figure(figno)
plot(highdose(1).time(timespan_indices),...
     highdose(1).TC_V(timespan_indices,tc_neuron),...
     'r','LineWidth',le_linewidth)
hold on
plot(baseline(1).time(timespan_indices),...
     baseline(1).TC_V(timespan_indices,tc_neuron),...
     'Color',spindling_green,'LineWidth',1.5*le_linewidth)
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
% ylabel('Voltage in mV')
ylim([-90 60])
% title('TC#1 Voltage traces')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_vtrace.fig'))

figno = figno + 1;
figure(figno)
plot(baseline(1).time(timespan_indices),...
     baseline(1).TC_V(timespan_indices,tc_neuron),...
     'Color',spindling_green,'LineWidth',1.5*le_linewidth)
hold on
plot(highdose(1).time(timespan_indices),...
     highdose(1).TC_V(timespan_indices,tc_neuron),...
     'r','LineWidth',le_linewidth)
p = patch(u,v,'b','Marker','.','LineWidth',0.00001);
set(p,'FaceColor',cornflowerblue)
set(p,'FaceAlpha',0.3)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([-80 -60])
% title('Zoomed TC#1 Voltages')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoom_vtrace.fig'))

figno = figno + 1;
figure(figno)
plot(baseline(1).time(timespan_indices),...
     -baseline(1).TC_iHChing2010TCSwapped_IHChing2010TCSwapped(timespan_indices,tc_neuron),...
     'Color',spindling_green,'LineWidth',1.5*le_linewidth)
hold on
plot(highdose(1).time(timespan_indices),...
     -highdose(1).TC_iHChing2010TCSwapped_IHChing2010TCSwapped(timespan_indices,tc_neuron),...
     'r','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([0.1 0.2])
% title('TC#1 H-current -I_H')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_I_H.fig'))

figno = figno + 1;
figure(figno)
plot(baseline(1).time(timespan_indices),...
     baseline(1).TC_iTChing2010TC_hT(timespan_indices,tc_neuron),...
     'Color',spindling_green,'LineWidth',1.5*le_linewidth)
hold on
plot(highdose(1).time(timespan_indices),...
     highdose(1).TC_iTChing2010TC_hT(timespan_indices,tc_neuron),'r',...
     'LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
hold off
ylim([0 0.2])
% title('TC#1 T-current h_T')
% ylabel(sprintf('h_T (unitless)'))
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_h_T.fig'))
