run_name = 'p1d1w2f2s2to4a7_20170126_simple_propo_ans';
%{
# Analyse Simulation file
- DEPENDENCIES:
    - This requires `p1d1w2f2s2to3r1_20170106_simple_propo.m` to have run, and
    have its data accessible
    - Run this TWICE in a live session to get consistent plotting:
        - Because MATLAB plots plots differently depending on if the plot is already
        in memory (aka visible in a figure window in the current session),
        versus if the plot has NOT already been generated/loaded into memory, in
        order for the ticks/etc. to be consistent upon saves to disk, you should
        run any plot-saving scripts TWICE, including once where the figures are
        already visible on screen, just to make sure it's printing consistently.
        Did I mention I hate matlab?
- Project 1:
    - Propofol PAC investigation
- Direction 1:
    - Thalamus-only propofol PAC modeling
- Writing item 2:
    - Journal article on results of p1d1
- Figure 2:
    - Simple baseline silence vs highdose propofol sustained alpha example
    simulation
- Subfigures 2-4:
    - Voltage traces of TC #1 under baseline and highdose, and zoom of high dose
- Analysis 4:
  - Illustrating baseline-silent depolarization, highdose-oscillating, using
  depolarized state 0.3, higher gh 0.0032, including rastergrams, over less
  time, and adding zoom of high dose traces. Now, removing ALL ticklabels so
  they'll be added at the right size in Inkscape later, since MATLAB is not to
  be trusted. Now converting everything to "fig" files for Michelle to work
  with, and restoring ticklabels.
- Date Created:
    - 20170126
- Inherits from:
    - 'p1d1w2f2s2to3a1_20170106_simple_propo_ans.m'
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Figure specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all figure backgrounds to white
set(gcf,'color','w');
% Erase current figures, since otherwise they may get written on
clf
close all

le_fontsize = 22;
le_linewidth = 1.5;
tc_neuron = 2;

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

% - Note that run-data directory name is different from above run-analysis name
data_name = 'p1d1w2f2s2to3r1_20170106_simple_propo';
data_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                  data_name);

if ~exist('baseline','var')
  baseline = load(strcat(data_dir, '/data/', 'study_sim1_data.mat'));
  highdose = load(strcat(data_dir, '/data/', 'study_sim2_data.mat'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Declare ranges to plot over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01;
downsample_factor = 6;
% Beginning of timespan we're interested in showing, in ms (after the initial
%   condition transient)
timebegin = 1;
% End of timespan we're interested in showing, in ms
timelimit = 500;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];

number_bins = 1000;

y_maximum = 60;
y_minimum = -90;

figno = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Plot baseline
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TC Trace
figno = figno + 1
h1_base = figure(figno)
plot(baseline(1).time(timespan_indices),...
     baseline(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% Not printing axes labels and titles in MATLAB because MATLAB is stupid and
%   resizes EVERYTHING differently even for almost IDENTICAL plots!!!
% ylabel('TC #1 Voltage [mV]')
% xlabel('Time [ms]')
% title('Baseline')
ylim([y_minimum y_maximum])
xlim([timebegin timelimit])
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_baseline_tc_vtrace.fig'))

fig1_base = get(handle(h1_base));
ax1_base = findobj(fig1_base.Children,'type','axes');
posn1 = get(ax1_base(1),'Position');

% Prepare rastergram plot
figno = figno + 1;
h1 = PlotData(baseline, 'variable','V',...
        'plot_type','rastergram');
fig1 = get(handle(h1));
ax1 = findobj(fig1.Children,'type','axes');

%
%% Grab TC rastergram Axes object to put it in its own figure
%
figno = figno + 1;
h10 = figure(figno);
% You cannot reset the (Axes) Child object of a Figure object to be an Axes
% object that is not already a Child. Instead, you must reset the Parent of the
% Axes object itself to be the Figure you want to own the Axes object. This is
% said in the MATLAB help Figure Properties entry for Children. UGH
% This is what moves the Axes object into a different figure hooray!
%
% Also, confusingly, `ax1(2)` is the TC rastergram, and `ax1(1)` is the RE
% rastergram.
set(ax1(2), 'Parent', h10)
set(ax1(2), 'Position', posn1)
set(gca,'FontSize',le_fontsize)

% set(ax1(2),'YTickLabel', '')
% ylabel(ax1(2),'')
% 
% set(ax1(2),'XTickLabel', '')
% set(ax1(2),'XTick', get(ax1_base(1),'XTick'))
% xlabel(ax1(2),'')
xlim([timebegin timelimit])

set(ax1(2),'TickDir','in')
set(ax1(2),'Box', 'on')
set(ax1(2),'BoxStyle', 'back')
set(ax1(2),'Layer', 'top')

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_highdose_tc_rastergram.fig'))

% RE Trace
figno = figno + 1;
figure(figno)
plot(baseline(1).time(timespan_indices),...
     baseline(1).RE_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% Not printing axes labels and titles in MATLAB because MATLAB is stupid and
%   resizes EVERYTHING differently even for almost IDENTICAL plots!!!
% ylabel('RE #1 Voltage [mV]')
% xlabel('Time [ms]')
% title('Baseline')
ylim([y_minimum y_maximum])
xlim([timebegin timelimit])
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_baseline_re_vtrace.fig'))

%
%% Grab RE rastergram Axes object to put it in its own figure
%
figno = figno + 1;
h11 = figure(figno);
set(ax1(1), 'Parent', h11)
set(ax1(1), 'Position', posn1)
set(gca,'FontSize',le_fontsize)

% set(ax1(1),'YTickLabel','')
% ylabel(ax1(1),'')
% 
% set(ax1(1),'XTickLabel', '');
% set(ax1(1),'XTick', get(ax1_base(1),'XTick'));
% set(ax1(1),'XColor', 'k')
% xlabel(ax1(1),'')
xlim([timebegin timelimit])

set(ax1(1),'TickDir','in')
set(ax1(1),'Box', 'on')
set(ax1(1),'BoxStyle', 'back')
set(ax1(1),'Layer', 'top')

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_highdose_re_rastergram.fig'))

% Close annoying big central rastergram from PlotData
close(h1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Plot highdose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TC Trace
figno = figno + 1;
h2_base = figure(figno)
plot(highdose(1).time(timespan_indices),...
     highdose(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% ylabel('TC #1 Voltage [mV]')
% xlabel('Time [ms]')
% title('High-dose')
ylim([y_minimum y_maximum])
xlim([timebegin timelimit])
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_highdose_tc_vtrace.fig'))

fig2_base = get(handle(h2_base));
ax2_base = findobj(fig2_base.Children,'type','axes');
posn2 = get(ax2_base(1),'Position');

% Create rastergram from PlotData
figno = figno + 1;
h2 = PlotData(highdose, 'variable','V',...
        'plot_type','rastergram');
fig2 = get(handle(h2));
ax2 = findobj(fig2.Children,'type','axes');

%
%% Grab TC rastergram Axes object to put it in its own figure
%
figno = figno + 1;
h20 = figure(figno);
% You cannot reset the (Axes) Child object of a Figure object to be an Axes
% object that is not already a Child. Instead, you must reset the Parent of the
% Axes object itself to be the Figure you want to own the Axes object. This is
% said in the MATLAB help Figure Properties entry for Children. UGH
% This is what moves the Axes object into a different figure hooray!
%
% Also, confusingly, `ax2(2)` is the TC rastergram, and `ax2(1)` is the RE
% rastergram.
set(ax2(2), 'Parent', h20)
set(ax2(2), 'Position', posn2)
set(gca,'FontSize',le_fontsize)

% set(ax2(2),'YTickLabel','')
% ylabel(ax2(2),'')
% 
% set(ax2(2),'XTickLabel', '')
% set(ax2(2),'XTick', get(ax2_base(1),'XTick'))
% xlabel(ax2(2),'')
xlim([timebegin timelimit])

set(ax2(2),'TickDir','in')
set(ax2(2),'Box', 'on')
set(ax2(2),'BoxStyle', 'back')
set(ax2(2),'Layer', 'top')

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_highdose_tc_rastergram.fig'))

% RE Trace
figno = figno + 1;
figure(figno)
plot(highdose(1).time(timespan_indices),...
     highdose(1).RE_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% ylabel('RE #1 Voltage [mV]')
% xlabel('Time [ms]')
% title('High-dose')
ylim([y_minimum y_maximum])
xlim([timebegin timelimit])
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_highdose_re_vtrace.fig'))

%
%% Grab RE rastergram Axes object to put it in its own figure
%
figno = figno + 1;
h21 = figure(figno);
set(ax2(1), 'Parent', h21)
set(ax2(1), 'Position', posn2)
set(gca,'FontSize',le_fontsize)

% set(ax2(1),'YTickLabel','')
% ylabel(ax2(1),'')
% 
% set(ax2(1),'XTickLabel', '')
% set(ax2(1),'XTick', get(ax2_base(1),'XTick'));
% set(ax2(1),'XColor', 'k')
% xlabel(ax2(1),'')
xlim([timebegin timelimit])

set(ax2(1),'TickDir','in')
set(ax2(1),'Box', 'on')
set(ax2(1),'BoxStyle', 'back')
set(ax2(1),'Layer', 'top')

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_highdose_re_rastergram.fig'))

% Close annoying big central rastergram from PlotData
close(h2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. Plot zoomed highdose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Beginning of timespan we're interested in showing, in ms (after the initial
%   condition transient)
timebegin = 440;
% End of timespan we're interested in showing, in ms
timelimit = 500;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];

% TC Trace
figno = figno + 1;
h3_base = figure(figno);
plot(highdose(1).time(timespan_indices),...
     highdose(1).TC_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% ylabel('TC #1 Voltage [mV]')
% xlabel('Time [ms]')
% title('High-dose')
ylim([y_minimum y_maximum])
xlim([timebegin timelimit])
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoomed_highdose_tc_vtrace.fig'))

fig3_base = get(handle(h3_base));
ax3_base = findobj(fig3_base.Children,'type','axes');
posn3 = get(ax3_base(1),'Position');

% Create rastergram from PlotData
figno = figno + 1;
h3 = PlotData(highdose, 'variable','V',...
        'plot_type','rastergram');
fig3 = get(handle(h3));
ax3 = findobj(fig3.Children,'type','axes');

%
%% Grab TC rastergram Axes object to put it in its own figure
%
figno = figno + 1;
h30 = figure(figno);
% You cannot reset the (Axes) Child object of a Figure object to be an Axes
% object that is not already a Child. Instead, you must reset the Parent of the
% Axes object itself to be the Figure you want to own the Axes object. This is
% said in the MATLAB help Figure Properties entry for Children. UGH
% This is what moves the Axes object into a different figure hooray!
%
% Also, confusingly, `ax3(2)` is the TC rastergram, and `ax3(1)` is the RE
% rastergram.
set(ax3(2), 'Parent', h30)
set(ax3(2), 'Position', posn3)
set(gca,'FontSize',le_fontsize)

% set(ax3(2),'YTickLabel','')
% ylabel(ax3(2),'')
% 
% set(ax3(2),'XTickLabel', '')
% set(ax3(2),'XTick', get(ax3_base(1),'XTick'))
% xlabel(ax3(2),'')
xlim([timebegin timelimit])

set(ax3(2),'TickDir','in')
set(ax3(2),'Box', 'on')
set(ax3(2),'BoxStyle', 'back')
set(ax3(2),'Layer', 'top')

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoomed_highdose_tc_rastergram.fig'))

% RE Trace
figno = figno + 1;
figure(figno)
plot(highdose(1).time(timespan_indices),...
     highdose(1).RE_V(timespan_indices,tc_neuron),...
     'k','LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
% ylabel('RE #1 Voltage [mV]')
% xlabel('Time [ms]')
% title('High-dose')
ylim([y_minimum y_maximum])
xlim([timebegin timelimit])
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoomed_highdose_re_vtrace.fig'))

%
%% Grab RE rastergram Axes object to put it in its own figure
%
figno = figno + 1;
h31 = figure(figno);
set(ax3(1), 'Parent', h31)
set(ax3(1), 'Position', posn3)
set(gca,'FontSize',le_fontsize)

% set(ax3(1),'YTickLabel','')
% ylabel(ax3(1),'')
% 
% set(ax3(1),'XTickLabel', '')
% set(ax3(1),'XTick', get(ax3_base(1),'XTick'));
% set(ax3(1),'XColor', 'k')
% xlabel(ax3(1),'')
xlim([timebegin timelimit])

set(ax3(1),'TickDir','in')
set(ax3(1),'Box', 'on')
set(ax3(1),'BoxStyle', 'back')
set(ax3(1),'Layer', 'top')

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_zoomed_highdose_re_rastergram.fig'))

% Close annoying big central rastergram from PlotData
close(h3)
