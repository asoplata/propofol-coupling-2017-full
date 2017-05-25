run_name = 'p1d1w2f7to8a5_20170126_triple_switch';
%{
# Analyse Simulation file
- DEPENDENCIES:
    - Requires 'p1d1w2f7to8r1_20170109_triple_switch.m' to have run
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
- Analysis 5:
    - Plotting single trough-max and peak-max voltage traces for fig7, and whole
    voltage trace and rastergram for fig8. Now also removing all text and
    ticklabels. Now restore ticklabels to make FIG files for michelle
- Date Created:
    - 20170109
- Inherits from:
    - 'p1d1w2f5s1to6a1_20170107_spindling_comparison_ans.m'
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Figure specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all figure backgrounds to white
set(gcf,'color','w');

le_fontsize = 22;
le_linewidth = 1.5;
cornflowerblue = [0.4 0.6 0.92]; % because fight club
tc_neuron = 4;

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

% data_name = 'p1d1w2f7to8r1_20170109_triple_switch';
data_name = 'p1d1w2f7to8r2_20170109_triple_switch';
% data_name = 'p1d1q36c5i3_20161106_triple_switch_sfn2016poster';
data_dir = strcat('/projectnb/crc-nak/asoplata/dynasim_data/',...
                  data_name);

if ~exist('result','var')
  load(strcat(data_dir, '/plots/', 'study_sim1_plot3_CalcPAC.mat'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Plot relevant/interesting raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dt = 0.01;
downsample_factor = 1;

number_bins = 1000;
figno = 0;

% From 'p1d1q42c1i1to4a1_20160919_cOFF_f3_repro_all_ans':
%   Judging from visual inspection, we can safely say the EFFECTIVE (as in, very
%     high probability of leading to burst) T-current window of
%     activation is limited to
tcurr_window.ceiling = -74;
tcurr_window.floor = -80;

% Beginning of timespan we're interested in showing, in ms (after the initial
%   condition transient)
timebegin = 1000;
% End of timespan we're interested in showing, in ms
timelimit = 6000;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];

figno = figno + 1;
figure(figno)
plot(result(1).time(timespan_indices),...
     result(1).TC_V(timespan_indices,tc_neuron),'k',...
     'LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
ylim([-90 60])
xlim([1000 6000])
% title('TC#1 Voltage traces')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_tmax_vtrace.fig'))

timebegin =  6000;
timelimit = 11000;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];
figno = figno + 1;
figure(figno)
plot(result(1).time(timespan_indices),...
     result(1).TC_V(timespan_indices,tc_neuron),'k',...
     'LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
ylim([-90 60])
xlim([6000 11000])
% title('TC#1 Voltage traces')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_pmax_vtrace.fig'))

timebegin =  1000;
timelimit = 11000;
timespan_indices = [round(timebegin/(dt*downsample_factor)):...
                    round(timelimit/(dt*downsample_factor))];
figno = figno + 1;
figure(figno)
plot(result(1).time(timespan_indices),...
     result(1).TC_V(timespan_indices,tc_neuron),'k',...
     'LineWidth',le_linewidth)
set(gca,'FontSize',le_fontsize)
% set(gca,'XTickLabel','')
% set(gca,'YTickLabel','')
ylim([-90 60])
xlim([1000 11000])
% title('TC#1 Voltage traces')
% ylabel('Voltage in mV')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_whole_vtrace.fig'))

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% 2. Plot coupling-gram
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colormap jet
% 
% clims = [0.02 0.20];
% 
% %% Rearrange modulogram matrix so that it looks cleaner (since SWO phase is
% %  relative/where the modulo starts doesn't matter).
% %
% % Using MUA Modulograms
% 
% n_bins = 10;
% offset = 1;
% 
% % for TC_V_Modulograms_MUA, "y"/rows are time, "x"/columns are phases
% rearr_matrix = zeros(size(result(1).TC_V_Modulograms_MUA));
% rearr_matrix(:,1:offset) = result(1).TC_V_Modulograms_MUA(:,((n_bins-offset)+1):n_bins);
% rearr_matrix(:,(offset+1):n_bins) = result(1).TC_V_Modulograms_MUA(:,1:(n_bins-offset));
% 
% figno = figno + 1;
% h1 = figure(figno)
% % subplot 211
% % imagesc(result(1).TC_V_Modulation_Time_Points',...
% %         result(1).TC_V_Modulogram_Bin_Locations,...
% %         result(1).TC_V_Modulograms_MUA',...
% %         clims)
% imagesc(linspace(0,11,19),...
%         result(1).TC_V_Modulogram_Bin_Locations,...
%         rearr_matrix',...
%         clims)
%         % result(1).TC_V_Modulograms_MUA',...
% cb1 = colorbar
% set(gca,'FontSize',le_fontsize)
% % cb1.Label.String = 'Alpha Power Modulation Index'; % seriously matlab?
% % xlabel('Time in sec')
% % ylabel('SWO phase in radians')
% % subplot 212
% % plot(result(1).time, result(1).TC_V)
% % print(h1,'TC_V_Modulograms_MUA','-dpng')
% saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
%        'TC_V_Modulograms_MUA.fig'))
% 
% % figno = figno + 1;
% % h2 = figure(figno)
% % colormap jet
% % % subplot 211
% % % imagesc(result(1).TC_V_Modulation_Time_Points,...
% % %         result(1).TC_V_Modulogram_Bin_Locations,...
% % %         result(1).TC_V_Modulograms_SUA.modulograms_single(1).matrix',...
% % %         clims)
% % % Hard-coding the time axes since my CalcPAC is bugged and producing INCORRECT
% % %   `TC_V_Modulation_Time_Points`
% % imagesc(linspace(0,11,19),...
% %         result(1).TC_V_Modulogram_Bin_Locations,...
% %         result(1).TC_V_Modulograms_SUA.modulograms_single(1).matrix',...
% %         clims)
% % cb2 = colorbar
% % cb2.Label.String = 'Modulation Index'; % seriously matlab?
% % xlabel('Time in sec')
% % ylabel('SWO phase in radians')
% % % subplot 212
% % % plot(result(1).time, result(1).TC_V)
% % % print(h2,'TC_V_Modulograms_SUA','-dpng')
% % saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
% %        'TC_V_Modulograms_SUA.fig'))

%% Rastergram
figno = figno + 1;
h1 = PlotData(result(1), 'variable','V',...
        'plot_type','rastergram',...
        'format','png');
fig1 = get(handle(h1));
ax1 = findobj(fig1.Children,'type','axes');
% Can NOT see in Workspace the property 'YLabel' in e.g. 'ax1(1)', even though
% it's there

% Finally!
% Removed `text_string` printed onto rastergram by PlotData via commenting out
% uses of that actual variable in PlotData code
% ylabel(ax1(1),'')
% ylabel(ax1(2),'')
% xlabel(ax1(1),'')
xlim(ax1(1),[timebegin timelimit])
xlim(ax1(2),[timebegin timelimit])

% set(ax1(1),'YTickLabel','')
% set(ax1(1),'XTickLabel','')
% set(ax1(1),'TickDir','both')
% 
% set(ax1(2),'YTickLabel','')
% set(ax1(2),'TickDir','both')

set(gcf,'color','w');

saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_rastergram.fig'))
