run_name = 'p1d1w2f6s1a1_20170108_tau_t_ans';
%{
# Analyse Simulation file
- DEPENDENCIES:
    - Requires 'p1d1w2f5s14r1_20170108_tau_t.m' to have run, AND you to have
    visually inspected firing frequencies from those simulations, as below
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
- Analysis 1:
    - Just plotting PRE-analyzed data
- Date Created:
    - 20170108
- Inherits from:
    -'p1d1q43c1i1a1_20160919_cOFF_f3_TcurrTau_ans.m'
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 0. Figure specifications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set all figure backgrounds to white
set(gcf,'color','w');

le_fontsize = 22;
le_linewidth = 2;
cornflowerblue = [0.4 0.6 0.92]; % because fight club
tc_neuron = 1;

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
%% 1. Plot TC frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From visual inspection of related sims
%   'p1d1w2f6s1r1_20170108_tau_t'

stats.tauHproportion = [0.25, 0.5, 0.75,  1, 1.25, 1.5, 1.75, 2, 3, 5, 10];
stats.frequency =      [   0,   0,   12, 13,   12,  11,   10, 0, 8, 8,  5];

figno = 1;
figure(figno)
semilogx(stats.tauHproportion, stats.frequency,'*-k','LineWidth',le_linewidth)

hold on
p1 = semilogx(  ones(20,1), linspace(0,14,20), 'r','LineWidth',le_linewidth)
hold off

legend([p1],{sprintf('Default')},'Location','northwest')

ylim([0 14])
set(gca,'FontSize',le_fontsize)
% title(sprintf('Network frequency across Tau_h_T proportion'))
% xlabel('Tau_h_T proportion multiplier')
% ylabel('Hz')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_tauH_freqs.fig'))
