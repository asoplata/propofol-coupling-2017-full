run_name = 'p1d1w2f6s2a1_20170108_spm_freqs_ans';
%{
# Analyse Simulation file
- DEPENDENCIES:
    - Requires 'p1d1w2f5s15r1_20170108_spm_freqs.m' to have run, AND you to have
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
- Subfigure 2:
    - Plot of network frequency across spm changes
- Date Created:
    - 20170108
- Inherits from:
    - 'p1d1q43c1i2a1_20160919_cOFF_f3_spm_ans.m'
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
%   ('p1d1w2f6s2r1_20170108_spm_freqs')

stats.spm =       [ 1,  2,  3,  4,  5,  6,  7,  8];
stats.frequency = [ 0,  11, 12, 11, 9,  8,  7.5,  7];

figno = 1;
figure(figno)
p0 = plot(stats.spm, stats.frequency,'*-k','LineWidth',le_linewidth)

hold on
p1 = plot(  ones(20,1), linspace(0,12,20), 'm','LineWidth',le_linewidth)
p2 = plot(2*ones(20,1), linspace(0,12,20), 'b','LineWidth',le_linewidth)
p3 = plot(3*ones(20,1), linspace(0,12,20), 'r','LineWidth',le_linewidth)
hold off

legend([p1 p2 p3],{'baseline', 'low-dose','high-dose'},...
    'Location','southeast')
ylim([0 12])
set(gca,'FontSize',le_fontsize)
% title(sprintf('Network frequency across extreme\npropofol GABA_A potentiation'))
% xlabel('Propofol GABA_A multiplier')
% ylabel('Hz')
saveas(gcf, strcat(outgoing_plots_dir, run_name, '_subfig', int2str(figno),...
       '_spm_freqs.fig'))
