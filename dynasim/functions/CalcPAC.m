function data = CalcPAC(data,varargin)
%% data = CalcPAC(data,'option',value)
% Inputs:
%   data - DynaSim data structure (see CheckData)
%   required key-value pair options:
%     'variable' - name of field containing data you want to calculate phase
%                  amplitude coupling for, both over individual cells and
%                  averaged across all cells. Note that this can be either a
%                  simple, common variable like "V", for which the variables in
%                  the data will be searched using this as a regular expression,
%                  or it can be a properly-namespaced variable like "TC_V".
%   other key-value pair options:
%     'peak_area_width' - Hz, size of frequency bin (centered on peak) over which to calculate area under spectrum (default: 5)
%     'exclude_data_flag' - whether to remove simulated data from result structure (default: 0)
%
% Outputs:
%   data: data structure, same as the input structure, but with added fields:
%     data(which_struct).VARIABLE_Modulation_Indices_SUA.modulation_indices_single:
%       matrix of (Modulation Index across each window) x (Number of cells)
%     data(which_struct).VARIABLE_Modulograms_SUA.modulograms_single(which_cell).matrix
%     data(which_struct).VARIABLE_Modulation_Indices_MUA.modulation_indices_mean
%     data(which_struct).VARIABLE_Modulograms_MUA.modulograms_mean

%         data.VARIABLE_Power_SUA.PeakFreq: frequency of spectral power (one value per cell)
%         data.VARIABLE_Power_SUA.PeakArea: area under spectrum around peak (one value per cell)
%   NOTE: for populations: spectrum of the mean waveform is stored in 
%         data.VARIABLE_Power_MUA.Pxx. population mean spectrum of the individual
%         waveforms can be calculated as mean(data.VARIABLE_Power_MUA.Pxx,2).
%
% organization scheme for spectral results:
%
% note:
% "variable" can be specified as the name of a variable listed in
% data.labels, a cell array of string listing variable names, or as a
% regular expression pattern for identifying variables to process.
% See SelectVariables for more info on supported specifications.
%
% Examples:
% s=[];
% s.populations(1).name='E';
% s.populations(1).equations='dv[2]/dt=@current+10; {iNa,iK}; v(0)=-65';
% s.populations(2).name='I';
% s.populations(2).equations='dv/dt=@current+10; {iNa,iK}; v(0)=-65';
% data=SimulateModel(s,'tspan',[0 1000]);
% data=CalcPAC(data,'variable','v');
% % Plot the spectrum of the E-cell average population voltage
% figure; plot(data.E_v_Power_MUA.frequency,data.E_v_Power_MUA.Pxx); 
% xlabel('frequency (Hz)'); ylabel('power'); xlim([0 200]);
%
% See also: PlotPower, AnalyzeStudy, SimulateModel, CheckData, SelectVariables

%% 1.0: Check inputs
options=CheckOptions(varargin,{...
  'variable',[],[],...
  'exclude_data_flag',0,{0,1},...
  'slow_freq_range',[0.5 1.5],[],...
  'fast_half_power_freqs',[5 16],[],...
  'window_time_length', 2.0,[],...
  'window_time_overlap',1.5,[],...
  'n_bins',18,[],...
  },false);

% Note: calling CheckData() at beginning enables analysis function to
%   accept data matrix [time x cells] in addition to DynaSim data structure.
data = CheckData(data);

% use AnalyzeStudy to recursively call CalcPAC on each data set
if numel(data)>1
  data = AnalyzeStudy(data,@CalcPAC,varargin{:});
  return;
end

if ~isfield(data,'results')
  data.results={};
end

%% 2.0: Initialize values
% time parameters
dt = data.time(2)-data.time(1); % time step

% frequency parameters
Fs = fix(1/(dt/1000)); % effective sampling frequency

% Compute window properties
window_properties.window_space_length = ceil(options.window_time_length*Fs);
window_properties.window_space_overlap = ceil(options.window_time_overlap*Fs);
window_properties.n_windows = fix((size(data.time,1) - window_properties.window_space_overlap) / ...
                                  (window_properties.window_space_length - window_properties.window_space_overlap));

% 
modulogram_bin_locations = zeros(options.n_bins,1);
bin_size = 2*pi/options.n_bins;
for kk = 1:options.n_bins
    modulogram_bin_locations(kk) = -pi+(kk-1)*bin_size;
end

modulation_time_points = (options.window_time_length / 2) + ...
                         [0:(fix(length(data.time)/Fs)-1)]*...
                           (options.window_time_length - options.window_time_overlap);
% From given simple variable names (e.g. "V"), find properly-namespaced "true
%   names" of the variables in the data that contain those variables (e.g. "TC_V")
options.variable=SelectVariables(data(1).labels,options.variable);

%% 3.0: Create faster frequency filter
% Only the filter for the faster frequency range (default alpha 8~13 Hz) is
%   created here. This filter has been tested with both realistic- and ideal-PAC
%   data, and seems to do the job well. Using `designfilt` like this appears to
%   be something close to the "recommended" way to construct filters in MATLAB.
%
% Note that you must supply the Half Power Frequencies, NOT the cutoff
%   frequencies of whichever band you're interested in. Using cutoff frequencies
%   so low on the frequency domain, combined with such a high sampling rate,
%   seemed to corrupt the filter, and so setting the Half Power Frequencies
%   instead seemed to work much better. Practically speaking, this just means
%   that you can set the Half Power Frequencies to be a few Hz broader than the
%   band you're primarily interested in, and you'll probably have a good enough
%   filter.
%
% Note that you can visually inspect the filter if you want by calling
%   `fvtool(fast_filter)`.
%
% Unfortunately, the same cannot be said for my slower frequency range
%   (default Slow Wave, 0.1-1.5 Hz) filter far below. Surprisingly, I've never
%   been able to construct an F/IIR filter using this method, bandpass or
%   lowpass, for that frequency range at the high sampling rate we use that
%   didn't fail catastrophically (this includes trying to build a Slow Wave
%   filter in Python!). However, the use of `idealfilter` in the child function
%   `compute_single_pac()` below is the ONLY thing I could get to work for the
%   slower frequency range -- that said, it works GREAT. It may require
%   additional MATLAB toolboxes, however. The slower frequency `idealfilter`
%   also does not need to be designed before it is applied, and so it is only
%   used in the single call inside `compute_single_pac`, rather than being
%   created here prior to use.
%
% NOTE: Frequencies used are in Hertz, even though MATLAB does not explicitly
%   say that!
fast_filter = designfilt('bandpassiir','FilterOrder', 10, ...
                    'HalfPowerFrequency1',options.fast_half_power_freqs(1),...
                    'HalfPowerFrequency2',options.fast_half_power_freqs(2),...
                    'SampleRate', Fs);

%% 4.0: Main loop!
for vv=1:length(options.variable)
  %% 4.1: Isolate the variable-specific data set
  var = options.variable{vv};
  dat = data.(var);
  n_cells = size(dat,2);

  %% 4.2: Pre-allocate result data structures
  %   Supposedly this http://stackoverflow.com/a/13664481 is fast, easy, and good
  modulograms_single = repmat(struct('matrix', zeros(options.n_bins, window_properties.n_windows)),...
                              n_cells, 1);
  modulation_indices_single = repmat(struct('time_series', zeros(window_properties.n_windows, 1)),...
                              n_cells, 1);

  % -----------------------------------------------------
  %% 4.3 Calculate Modulation Index and Modulograms for Single Unit Activity, looping over cells:
  for jj = 1:n_cells
    % Note: pass in time in SECONDS, not ms
    [modulograms_single(jj).matrix, ...
     modulation_indices_single(jj).time_series] = ...
         compute_single_pac(dat(:,jj), (data.time)./1000, options,...
         window_properties, fast_filter);
  end
  % -----------------------------------------------------
  %% 4.4 Calculate Modulation Index and Modulograms for Multi-Unit Activity, looping over cells:
  if n_cells == 1
    % same as SUA
    modulograms_mean = modulograms_single(1).matrix;
    modulation_indices_mean = modulation_indices_single(1).time_series;
  else
    % calculate MUA
    [modulograms_mean, ...
     modulation_indices_mean] = ...
         compute_single_pac(mean(dat,2), (data.time)./1000, options,...
         window_properties, fast_filter);
  end

  %% 4.5 Add results to data structure
  data.([var '_Modulation_Indices_SUA']).modulation_indices_single = modulation_indices_single;
  data.([var '_Modulograms_SUA']).modulograms_single = modulograms_single;
  data.([var '_Modulation_Indices_MUA']) = modulation_indices_mean;
  data.([var '_Modulograms_MUA']) = modulograms_mean;
  data.([var '_Modulation_Time_Points']) = modulation_time_points;
  data.([var '_Modulogram_Bin_Locations']) = modulogram_bin_locations;

  if ~ismember([var '_Modulation_Indices_SUA'],data.results)
    data.results{end+1}=[var '_Modulation_Indices_SUA'];
  end
  if ~ismember([var '_Modulograms_SUA'],data.results)
    data.results{end+1}=[var '_Modulograms_SUA'];
  end
  if ~ismember([var '_Modulation_Indices_MUA'],data.results)
    data.results{end+1}=[var '_Modulation_Indices_MUA'];
  end
  if ~ismember([var '_Modulograms_MUA'],data.results)
    data.results{end+1}=[var '_Modulograms_MUA'];
  end
  if ~ismember([var '_Modulation_Time_Points'],data.results)
    data.results{end+1}=[var '_Modulation_Time_Points'];
  end
  if ~ismember([var '_Modulogram_Bin_Locations'],data.results)
    data.results{end+1}=[var '_Modulogram_Bin_Locations'];
  end

  % Delete data if desired
  if options.exclude_data_flag
    for oo = 1:length(data.labels)
      data=rmfield(data,data.labels{oo});
    end
  end
end
end % function

function [modulogram_matrix, modulation_index_timeseries] = ...
    compute_single_pac(data_vector, time_vector, opts, ...
    window, faster_filter)

  % Pre-allocate result data structures
  modulogram_matrix = zeros(opts.n_bins, window.n_windows);
  modulation_index_timeseries = zeros(1, window.n_windows);

  try
    %% 1.0: Filter the data
    %% 1.1: Create and apply slower frequency filter simultaneously
    % NOTE: `slow_data` is a TIMESERIES objects, not a regular array.
    % NOTE: Need to convert time series into Seconds, since that's what it wants
    slow_data = idealfilter(timeseries(data_vector,time_vector),...
                            [opts.slow_freq_range(1),...
                             opts.slow_freq_range(2)],...
                            'pass');
    post_filter_slow = slow_data.Data;

    %% 1.2: Apply faster frequency filter
    % `filtfilt` NEEDS the data to be of the "double" type instead of
    % "single", which you have to do manually, and even though `filter`
    % doesn't care? Am I taking crazy pills???
    post_filter_fast = filtfilt(faster_filter, double(data_vector));

    % todo: HERE is probably where you'd want to use any of the available PAC
    % libraries, like the Eden/Kramer GLMCFC, since the signal has been
    % filtered, but not transformed into the analytic signal, and not windowed
    % yet (to get a phase-amplitude coupling "time-series" so to speak).

    %% 2.0 Get slower frequency phase angle and faster frequency amplitude from the
    %  respective signals via Hilbert Transform
    phi = angle(hilbert(post_filter_slow));
    amp = abs(hilbert(post_filter_fast));

    %% 3.0 Construct the windows across the timeseries, and find their
    %  corresponding data
    %    Adapted & taken from Angela Onslow's 'pac_code_best/window_data.m' of
    %    her MATLAB Toolbox for Estimating Phase-Amplitude Coupling from
    %
    %    http://www.cs.bris.ac.uk/Research/MachineLearning/pac/
    %
    % `idx` are each set of the INDICES we will use to "grab" the values
    %   in each of the windows
    idx = bsxfun(@plus,...
                (1:window.window_space_length)',...
                1+(0:(window.n_windows - 1))*(window.window_space_length - window.window_space_overlap)) - 1;

    % Note: If you're running out of memory, look here -- you can use 'idx' as
    %   an indexing matrix in the major iteration below, instead of lazily
    %   re-copying all the data several times to their individual windows here as
    %   I've done. This is lazy, but also makes the code further on much easier
    %   to understand.
    amp_window = amp(horzcat(idx(:, 1:size(idx,2))));
    phi_window = phi(horzcat(idx(:, 1:size(idx,2))));

    %% 4.0 Bin the faster frequency's amplitude in the slower's phase bins (this
    %  is the important part)
    %    Adapted & taken from Adriano Tort's
    %    'Neurodynamics-master/16ch/Comodulation/ModIndex_v1.m' of the
    %    'Neurodynamics-Toolbox' repo on Github, at
    %
    %    https://github.com/cineguerrilha/Neurodynamics
    %
    % Define the beginning (not the center) of each bin (in rads)
    phi_bin_beginnings = zeros(opts.n_bins,1);
    bin_size = 2*pi/opts.n_bins;
    for kk = 1:opts.n_bins
        phi_bin_beginnings(kk) = -pi+(kk-1)*bin_size;
    end

    % Iterate over each individual window
    for mm = 1:size(idx,2)
        % Now we compute the mean amplitude in each phase bin:
        amp_means = zeros(1,opts.n_bins);
        for nn = 1:opts.n_bins
            amp_means(nn) = ...
                mean(amp_window((phi_window(:,mm) >= phi_bin_beginnings(nn)) & ...
                                (phi_window(:,mm) <  phi_bin_beginnings(nn)+bin_size),mm));
        end
        modulogram_matrix(:,mm) = (amp_means/sum(amp_means))';

        % Quantify the amount of amp modulation by means of a normalized
        %   entropy index (Tort et al PNAS 2008):
        modulation_index_timeseries(:,mm) = ...
            (log(opts.n_bins)-(-sum((amp_means/sum(amp_means)).*...
             log((amp_means/sum(amp_means))))))/log(opts.n_bins)';
    end

    % So each row corresponds to a time point
    modulogram_matrix = modulogram_matrix';

    % % Uncomment this if you want some helpful debugging plots. These show the
    % % most important information for a single call of this function.
    % figure(10)
    % subplot 511
    % plot(time_vector, data_vector)
    % title('pre-filtered data')
    % subplot 512
    % plot(time_vector, post_filter_fast)
    % title('POST-filtered alpha data')
    % subplot 513
    % plot(time_vector, post_filter_slow)
    % title('POST-filtered SWO data')
    % subplot 514
    % plot(time_vector, amp)
    % title('alpha amplitude from hilbert xform')
    % subplot 515
    % plot(time_vector, phi)
    % title('SWO phase from hilbert xform')
    % xlabel('time')

    % figure(11)
    % subplot 211
    % imagesc(modulogram_matrix)
    % title('modulogram matrix')
    % colorbar
    % subplot 212
    % plot(modulation_index_timeseries)
    % title('modulation index MI timeseries')
    % xlabel('time')

  catch err
    fprintf('\n!\n!\n!\n whoa somethings wrong, DEBUG\n!\n!\n!\n')
    DisplayError(err)
  end
end
