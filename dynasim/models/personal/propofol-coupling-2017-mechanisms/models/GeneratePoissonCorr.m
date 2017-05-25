function psps = GeneratePoissonCorr(no_cells, inputs_per_cell, rate, tau_i, tau_1, tau_d, tau_r, T, dt, jitter_sigma)
%{
# DESCRIPTION
This code, similar to `GeneratePoissonUncorr.m`, produces a set of spiketrains
that are convolved with some EPSP waveform, the only difference being that
you can control how similar/correlated the spiketrains are.

It does this by producing one random spiketrain, copying it for each cell,
and then applying some amount of 'jitter', drawn from a Normal/Gaussian
distribution with mean 0 and standard deviation 'jitter_sigma', to every
spike. Increasing the jitter will decrease the correlation between the
spiketrains.

# INPUTS
jitter_sigma : the standard deviation of the Gaussian distribution you want
applied to the jitter of individual spikes, in units of whatever your 'dt'
is. So if your 'dt' is 0.01 ms, setting this to 500 means the spike jitter
will have a standard deviation of 5 ms.
%}

t = 0:dt:T;

%% Generate the main spiketrain that others will be based off of
spikes = rand(inputs_per_cell, size(t,2));
spikes = spikes < rate*dt/1000;

if inputs_per_cell > 1
    combined_inputs = sum(spikes);
else
    combined_inputs = spikes;
end
combined_times  = find(combined_inputs);

%% Jitter
% Jitter the spike times a bit
identical_times = repmat(combined_times, no_cells, 1);
jitter_matrix = normrnd(0, jitter_sigma, size(identical_times));
jittered_times = identical_times + jitter_matrix;

% Correct any out-of-time-bounds jittered times, and integer-ize them,
% sinced they will be used as indices
jittered_times(jittered_times < 1) = 1;
jittered_times(jittered_times > (T/dt)) = T/dt;
jittered_times = round(jittered_times);

%% Transform the jittered spike times back into spike trains
% Because MATLAB apparently doesn't like slicing...
jittered_neuron_indices = zeros(size(jittered_times));
for ii=1:no_cells
    jittered_neuron_indices(ii,:) = ii;
end

% Write the spike times into spike trains
spike_arrivals= zeros(no_cells, length(combined_inputs));
spike_arrivals(sub2ind(size(spike_arrivals), ...
               reshape(jittered_neuron_indices',[],1), ...
               reshape(jittered_times',[],1))) = 1;

%% Shape the PSPs
psps = nan(size(spike_arrivals)); % Initialize the EPSPs experienced by each cell.

% Shape the actual EPSP for spikes at time t = 0.
psp = tau_i*(exp(-max(t - tau_1,0)/tau_d) - exp(-max(t - tau_1,0)/tau_r))/(tau_d - tau_r);
psp = psp(psp > eps);    %?
psp = [zeros(1,length(psp)) psp]; %?

%% Convolve everything!
for c = 1:no_cells
    psps(c,:) = conv(spike_arrivals(c,:),psp,'same');
end
