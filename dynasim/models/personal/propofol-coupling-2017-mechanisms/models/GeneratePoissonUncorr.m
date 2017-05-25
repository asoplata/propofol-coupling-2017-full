function psps = GeneratePoissonUncorr(no_cells, inputs_per_cell, rate, tau_i, tau_1, tau_d, tau_r, T, dt)
% Derived from Ben Poletta's "multi_Poisson.m"

t = 0:dt:T;

% EPSP for spikes at time t = 0.
psp = tau_i*(exp(-max(t - tau_1,0)/tau_d) - exp(-max(t - tau_1,0)/tau_r))/(tau_d - tau_r);
psp = psp(psp > eps);    %?
psp = [zeros(1,length(psp)) psp]; %?

no_inputs = inputs_per_cell*no_cells;

C = repmat(eye(no_cells), 1, no_inputs/no_cells);

spikes = rand(no_inputs, size(t,2));
spikes = spikes < rate*dt/1000;

spike_arrivals = C*spikes; % Calculating presynaptic spikes for each cell.

psps = nan(size(spike_arrivals)); % Calculating EPSP experienced by each cell.
for c = 1:no_cells
    psps(c,:) = conv(spike_arrivals(c,:),psp,'same');
end
