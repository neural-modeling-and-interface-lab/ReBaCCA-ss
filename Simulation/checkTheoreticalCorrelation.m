clc,clear, close all

rng(0)

% Parameters
T = 10;          % Total time (seconds)
lambda1 = 1;     % Firing rate (spikes/second)
lambda2 = 1;
dt = 0.001;      % Time bin size (seconds)
sigma_vec = logspace(-3, -0.3, 20);  % Range of kernel bandwidths (0.01 to 1 s)

totalRepeat = 1000;
% Smooth and compute correlation for each bandwidth
for iRepeat = 1:totalRepeat
    % Generate two independent Poisson spike trains
    N1 = poissrnd(lambda1 * T);
    spike_times1 = sort(rand(N1, 1) * T);
    N2 = poissrnd(lambda2 * T);
    spike_times2 = sort(rand(N2, 1) * T);
    
    % Bin the spike times
    edges = 0:dt:T;
    n1 = histcounts(spike_times1, edges);
    n2 = histcounts(spike_times2, edges);
    
    % Compute rate vectors
    r1 = n1 / dt;
    r2 = n2 / dt;
    
    % Initialize correlation vector
    corr_vec = zeros(size(sigma_vec));
    for i = 1:length(sigma_vec)
        sigma = sigma_vec(i);
        % Define Gaussian kernel
        M = ceil(3 * sigma / dt);         % Kernel support: +/- 3*sigma
        t_kernel = (-M:M) * dt;
        kernel = (1 / (sigma * sqrt(2 * pi))) * exp(-t_kernel.^2 / (2 * sigma^2));
        kernel = kernel / sum(kernel);    % Normalize to sum to 1
        % Smooth the rate vectors
        smoothed_r1 = conv(r1, kernel, 'same');
        smoothed_r2 = conv(r2, kernel, 'same');
        % Compute correlation
        corr_vec(i) = abs(corr(smoothed_r1(:), smoothed_r2(:)));
    end
    total_corr(iRepeat, :)=corr_vec;
end

%% Plot the results
figure;

% varRho = sigma_vec * sqrt(2*pi) / T;
% y = sqrt(varRho * 2 / pi);

y = sqrt(sigma_vec / T) * (8 / pi)^0.25;
hold on
plot(sigma_vec * 1000, y, 'r', 'LineWidth', 1.5)

plot(sigma_vec * 1000, mean(total_corr), '-bo', 'LineWidth', 1.5);
xlabel('Kernel Bandwidth \sigma (ms)');
ylabel('Correlation');
title({['CCA of Two Independent Poisson Spike Trains'] ; [num2str(lambda1) 'Hz and ' num2str(lambda2) 'Hz']});
grid on
xlim([0 max(sigma_vec * 1000)])
ylim([0 0.3])

legend('Theoretical value', 'Simulation', 'Location', 'best')

export_pdf_figure(1, ['Fig/' num2str(lambda1) 'Hz and ' num2str(lambda2) 'Hz'], 0)