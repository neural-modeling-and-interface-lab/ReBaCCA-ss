clc, clear, close all

rng(0)  % Set random seed for reproducibility

% Parameters
totalNeuron1 = 1;
totalNeuron2 = 1;
T = 10;          % Total time (seconds)
lambda1 = 4 * ones(1, totalNeuron1);     % Firing rate (spikes/second) for first spike train
lambda2 = 4 * ones(1, totalNeuron2);     % Firing rate (spikes/second) for second spike train
dt = 0.001;      % Time bin size (seconds)
sigma_vec = logspace(-3, log10(0.15), 10);  % Range of kernel bandwidths
totalRepeat = 500;  % Number of repeats
edges = 0:dt:T;  % Time bin edges
num_bins = length(edges) - 1;  % Number of time bins (10000)

% Step 2: Smooth and compute correlation for each bandwidth
total_vaccc = zeros(totalRepeat, length(sigma_vec));
total_corr1 = zeros(totalRepeat, length(sigma_vec));  % Store correlations
total_corr2 = zeros(totalRepeat, length(sigma_vec));  % Store correlations
total_corr = zeros(totalRepeat, length(sigma_vec));
for i = 1:length(sigma_vec)

    sigma = sigma_vec(i);
    % Define Gaussian kernel
    M = ceil(3 * sigma / dt);         % Kernel support: +/- 3*sigma
    t_kernel = (-M:M) * dt;
    kernel = (1 / (sigma * sqrt(2 * pi))) * exp(-t_kernel.^2 / (2 * sigma^2));
    kernel = kernel / sum(kernel);    % Normalize kernel to sum to 1
    
    for repeat = 1:totalRepeat
        % Step 1: Generate two spike matrices
        spike_matrix1 = zeros(totalNeuron1, num_bins);  % Matrix for first spike train
        spike_matrix2 = zeros(totalNeuron2, num_bins);  % Matrix for second spike train
        
        for iNeuron = 1:totalNeuron1
            % Generate spike times for first spike train
            N1 = poissrnd(lambda1(iNeuron) * T);
            spike_times1 = sort(rand(N1, 1) * T);
            n1 = histcounts(spike_times1, edges);
            spike_matrix1(iNeuron, :) = n1;  % Store binned counts in matrix
        end
        for iNeuron = 1:totalNeuron2
            % Generate spike times for second spike train
            N2 = poissrnd(lambda2(iNeuron) * T);
            spike_times2 = sort(rand(N2, 1) * T);
            n2 = histcounts(spike_times2, edges);
            spike_matrix2(iNeuron, :) = n2;  % Store binned counts in matrix
        end

        smooth_matrix1 = zeros(totalNeuron1, num_bins); 
        for iNeuron = 1:totalNeuron1
            smooth_matrix1(iNeuron, :) = conv(spike_matrix1(iNeuron, :), kernel, 'same');
        end
        smooth_matrix2 = zeros(totalNeuron2, num_bins); 
        for iNeuron = 1:totalNeuron2
            smooth_matrix2(iNeuron, :) = conv(spike_matrix2(iNeuron, :), kernel, 'same');
        end
        tol = 1e-3;
        maxIter = 100;
        percentVar = 1 - 1e-3;
        nComponents = min([totalNeuron1 totalNeuron2]);
            
        [A, B, U, V, P, Q, ~, objComp, objCompTrue, varExplained] = ...
            continuumRegressionMulti(smooth_matrix1', smooth_matrix2', 0.5, nComponents, tol, maxIter, percentVar);
        total_vaccc(repeat, i) =  varExplained*diag(corr(U,V));
        
        [Acca, Bcca, r, Ucca, Vcca] = canoncorr(smooth_matrix1', smooth_matrix2');

        % Compute absolute correlation
        % total_corr1(repeat, i) = r(1);
        % total_corr2(repeat, i) = r(2);
        total_corr(repeat, i) = mean(r);
    end
end

%% Plot the spike rasters
figure
plotSpikeRaster(spike_matrix1, size(spike_matrix1,1),T, '', size(spike_matrix1,1));
xticks([]);yticks([])
box on
xlabel('');ylabel('')
export_pdf_figure(gcf, ['Temp Fig/' num2str(totalNeuron1) ' neurons vs ' num2str(totalNeuron2) ' neurons spike matrix 1'], 0)

plotSpikeRaster(spike_matrix2, size(spike_matrix2,1),T, '', size(spike_matrix2,1));
xticks([]);yticks([])
box on
xlabel('');ylabel('')
export_pdf_figure(gcf, ['Temp Fig/' num2str(totalNeuron1) ' neurons vs ' num2str(totalNeuron2) ' neurons spike matrix 2'], 0)

%% Plot the results
figure;

% varRho = mean([totalNeuron1 totalNeuron2]) * sigma_vec * sqrt(2*pi) / T ;
% y = sqrt(varRho * 2 / pi);

y = (8 / pi) ^ 0.25 * sqrt(sigma_vec / T * mean([totalNeuron1 totalNeuron2]));

hold on
plot(sigma_vec * 1000, y, 'r', 'LineWidth', 1.5)

% plot(sigma_vec * 1000, mean(total_corr1), '-bo', 'LineWidth', 1.5);
% plot(sigma_vec * 1000, mean(total_corr2), '-bo', 'LineWidth', 1.5);
% plot(sigma_vec * 1000, mean([total_corr1;total_corr2]), '-ob', 'LineWidth', 1.5);
plot(sigma_vec * 1000, mean(total_corr), '-ob', 'LineWidth', 1.5);
% plot(sigma_vec * 1000, mean(total_vaccc), '-go', 'LineWidth', 1.5);

xlabel('Kernel Bandwidth \sigma (ms)');
ylabel('Correlation');
title({'CCA of Two Independent Poisson Spike Matrices' ; [num2str(totalNeuron1) ' neurons vs ' num2str(totalNeuron2) ' neurons']});
grid on
xlim([0 max(sigma_vec * 1000)])
ylim([0 1])

legend('Theoretical value', 'Simulated value', 'Location', 'best')

export_pdf_figure(gcf, ['Temp Fig/' num2str(totalNeuron1) ' neurons vs ' num2str(totalNeuron2) ' neurons'], 0)