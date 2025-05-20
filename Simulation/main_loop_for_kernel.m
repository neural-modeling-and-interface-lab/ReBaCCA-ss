clc, clear, close all
tic;

totalRepeat = 8;  % 2025-03-16 Original 10
type = 'multiply';

kernel_width_pool = logspace(log10(1), log10(100), 24);
alphaPool = 0.5;
tol = 1e-3;
maxIter = 100;
percentVar = 1 - 1e-3;

% **Preallocate result arrays**
nWidths = length(kernel_width_pool);
totalIndex1 = zeros(nWidths, totalRepeat);
totalIndex2 = zeros(nWidths, totalRepeat);
totalIndex_perm1 = zeros(nWidths, totalRepeat);
totalIndex_perm2 = zeros(nWidths, totalRepeat);
totalError1_repeat = zeros(totalRepeat, nWidths, 2);
totalError2_repeat = zeros(totalRepeat, nWidths, 2);
totalExplained1_repeat = zeros(totalRepeat, nWidths, 2);
totalExplained2_repeat = zeros(totalRepeat, nWidths, 2);
totalCorr1_repeat = zeros(totalRepeat, nWidths, 2);
totalCorr2_repeat = zeros(totalRepeat, nWidths, 2);
totalCorrCCA_repeat = zeros(totalRepeat, nWidths);
totalCorrCCA_perm_repeat = zeros(totalRepeat, nWidths);
totalRA_repeat = zeros(totalRepeat, nWidths);
totalRA_perm_repeat = zeros(totalRepeat, nWidths);
objComponent_repeat = cell(totalRepeat, nWidths);
objComponent_repeat_perm = cell(totalRepeat, nWidths);
A_total_repeat = cell(totalRepeat, nWidths);
B_total_repeat = cell(totalRepeat, nWidths);
U_total_repeat = cell(totalRepeat, nWidths);
V_total_repeat = cell(totalRepeat, nWidths);

for alpha = alphaPool
    for iRepeat = 1:totalRepeat
        % Load data once per repeat
        load(['Data/Simulated data_' num2str(iRepeat)])
        
        nComponents = min(size(concatenated_spike_trains1,1), size(concatenated_spike_trains2,1));
        extraction_start = 1500;
        extraction_end = 3500;
        extraction_trial_length = extraction_end - extraction_start + 1;
        
        extractIdx = [];
        for iTrial = 1:num_trials
            start_idx = extraction_start + (iTrial - 1) * trial_length;
            end_idx = extraction_end + (iTrial - 1) * trial_length;
            extractIdx = [extractIdx start_idx:end_idx];
        end
        
        legendEntry = {'Projected Intensity Function'};
        
        % **Parallelize over kernel_width_pool**
        parfor idx = 1:length(kernel_width_pool)
            kernel_width = kernel_width_pool(idx);
                        
            % Smoothing
            smoothed_spike_trains1 = gaussian_kernel_smoothing(concatenated_spike_trains1, kernel_width);
            smoothed_permuted_spike_trains1 = gaussian_kernel_smoothing(concatenated_permuted_spike_trains1, kernel_width);
            smoothed_spike_trains1 = smoothed_spike_trains1(:, extractIdx);
            smoothed_permuted_spike_trains1 = smoothed_permuted_spike_trains1(:, extractIdx);
            smoothed_spike_trains2 = gaussian_kernel_smoothing(concatenated_spike_trains2, kernel_width);
            smoothed_permuted_spike_trains2 = gaussian_kernel_smoothing(concatenated_permuted_spike_trains2, kernel_width);
            smoothed_spike_trains2 = smoothed_spike_trains2(:, extractIdx);
            smoothed_permuted_spike_trains2 = smoothed_permuted_spike_trains2(:, extractIdx);
            
            % PCA
            [coeffs1, scores1, explained1] = pca(smoothed_spike_trains1');
            [coeffs_perm1, scores_perm1, explained_perm1] = pca(smoothed_permuted_spike_trains1');
            [coeffs2, scores2, explained2] = pca(smoothed_spike_trains2');
            [coeffs_perm2, scores_perm2, explained_perm2] = pca(smoothed_permuted_spike_trains2');
            
            % CCA
            index1 = find(cumsum(explained1/(sum(explained1))) >= percentVar, 1);
            index2 = find(cumsum(explained2/(sum(explained2))) >= percentVar, 1);
            index_perm1 = find(cumsum(explained_perm1/(sum(explained_perm1))) >= percentVar, 1);
            index_perm2 = find(cumsum(explained_perm2/(sum(explained_perm2))) >= percentVar, 1);
            [A, B, r] = canoncorr(scores1(:,1:index1), scores2(:,1:index2));
            [A_perm, B_perm, r_perm] = canoncorr(scores_perm1(:,1:index_perm1), scores_perm2(:,1:index_perm2));
            
            % Continuum Regression
            [A, B, U, V, P, Q, ~, objComponent, ~, varExplained] = ...
                continuumRegressionMulti(smoothed_spike_trains1', smoothed_spike_trains2', ...
                alpha, nComponents, tol, maxIter, percentVar);
            [A_perm, B_perm, U_perm, V_perm, P_perm, Q_perm, ~, objComponent_perm, ~, varExplained_perm] = ...
                continuumRegressionMulti(smoothed_permuted_spike_trains1', smoothed_permuted_spike_trains2', ...
                alpha, nComponents, tol, maxIter, percentVar);
            
            R_total = varExplained * diag(corr(U, V));
            R_total_perm = varExplained_perm * diag(corr(U_perm, V_perm));
            
            % Error and correlation computations
            total_firing_rates1 = repmat(firing_rates1(:, :, 1), 1, num_trials);
            total_firing_rates2 = repmat(firing_rates2(:, :, 1), 1, num_trials);
            intensityAfterProjection1 = (total_firing_rates1(:, extractIdx) - ...
                mean(total_firing_rates1(:, extractIdx), 2))' * coeffs1;
            intensityAfterProjection2 = (total_firing_rates2(:, extractIdx) - ...
                mean(total_firing_rates2(:, extractIdx), 2))' * coeffs2;
            
            error1 = [rmse(intensityAfterProjection1(:,1), scores1(:,1)), ...
                      rmse(intensityAfterProjection1(:,2), scores1(:,2))];
            error2 = [rmse(intensityAfterProjection2(:,1), scores2(:,1)), ...
                      rmse(intensityAfterProjection2(:,2), scores2(:,2))];
            corr1 = [abs(corr(intensityAfterProjection1(:,1), scores1(:,1))), ...
                     abs(corr(intensityAfterProjection1(:,2), scores1(:,2)))];
            corr2 = [abs(corr(intensityAfterProjection2(:,1), scores2(:,1))), ...
                     abs(corr(intensityAfterProjection2(:,2), scores2(:,2)))];
            
            % **Store results directly in preallocated arrays**
            totalIndex1(idx, iRepeat) = index1;
            totalIndex2(idx, iRepeat) = index2;
            totalIndex_perm1(idx, iRepeat) = index_perm1;
            totalIndex_perm2(idx, iRepeat) = index_perm2;
            totalError1_repeat(iRepeat, idx, :) = error1;
            totalError2_repeat(iRepeat, idx, :) = error2;
            totalExplained1_repeat(iRepeat, idx, :) = explained1(1:2)' / sum(explained1);
            totalExplained2_repeat(iRepeat, idx, :) = explained2(1:2)' / sum(explained2);
            totalCorr1_repeat(iRepeat, idx, :) = corr1;
            totalCorr2_repeat(iRepeat, idx, :) = corr2;
            totalCorrCCA_repeat(iRepeat, idx) = mean(r);
            totalCorrCCA_perm_repeat(iRepeat, idx) = mean(r_perm);
            totalRA_repeat(iRepeat, idx) = R_total;
            totalRA_perm_repeat(iRepeat, idx) = R_total_perm;
            objComponent_repeat{iRepeat, idx} = objComponent;
            objComponent_repeat_perm{iRepeat, idx} = objComponent_perm;
            A_total_repeat{iRepeat, idx} = A;
            B_total_repeat{iRepeat, idx} = B;
            U_total_repeat{iRepeat, idx} = U;
            V_total_repeat{iRepeat, idx} = V;
       end
             
       disp(iRepeat)
    end
    
    elapsedTime = toc;
    display(elapsedTime);
    save(['Results\Loop_for_kernel_at_alpha_' num2str(alpha,'%.2f') '.mat']);
end