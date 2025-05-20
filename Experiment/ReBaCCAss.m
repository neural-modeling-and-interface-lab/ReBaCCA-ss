function [optimal_kernel, optimal_value, optimal_VAE, optimal_Corr, optimal_VAE_perm, optimal_Corr_perm] = ReBaCCAss(...
    spike_trains1, spike_trains2, kernel_width_pool, alpha, tol, maxIter, percentVar, extractIdx)

totalPermuteTime = 4;
totalRA_kernel = zeros(length(kernel_width_pool), 1);
totalRA_perm_kernel = zeros(length(kernel_width_pool), 1);
totalDimCorr = cell(length(kernel_width_pool), 1);
totalDimVAE = cell(length(kernel_width_pool), 1);
totalDimCorr_perm = cell(length(kernel_width_pool), 1);
totalDimVAE_perm = cell(length(kernel_width_pool), 1);

parfor idx = 1:length(kernel_width_pool)

    kernel_width = kernel_width_pool(idx);

    smoothed_spike_trains1 = gaussian_kernel_smoothing(spike_trains1, kernel_width);
    smoothed_spike_trains2 = gaussian_kernel_smoothing(spike_trains2, kernel_width);

    smoothed_spike_trains1 = smoothed_spike_trains1(:, extractIdx);
    smoothed_spike_trains2 = smoothed_spike_trains2(:, extractIdx);

    nComponents = min([rank(spike_trains1), rank(spike_trains2)]);

    [A, B, U, V, P, Q, ~, objComponent, ~, varExplained,varExplainedX,varExplainedY] = ...
        continuumRegressionMulti(smoothed_spike_trains1', smoothed_spike_trains2', ...
        alpha, nComponents, tol, maxIter, percentVar);
    R_total = varExplained * diag(corr(U, V));
    totalDimCorr{idx} = diag(corr(U,V));
    totalDimVAE{idx} = [varExplainedX ; varExplainedY];

    R_total_perm = zeros(totalPermuteTime, 1);
    % Just for recording because sometimes after permutation, the rank is
    % different
    nComponents_rec = nComponents - 5;
    corr_total = zeros(nComponents_rec, totalPermuteTime);
    VAE_X_total = zeros(totalPermuteTime, nComponents_rec);
    VAE_Y_total = zeros(totalPermuteTime, nComponents_rec);
    for iPermute = 1:totalPermuteTime
        
        permuted_spike_trains1 = zeros(size(spike_trains1));
        for iNeuron = 1:size(spike_trains1, 1)
            permuted_spike_trains1(iNeuron, :) = spike_trains1(iNeuron, randperm(size(spike_trains1, 2)));
        end
        permuted_spike_trains2 = zeros(size(spike_trains2));
        for iNeuron = 1:size(spike_trains2, 1)
            permuted_spike_trains2(iNeuron, :) = spike_trains2(iNeuron, randperm(size(spike_trains2, 2)));
        end

        smoothed_permuted_spike_trains1 = gaussian_kernel_smoothing(permuted_spike_trains1, kernel_width);
        smoothed_permuted_spike_trains2 = gaussian_kernel_smoothing(permuted_spike_trains2, kernel_width);

        smoothed_permuted_spike_trains1 = smoothed_permuted_spike_trains1(:, extractIdx);
        smoothed_permuted_spike_trains2 = smoothed_permuted_spike_trains2(:, extractIdx);

        nComponents = min([rank(smoothed_permuted_spike_trains1), rank(smoothed_permuted_spike_trains2)]);

        [A_perm, B_perm, U_perm, V_perm, P_perm, Q_perm, ~, objComponent_perm, ~, varExplained_perm,...
            varExplainedX_perm,varExplainedY_perm] = ...
                continuumRegressionMulti(smoothed_permuted_spike_trains1', smoothed_permuted_spike_trains2', ...
                alpha, nComponents, tol, maxIter, percentVar);
        R_total_perm(iPermute) = varExplained_perm * diag(corr(U_perm, V_perm));
        
        tempCorr = diag(corr(U_perm,V_perm));
        corr_total(:, iPermute) = tempCorr(1:nComponents_rec);
        VAE_X_total(iPermute, :) = varExplainedX_perm(1:nComponents_rec);
        VAE_Y_total(iPermute, :) = varExplainedY_perm(1:nComponents_rec);
    end

    totalDimCorr_perm{idx} = mean(corr_total, 2);
    totalDimVAE_perm{idx} = [mean(VAE_X_total, 1) ; mean(VAE_Y_total, 1)];

    totalRA_kernel(idx) = R_total;
    totalRA_perm_kernel(idx) = mean(R_total_perm);
end

y = totalRA_kernel - totalRA_perm_kernel;
[~, maxIndex] = max(y);
optimal_kernel = kernel_width_pool(maxIndex);
optimal_value = y(maxIndex);
optimal_VAE = totalDimVAE{maxIndex};
optimal_Corr = totalDimCorr{maxIndex};
optimal_VAE_perm = totalDimVAE_perm{maxIndex};
optimal_Corr_perm = totalDimCorr_perm{maxIndex};

end