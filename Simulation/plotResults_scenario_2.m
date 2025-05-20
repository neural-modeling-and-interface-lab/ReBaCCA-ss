clc,clear,close all

load('Results/Loop_for_kernel_at_alpha_0.50_scenario 2.mat')

%% Calculate the real value
% figure
% blue_key = [0.6, 0.7, 1;    % Light blue
%             0.6, 0.6, 1;    % Medium blue
%             0, 0, 0.5];     % Dark blue
% pos = [0, 0.5, 1];         % Positions for interpolation (start, middle, end)
% colorsX = interp1(pos, blue_key, linspace(0, 1, length(kernel_width_pool)));
% orange_key = [1, 0.7, 0.5;  % Light orange
%               1, 0.7, 0.3;  % Medium orange
%               1, 0.3, 0];   % Dark orange
% pos = [0, 0.5, 1];         % Positions for interpolation
% colorsY = interp1(pos, orange_key, linspace(0, 1, length(kernel_width_pool)));
% 
% for iRepeat = 1:totalRepeat
%     for iKernel = 1:length(kernel_width_pool)
%         objComp = objComponent_repeat{iRepeat, iKernel};
%         varEx = objComp(1,:).^0.5 .* objComp(3,:).^0.5;
%         cutoffIdx = find(cumsum(varEx)/sum(varEx) > 0.8, 1);
%         cutoffIdx = 20;
%         totalCutOff(iKernel) = cutoffIdx;
%         test_repeat(iRepeat, iKernel) = sum(prod(objComp(:, 1:cutoffIdx).^0.5));
%         objComp_perm = objComponent_repeat_perm{iRepeat, iKernel}(:, 1:cutoffIdx);
%         varEx_perm = objComp_perm(1,:).^0.5 .* objComp_perm(3,:).^0.5;
%         test_repeat_perm(iRepeat, iKernel) = sum(prod(objComp_perm.^0.5));
% 
% 
%         plot(cumsum(varEx), 'Color', colorsX(iKernel, :), 'LineWidth', 1);
%         hold on
%         plot(cumsum(varEx_perm), 'Color', colorsY(iKernel, :), 'LineWidth', 1)
%     end
% end
% figure
% plot(kernel_width_pool,mean(test_repeat,1),'-','Color',[0 0.4470 0.7410],'LineWidth', 2)
% hold on
% plot(kernel_width_pool,mean(test_repeat_perm,1),'--','Color',[0 0.4470 0.7410],'LineWidth', 2)
% plot(kernel_width_pool,mean(test_repeat - test_repeat_perm,1),'-r','LineWidth', 2)
% legend('ReBaCCA','ReBaCCA-perm','ReBaCCA-ss','Location','best')
% xlabel('Kernel width (ms)')
% 
% figure;
% plot(totalCutOff,'-o');


%% Plot firing rate and spike raster
figure
subplot(221)
plot((extraction_start:extraction_end)-extraction_start,firing_rates1(1,extraction_start:extraction_end, 1),'b','LineWidth',2);
hold on
plot((extraction_start:extraction_end)-extraction_start,firing_rates1(end,extraction_start:extraction_end, 1),'r','LineWidth',2);
xlim([0 extraction_end-extraction_start])
xlabel('');xticklabels('');
title('Firing probability')
legend('Type 1','Type 2','Location','best')

subplot(222)
plot((extraction_start:extraction_end)-extraction_start,firing_rates2(1,extraction_start:extraction_end, 1),'b','LineWidth',2);
hold on
plot((extraction_start:extraction_end)-extraction_start,firing_rates2(end,extraction_start:extraction_end, 1),'r','LineWidth',2);
xlim([0 extraction_end-extraction_start])
xlabel('');xticklabels('');
title('Firing probability')

subplot(2,2,3);
plotSpikeRaster(spike_trains1(:, :, 2),num_neurons1,trial_length/1000,...
    'Spike Raster',split_index1)
xlabel('Time (ms)')
xlim([extraction_start extraction_end]/1000)
ax=gca;
xt=ax.XTick;
new_xt = xt*1000-extraction_start;
ax.XTickLabel = arrayfun(@num2str, new_xt, 'UniformOutput', false);

subplot(2, 2, 4);
plotSpikeRaster(spike_trains2(:, :, 2),num_neurons2,trial_length/1000,...
    'Spike Raster',split_index2)
xlabel('Time (ms)')
xlim([extraction_start extraction_end]/1000)
ax=gca;
xt=ax.XTick;
new_xt = xt*1000-extraction_start;
ax.XTickLabel = arrayfun(@num2str, new_xt, 'UniformOutput', false);

% Copy an existing figure only leave content
oldFig = figure(1); % Assuming your figure is figure 1
newFig = figure; % Create a new figure
% Find all axes in the old figure and copy them
oldAxes = findall(oldFig, 'type', 'axes'); % Find axes in the old figure
newAxes = copyobj(oldAxes, newFig); % Copy axes and their content to the new figure
% Loop over the copied axes to modify them
% Find and store the legend, if it exists
for ax = newAxes'  
    xlabel(ax, ''); % Remove X-axis label
    ylabel(ax, ''); % Remove Y-axis label
    title(ax, '');  % Remove title
    set(ax, 'XTick', [], 'YTick', []); % Remove tick marks
    box(ax, 'off'); % Optionally remove the box around the plot
end

% Copy an existing figure only leave axis
newFig = figure; % Create a new figure
% Find all axes in the old figure and copy them
oldLeg = findobj(oldFig, 'Type', 'Legend');
oldAxes = findall(oldFig, 'type', 'axes'); % Find axes in the old figure
newAxes = copyobj([oldAxes;oldLeg], newFig); % Copy axes and their content to the new figure
for ax = newAxes(1:2)
    % Remove plot content but keep the axes
    delete(findall(ax, 'type', 'line'));   % Delete all line plots
    delete(findall(ax, 'type', 'patch'));  % Delete patches (e.g., bar plots, fill plots)
    delete(findall(ax, 'type', 'surface')); % Delete surfaces (e.g., 3D plots)
    delete(findall(ax, 'type', 'image'));  % Delete images (if applicable)
end

kernel_width = 50;
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

figure
h5=subplot(3, 2,1);
avgPlot=reshape(scores1,[extraction_trial_length size(scores1,1)/extraction_trial_length size(scores1,2) ]);
avgPlot=squeeze(mean(avgPlot,2));
plot_pca(extraction_start:extraction_end,avgPlot, explained1, 'PCA');
xticklabels('')
title('Dataset 1')
xlim([extraction_start extraction_end])
ylabel('PCA')

h6=subplot(3, 2, 2);
avgPlot=reshape(scores2,[extraction_trial_length size(scores2,1)/extraction_trial_length size(scores2,2) ]);
avgPlot=squeeze(mean(avgPlot,2));
plot_pca(extraction_start:extraction_end,avgPlot, explained2, 'PCA');
xticklabels('')
title('Dataset 2')
legend('PC1','PC2','Location','best')
xlim([extraction_start extraction_end])

% CCA
index1 = find(cumsum(explained1/(sum(explained1))) >= percentVar, 1);
index2 = find(cumsum(explained2/(sum(explained2))) >= percentVar, 1);
index_perm1 = find(cumsum(explained_perm1/(sum(explained_perm1))) >= percentVar, 1);
index_perm2 = find(cumsum(explained_perm2/(sum(explained_perm2))) >= percentVar, 1);

% [A, B, r, U, V] = canoncorr(scores1(:,1:index1), scores2(:,1:index2));
[A, B, r, U, V] = canoncorr(smoothed_spike_trains1', smoothed_spike_trains2');
avgPlotU=reshape(U,[extraction_trial_length size(U,1)/extraction_trial_length size(U,2) ]);
avgPlotU=squeeze(mean(avgPlotU,2));
avgPlotV=reshape(V,[extraction_trial_length size(V,1)/extraction_trial_length size(V,2) ]);
avgPlotV=squeeze(mean(avgPlotV,2));

% S1 = scores1(:, 1:index1);
% S2 = scores2(:, 1:index2);
S1 = smoothed_spike_trains1';
S2 = smoothed_spike_trains2';
u = S1 * A(:, 1) / norm(A(:, 1));
varEx = u' * u / trace(S1' * S1);
subplot(323)
plot((extraction_start:extraction_end)-extraction_start,avgPlotU(:,1),'LineWidth',2);
xticklabels('')
title([num2str(varEx*100,'%.0f') '% of total variance'])
ylabel('CCA')

u = S2 * B(:, 1) / norm(B(:, 1));
varEx = u' * u / trace(S2' * S2);
subplot(324)
plot((extraction_start:extraction_end)-extraction_start,avgPlotV(:,1),'LineWidth',2);
xticklabels('')
title([num2str(varEx*100,'%.0f') '% of total variance'])
legend('CCA1','Location','best')

u = S1 * A(:, 2) / norm(A(:, 2));
varEx = u' * u / trace(S1' * S1);
subplot(325)
plot((extraction_start:extraction_end)-extraction_start,avgPlotU(:,2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
xlabel('Time (ms)');
title([num2str(varEx*100,'%.0f') '% of total variance'])
ylabel('CCA')

u = S2 * B(:, 2) / norm(B(:, 2));
varEx = u' * u / trace(S2' * S2);
subplot(326)
plot((extraction_start:extraction_end)-extraction_start,avgPlotV(:,2),'LineWidth',2,'Color',[0.8500 0.3250 0.0980]);
xlabel('Time (ms)');
title([num2str(varEx*100,'%.0f') '% of total variance'])
legend('CCA2','Location','best')

linkaxes([h5 h6],'y')



%% Plot the ReBaCCA, CCA, and permuted ReBaCCA results
x = kernel_width_pool;
y=mean(totalRA_repeat-totalRA_perm_repeat,1);
% y = mean(test_repeat - test_repeat_perm,1);
[maxValue, maxIndex] = max(y);
maxX = x(maxIndex);
maxY = y(maxIndex);

figure
plot(kernel_width_pool,mean(totalRA_repeat,1),'-','Color',[0 0.4470 0.7410],'LineWidth', 2)
hold on
plot(kernel_width_pool,mean(totalRA_perm_repeat,1),'-.','Color',[0 0.4470 0.7410],'LineWidth', 2)
line([maxX maxX], [0 1], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
% text(maxX, 0, sprintf('%.1f', maxX),'Color','r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
xlabel('Kernel width (ms)')
ylim([0 1])
ylabel('Similarity')
legend('$V_{\alpha}(S_1,S_2)$','$V_{\alpha}(\tilde{S}_1,\tilde{S}_2)$','Location','best','Interpreter','latex')
% set(gca,'XScale','log')
xlim([0 100])

% Plot the objective function (difference between original and permuted ReBaCCA)
figure
xlabel('Kernel width (ms)')
title('ReBaCCA-ss')
hold on
% Find the maximum value and its index
plot(kernel_width_pool,mean(totalRA_repeat-totalRA_perm_repeat,1),'-k','Linewidth',2)
plot(maxX, maxY, 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
line([maxX maxX], [0 maxY], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
% text(maxX, 0, sprintf('%.1f', maxX),'Color','r', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
% set(gca,'XScale','log')
xlim([0 100])
ylim([0 1])

% Plot the number of dimensions that capture 80% variance for original and permuted data
figure
plot(kernel_width_pool,mean(totalIndex1,2),'-ro')
hold on
plot(kernel_width_pool,mean(totalIndex_perm1,2),'-bo')
xlabel('Kernel width (ms)')
title('Number of dimension that contains 80% of variance')
legend('Original spike data','Surrogated data','Location','Best')


% Plot the correlation between reconstructed latent dynamics and ground truth
figure
errorToPlot1=squeeze(mean(totalCorr1_repeat,1));
errorToPlot2=squeeze(mean(totalCorr2_repeat,1));
y=(sum(errorToPlot1,2)+sum(errorToPlot2,2))/4;
plot(kernel_width_pool,y,'-k','LineWidth',2);
xlabel('Kernel width (ms)')
title('Correlation between reconstruction and ground truth')
hold on
% Plot the maximum point
plot(maxX, y(maxIndex), 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
line([maxX maxX], [0 y(maxIndex)], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
% text(maxX, 0, sprintf('%.1f', maxX),'Color','k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
% set(gca,'XScale','log')
xlim([0 100])

% Plot the averaged correlation from CCA for original and permuted data
figure
plot(kernel_width_pool,mean(totalCorrCCA_repeat,1),'-ro')
hold on
plot(kernel_width_pool,mean(totalCorrCCA_perm_repeat,1),'-bo')
xlabel('Kernel width (ms)')
title('Averaged correlation from CCA')
ylim([0 1])
xlim([0 100])
% set(gca,'XScale','log')

% Plot the root mean squared error (RMSE) between latent dynamics and ground truth
figure
errorToPlot1=squeeze(mean(totalError1_repeat,1));
errorToPlot2=squeeze(mean(totalError2_repeat,1));
y=(sum(errorToPlot1,2)+sum(errorToPlot2,2))/4;
plot(kernel_width_pool,y,'-k','LineWidth',2)
% plot(kernel_width_pool,errorToPlot2(:,1),'-ro')
xlabel('Kernel width (ms)')
title('Root mean squared error')
hold on
% Plot the maximum point
plot(maxX, y(maxIndex), 'ro', 'MarkerSize', 10, 'MarkerFaceColor','r');
line([maxX maxX], [0 y(maxIndex)], 'Color', 'r', 'LineStyle', '--', 'LineWidth', 1);
xlim([0 100])
% text(maxX, 0, sprintf('%.1f', maxX),'Color','k', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
% set(gca,'XScale','log')

% Plot the reconstruction results
% figure
% avgPlot=reshape(intensityAfterProjection1,[extraction_trial_length ...
%     size(intensityAfterProjection1,1)/extraction_trial_length size(intensityAfterProjection1,2) ]);
% avgPlot=squeeze(mean(avgPlot,2));
% plot(extraction_start:extraction_end, avgPlot(:,1),'b','Linewidth',2);
% hold on;
% plot(extraction_start:extraction_end, avgPlot(:,2),'r','Linewidth',2);
% xlabel('');xticklabels('');
% title('Firing Probability Projected on PCA - Dataset 1, Trial 2');
% xlim([0 trial_length])
% 
% h2=subplot(2,2,2);
% avgPlot=reshape(intensityAfterProjection2,[extraction_trial_length ...
%     size(intensityAfterProjection2,1)/extraction_trial_length size(intensityAfterProjection2,2) ]);
% avgPlot=squeeze(mean(avgPlot,2));
% plot(extraction_start:extraction_end, avgPlot(:,1),'b','Linewidth',2);
% hold on;
% plot(extraction_start:extraction_end, avgPlot(:,2),'r','Linewidth',2);
% xlabel('');xticklabels('');
% title('Firing Probability Projected on PCA - Dataset 2, Trial 2');
% xlim([0 trial_length])  

%% Output the figure
export_pdf_figure(2,'Fig/Firing and raster only content scenario 2',0)
export_pdf_figure(3,'Fig/Firing and raster only axis scenario 2',0)
export_pdf_figure(4,'Fig/PCA and CCA scenario 2',0)
export_pdf_figure(5,'Fig/ReBaCCA scenario 2',0)
export_pdf_figure(6,'Fig/ReBaCCA-ss scenario 2',0)
export_pdf_figure(8,'Fig/Correlation scenario 2',0)
export_pdf_figure(10,'Fig/RMSE scenario 2',0)
