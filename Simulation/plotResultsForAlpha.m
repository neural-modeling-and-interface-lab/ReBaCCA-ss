clc,clear,close all

load('Results\Loop for alpha at kernel_45_alpha_0.00.mat')

% Plot the details of different alpha
% alphaPoolName = {'0.00', '0.50', '0.95'};
% kernelWidthIdx = 1;
% objVsAlpha = zeros(totalRepeat, length(alphaPoolName));
% alphaValues = zeros(1, length(alphaPoolName));
% obj_alpha = zeros(1, length(alphaPoolName));
% A = A_total_repeat{1};
% B = B_total_repeat{1};
% 
% % figure
% % set(gcf, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.4]); % Wide figure
% % for iAlpha = 1:length(alphaPoolName)
% %     load(['Results\Loop for alpha at kernel_15_alpha_' num2str(alphaPoolName{iAlpha}, '%.2f') '.mat'])
% % 
% % 
% %     alphaValues(iAlpha) = alpha;
% %     objVsAlpha(:,iAlpha) = totalRA_repeat(:, kernelWidthIdx);
% %     Atotal{iAlpha} = A;
% %     Btotal{iAlpha} = B;
% % 
% %     objComponent = objComponent_repeat{1};
% % 
% %     % ---------------------
% %     % Main subplot: original values
% %     % ---------------------
% %     mainIdx = iAlpha; % This is the main line plot index.
% %     subplot(1, 4, mainIdx);
% %     plot(objComponent(2,:), '-o')
% %     hold on
% %     plot(sqrt(objComponent(1,:).*objComponent(3,:)), '-o');
% %     title(['Alpha=' alphaPoolName{iAlpha}])
% %     if iAlpha == 1
% %         legend('Corr', 'Var Explained');
% %         xlabel('Dimension')
% %     end
% % 
% %     % Get position of the current subplot for reference
% %     mainPos = get(gca, 'Position');
% % 
% %     % ---------------------
% %     % Small bar plot (for objVsAlpha)
% %     % ---------------------
% %     % Instead of using subplot again, create a new axes with a smaller box
% %     barWidth = mainPos(3)*0.25;  % 25% of the main subplot width
% %     barHeight = mainPos(4)*0.5;  % 50% of the main subplot height
% %     % Position the bar to the right side of the main plot
% %     barLeft = mainPos(1) + mainPos(3)*0.7; % slightly to the right within main subplot
% %     barBottom = mainPos(2) + (mainPos(4)-barHeight)/2; % centered vertically
% % 
% %     barPlot = axes('Position', [barLeft, barBottom, barWidth, barHeight]);
% %     bar(1, mean(objVsAlpha(:,iAlpha)), 0.5, 'FaceColor', [0.5, 0.5, 0.5], 'EdgeColor', 'none');
% %     title('Overall');
% %     ylim([-0.01, 1]);
% %     xticks([]);
% %     xlim([0.5, 1.5]);
% %     set(barPlot, 'Box', 'off');
% % end
% % 
% % export_pdf_figure(gcf, 'Fig/Details of VACCC', 1);


%% Explore different alpha
alphaPool = [0:0.1:0.95 0.95];
alphaPoolName = arrayfun(@(x) sprintf('%.2f', x), alphaPool, 'UniformOutput', false);
kernelWidthIdx = 1;
objVsAlpha = zeros(totalRepeat, length(alphaPoolName));
alphaValues = zeros(1, length(alphaPoolName));
obj_alpha = zeros(1, length(alphaPoolName));

Afirst = zeros(size(A_total_repeat{1}, 1), length(alphaPoolName));
Bfirst = zeros(size(B_total_repeat{1}, 1), length(alphaPoolName));
objCompFirst = zeros(3, length(alphaPoolName));
for iAlpha = 1:length(alphaPoolName)
    load(['Results\Loop for alpha at kernel_45_alpha_' num2str(alphaPoolName{iAlpha}, '%.2f') '.mat'])

    tempObjComp = objComponent_repeat{1};
    objCompFirst(:, iAlpha) = tempObjComp(:, 1);
    
    A = A_total_repeat{1};
    B = B_total_repeat{1};    alphaValues(iAlpha) = alpha;
    objVsAlpha(:,iAlpha) = totalRA_repeat(:, kernelWidthIdx);
    Atotal{iAlpha} = A;
    Btotal{iAlpha} = B;
    Afirst(:, iAlpha) = A(:, 1);
    Bfirst(:, iAlpha) = B(:, 1);
end

%% For plotting
figure
h(1) = subplot(211);
% Optional: Set renderer to 'opengl' for better transparency handling (if needed)
set(gcf, 'Renderer', 'painters');
% Enable hold to overlay multiple curves
hold on;
blue_key = [0.6, 0.7, 1;    % Light blue
            0.6, 0.6, 1;    % Medium blue
            0, 0, 0.5];     % Dark blue
pos = [0, 0.5, 1];         % Positions for interpolation (start, middle, end)
colorsX = interp1(pos, blue_key, linspace(0, 1, length(alphaPool)));

% Plot each column with the same color but different transparency
for i = 1:size(Afirst, 2)
    plot(Afirst(:, i), 'Color', colorsX(i,:), 'linewidth', 2);
end
title('Pattern 1');
xticks([]);
colormap(h(1), colorsX);
cbar1 = colorbar;
cbar1.Title.String = '\alpha';

h(2) = subplot(212);
orange_key = [1, 0.7, 0.5;  % Light orange
              1, 0.7, 0.3;  % Medium orange
              1, 0.3, 0];   % Dark orange
pos = [0, 0.5, 1];         % Positions for interpolation
colorsY = interp1(pos, orange_key, linspace(0, 1, length(alphaPool)));
hold on
% Plot each column with the same color but different transparency
for i = 1:size(Bfirst, 2)
    plot(Bfirst(:, i), 'Color', colorsY(i,:), 'linewidth', 2);
end
colormap(h(2), colorsY);
cbar2 = colorbar;
xlabel('Vector index')
title('Pattern 2');
sgtitle('Projection vector')
ylim([-1 0.5])
linkaxes(h,'xy')
export_pdf_figure(gcf, 'Fig/Projection vector vs alpha', 0);


% % Projection direction differences
% pairDis = zeros(length(alphaPoolName));
% principalAngle = zeros(length(alphaPoolName));
% for iAlpha = 1:length(alphaPoolName)
%     for jAlpha = 1:length(alphaPoolName)
%         pairDis(iAlpha, jAlpha) = norm(Atotal{iAlpha}(:, 1:2) - Atotal{jAlpha}(:, 1:2), "fro") / ...
%             norm(Atotal{iAlpha}(:, 1:2), "fro");
%         % pairDis(iAlpha, jAlpha) = norm(Btotal{iAlpha}(:, 1:2) - Btotal{jAlpha}(:, 1:2), "fro") / ...
%         %     norm(Btotal{iAlpha}(:, 1:2), "fro");
%     end
% end
% [Y, eigenvals] = cmdscale(pairDis);
% 
% % Plot the projection space and VACCC values
% % Scatter plot with color based on normalized alphaPoolName
% figure;
% cmap = flipud(copper);
% colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, alphaPool);
% scatter(Y(:,1), Y(:,2), 100, colors, 'filled'); % Use colors for scatter plot
% colorbarHandle = colorbar; % Add a colorbar
% colormap(cmap)
% ylabel(colorbarHandle, 'Alpha Value', 'FontSize', 12);
% % colormap(gray); % Use grayscale colormap
% clim([0 1]); % Match color scale to alphaPoolName range
% xlabel('MDS Dim 1');
% ylabel('MDS Dim 2');
% title('Low dimensional representation of projection spaces');
% export_pdf_figure(gcf, 'Fig/Low Rep of projection', 0);

%% Obj versus alpha
figure
alphaValues(end) = 1;  % For visualization purpose, alpha=1 cannot calculate inverse
b = bar(alphaValues, objCompFirst([1 3], :));
hold on
% Apply the gradient to b(1)
set(b(1), 'FaceColor', 'flat', 'CData', colorsX, 'EdgeColor', 'none');
% Apply the gradient to b(2)
set(b(2), 'FaceColor', 'flat', 'CData', colorsY, 'EdgeColor', 'none');

% Add mean line and finalize plot
p = plot(alphaValues, mean(objVsAlpha), 'o-', 'LineWidth', 1.5, 'Color', [0.4940 0.1840 0.5560]);
legend([b(1), b(2), p], {'VAE 1', 'VAE 2', 'ReBaCCA'}, 'Location','best');
xlabel('\alpha');
ylim([0 1]); % Zooming in on the relevant y-axis range
xlim([-0.05 1.05]); % Keeping the original x-axis range
pbaspect([3 1 1]); % Set the plot box aspect ratio so x-axis is twice as long as y-axis
export_pdf_figure(gcf, 'Fig/ReBaCCA vs alpha curve', 0);
