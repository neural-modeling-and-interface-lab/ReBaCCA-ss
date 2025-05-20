clc, clear, close all

rng(42)

% Load for across session
dataNameTotal = {{'1053_6','1053_1'},{'1053_6','1029_1'},{'1053_1','1029_1'}};
eventName = {'LN', 'LS', 'RN', 'RS'};
num_selected = 8;
D = ones(num_selected * length(dataNameTotal));
subSize = num_selected * length(eventName);

pairs = [2 1; 3 1; 3 2];
for iData = 1:length(dataNameTotal)
	dataNamePair = dataNameTotal{iData};
	load(['Results/Cross Session ReBaCCA for MDS ' dataNamePair{1} ' vs ' dataNamePair{2}]);
	rows = (pairs(iData,1) - 1) * subSize + 1 : pairs(iData,1) * subSize;
	cols = (pairs(iData,2) - 1) * subSize + 1 : pairs(iData,2) * subSize;
	D(rows, cols) = ReBaCCATotal;
end

dataNameTotal = {'1053_6', '1053_1', '1029_1'};
% Load for within session
for iData = 1:length(dataNameTotal)
	dataName = dataNameTotal(iData);
	load(['Results/Self ReBaCCA for MDS ' dataName{:}]);
	rows = (iData - 1) * subSize + 1 : iData * subSize;
	cols = (iData - 1) * subSize + 1 : iData * subSize;
	D(rows, cols) = ReBaCCATotal;
end

sideLen = 0.02;
eventLabel = repmat(repelem([1, 2, 3, 4], 8), 1, 3);
animalLabel = [ones(1, 64) 2*ones(1,32)];
eventTriangle = {[-1.5*sideLen 0 sideLen; -1.5*sideLen sideLen 0], [-1.5*sideLen 0 sideLen;1.5*sideLen -sideLen 0],...
    [1.5*sideLen -sideLen 0; -1.5*sideLen 0 sideLen], [1.5*sideLen -sideLen 0; 1.5*sideLen 0 -sideLen;]};

datasetColors = [
    [1 0 0];    % 'LEFT nonmatch' - dark red
    [1 0 1];  % 'RIGHT nonmatch' - light red (salmon)
    [0.3, 0.3, 1];  % 'RIGHT sample' - light blue
];
% Define dataset labels
datasetLabel = repelem(1:3, 32); % 1 repeated 32 times, 2 repeated 32 times, 3 repeated 32 times

S = tril(D) + triu(D',1);

% Step 2: Apply classical MDS
[Y, eigvals] = cmdscale(1 - S);

% Just for visualization, rotate and flip
% Y = Y(:,1:2) * [0 -1; 1 0];

% Step 3: Plot the MDS result
figure;
hold on;
tempY = Y(:,1:2);
% Step 3: Perform k-means clustering
[idx, C] = kmeans(tempY, 3, 'Start', 'plus', 'Replicates', 10);

% Step 4: Plot ellipses for each cluster
for k = 1:3
    % Extract points in the current cluster
    Y_cluster = tempY(idx == k, :);
    % Ensure there are enough points to compute covariance
    if size(Y_cluster, 1) > 1
        % Compute covariance matrix
        Sigma = 3 * cov(Y_cluster);
        % Get eigenvectors (V) and eigenvalues (D)
        [V, D] = eig(Sigma);
        % Generate points on a unit circle
        theta = linspace(0, 2*pi, 100);
        U = [cos(theta); sin(theta)];
        % Scale by standard deviations (square roots of eigenvalues)
        sqrtD = diag(sqrt(diag(D)));
        % Compute ellipse points
        ellipse_points = V * sqrtD * U;
        % Plot the ellipse
        plot(C(k,1) + ellipse_points(1,:), C(k,2) + ellipse_points(2,:), '--k', 'LineWidth', 0.5);
    end

end

for eventIdx = 1:4
    for datasetIdx = 1:3
        % Find indices for this event and dataset
        indices = find(eventLabel == eventIdx & datasetLabel == datasetIdx);
        for iPoint = 1:length(indices)
            xVerts = Y(indices(iPoint),1) + eventTriangle{eventIdx}(1,:);
            yVerts = Y(indices(iPoint),2) + eventTriangle{eventIdx}(2,:);
            if eventIdx == 1 || eventIdx == 3
                fill(xVerts, yVerts, datasetColors(datasetIdx,:),'Linewidth', 1.5,...
                    'EdgeColor',datasetColors(datasetIdx,:),'FaceAlpha',1);
            else
                fill(xVerts, yVerts, datasetColors(datasetIdx,:), 'Linewidth', 1.5, ...
                    'FaceColor','white','FaceAlpha',1, 'EdgeColor', datasetColors(datasetIdx,:));
            end
        end
        % % Plot with event color and dataset marker; set DisplayName only for first dataset
        % scatter(Y(indices,1), Y(indices,2), 100, datasetColors(eventIdx,:), markers{datasetIdx}, 'filled');
    end
end

axis equal
xlim([-0.33 0.33])
xticks([-0.2 0 0.2])
ylim([-0.33 0.33])
yticks([-0.2 0 0.2])
xlabel('MDS 1');ylabel('MDS 2')
title('MDS for 3 datasets', 'Interpreter','none');
hold off;

export_pdf_figure(gcf, 'Fig/MDS for 3 datasets', 0)

%% plot for same event different sessions
for eventIdx = 1:4
    figure;
    hold on;
    tempY = Y(eventLabel == eventIdx,1:2);
    % Step 3: Perform k-means clustering
    [idx, C] = kmeans(tempY, 3, 'Start', 'plus', 'Replicates', 10);

    % Step 4: Plot ellipses for each cluster
    for k = 1:3
        % Extract points in the current cluster
        Y_cluster = tempY(idx == k, :);
        % Ensure there are enough points to compute covariance
        if size(Y_cluster, 1) > 1
            % Compute covariance matrix
            Sigma = 3 * cov(Y_cluster);
            % Get eigenvectors (V) and eigenvalues (D)
            [V, D] = eig(Sigma);
            % Generate points on a unit circle
            theta = linspace(0, 2*pi, 100);
            U = [cos(theta); sin(theta)];
            % Scale by standard deviations (square roots of eigenvalues)
            sqrtD = diag(sqrt(diag(D)));
            % Compute ellipse points
            ellipse_points = V * sqrtD * U;
            % Plot the ellipse
            plot(C(k,1) + ellipse_points(1,:), C(k,2) + ellipse_points(2,:), '--k', 'LineWidth', 0.5);
        end

    end

    for datasetIdx = 1:3
        % Find indices for this event and dataset
        indices = find(eventLabel == eventIdx & datasetLabel == datasetIdx);
        for iPoint = 1:length(indices)
            xVerts = Y(indices(iPoint),1) + eventTriangle{eventIdx}(1,:);
            yVerts = Y(indices(iPoint),2) + eventTriangle{eventIdx}(2,:);
            if eventIdx == 1 || eventIdx == 3
                patch(xVerts, yVerts, datasetColors(datasetIdx,:), 'Linewidth', 1.5,...
                    'EdgeColor','none','FaceAlpha',1);
                % plot([xVerts xVerts(1)], [yVerts yVerts(1)], 'Linewidth', 1.5,...
                %     'Color',datasetColors(datasetIdx,:),'LineJoin','round');
            else
                patch(xVerts, yVerts, datasetColors(datasetIdx,:),'Linewidth', 1.5,...
                    'FaceColor','white','FaceAlpha',1, 'EdgeColor', datasetColors(datasetIdx,:));
            end
        end
    end
    
    
    axis equal
    xlim([-0.35 0.35])
    xticks([-0.2 0 0.2])
    ylim([-0.35 0.35])
    yticks([-0.2 0 0.2])
    title(eventName{eventIdx})
    box on
    xlabel('MDS 1');ylabel('MDS 2')

    export_pdf_figure(gcf, ['Fig/MDS for 3 datasets for event ' eventName{eventIdx}], 0)

    %% Plot pair-wise evaluation
    pairSim = S;
    within_session_dis = [];
    for datasetIdx = 1:3
        matrixPos = find(datasetLabel == datasetIdx & eventLabel == eventIdx);
        for i = 1:length(matrixPos)
            for j = i+1:length(matrixPos)
                within_session_dis = [within_session_dis pairSim(matrixPos(i), matrixPos(j))];
            end
        end
    end

    cross_session_dis = [];
    for datasetIdx = 1:3
        matrixPos1 = find(datasetLabel == datasetIdx & eventLabel == eventIdx);
        matrixPos2 = find(datasetLabel > datasetIdx & eventLabel == eventIdx);
        for i = 1:length(matrixPos1)
            for j = 1:length(matrixPos2)
               cross_session_dis = [cross_session_dis pairSim(matrixPos1(i), matrixPos2(j))];
            end
        end
    end

    % Calculate means and SEMs for each category
    mean_within = mean(within_session_dis);
    std_within = std(within_session_dis);

    mean_diff = mean(cross_session_dis);
    std_diff = std(cross_session_dis);

    % Prepare data for plotting
    means = [mean_within, mean_diff];
    stds = [std_within, std_diff];
    
    % Perform t-test
    [h, p] = ttest2(within_session_dis, cross_session_dis);
    
    % Plotting
    figure;
    set(gcf, 'Position', [100, 100, 600, 400]);
    hold on;
    bar(1:2, means, 0.6, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
    errorbar(1:2, means, stds, 'color', 'k', 'LineWidth', 1, 'CapSize', 10, 'LineStyle', 'none');
    set(gca, 'FontName', 'Arial', 'FontSize', 12);
    title('Pairwise Similarity', 'FontName', 'Arial', 'FontSize', 14);
    
    % Add line and star if significant
    if p < 0.05
        % Calculate the y-position for the line (10% above the highest bar with error bar)
        max_y = max(means + stds) * 1.05;
        
        % Draw the line between the bars
        plot([1, 2], [max_y, max_y], 'k-', 'LineWidth', 1.5);
        
        % Add a star above the line
        text(1.5, max_y + 0.002, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    hold off;
    ylim([0 0.14]);
    yticks([0 0.05 0.1]);
    xticks([1 2]);
    xticklabels({'Within\newlinesession', 'Cross\newlinesessions'});
    pbaspect([1 2 1]);
    export_pdf_figure(gcf, ['Fig/Similarity measure session comparison ' eventName{eventIdx}], 0);


    %% Cross animal comparison
    within_animal_dis = [];
    for animalIdx = 1:2
        matrixPos = find(animalLabel == animalIdx & eventLabel == eventIdx);
        for i = 1:length(matrixPos)
            for j = i+1:length(matrixPos)
                within_animal_dis = [within_animal_dis pairSim(matrixPos(i), matrixPos(j))];
            end
        end
    end

    cross_animal_dis = [];
    for animalIdx = 1:2
        matrixPos1 = find(animalLabel == animalIdx & eventLabel == eventIdx);
        matrixPos2 = find(animalLabel > animalIdx & eventLabel == eventIdx);
        for i = 1:length(matrixPos1)
            for j = 1:length(matrixPos2)
               cross_animal_dis = [cross_animal_dis pairSim(matrixPos1(i), matrixPos2(j))];
            end
        end
    end

    % Calculate means and SEMs for each category
    mean_within = mean(within_animal_dis);
    std_within = std(within_animal_dis);

    mean_diff = mean(cross_animal_dis);
    std_diff = std(cross_animal_dis);

    % Prepare data for plotting
    means = [mean_within, mean_diff];
    stds = [std_within, std_diff];
    
    % Perform t-test
    [h, p] = ttest2(within_animal_dis, cross_animal_dis);
    
    % Plotting
    figure;
    set(gcf, 'Position', [100, 100, 600, 400]);
    hold on;
    bar(1:2, means, 0.6, 'FaceColor', [0.7, 0.7, 0.7], 'EdgeColor', 'none');
    errorbar(1:2, means, stds, 'color', 'k', 'LineWidth', 1, 'CapSize', 10, 'LineStyle', 'none');
    set(gca, 'FontName', 'Arial', 'FontSize', 12);
    title('Pairwise Similarity', 'FontName', 'Arial', 'FontSize', 14);
    
    % Add line and star if significant
    if p < 0.05
        % Calculate the y-position for the line (10% above the highest bar with error bar)
        max_y = max(means + stds) * 1.05;
        
        % Draw the line between the bars
        plot([1, 2], [max_y, max_y], 'k-', 'LineWidth', 1.5);
        
        % Add a star above the line
        text(1.5, max_y + 0.002, '*', 'HorizontalAlignment', 'center', 'FontSize', 16, 'FontWeight', 'bold');
    end
    
    hold off;
    ylim([0 0.14]);
    yticks([0 0.05 0.1]);
    xticks([1 2]);
    xticklabels({'Within\newlineanimal', 'Cross\newlineanimals'});
    pbaspect([1 2 1]); % Set the plot box aspect ratio so x-axis is twice as long as y-axis
    export_pdf_figure(gcf, ['Fig/Similarity measure animal comparison ' eventName{eventIdx}], 0);

end
