clc,clear,close all

rng(42)
dataNamePairTotal = {{'1053_6', '1029_1'}, {'1053_6', '1053_1'}, {'1053_1', '1029_1'}};
for iDataFile = 1:length(dataNamePairTotal)
    dataNamePair = dataNamePairTotal{iDataFile};
    trialTypeName={'_LEFT_nonmatch','_LEFT_sample','_RIGHT_nonmatch','_RIGHT_sample'};
    
    spike_trains=cell(length(dataNamePair),length(trialTypeName));
    permuted_spike_trains=cell(length(dataNamePair),length(trialTypeName));
    
    selected_indices_together = cell(2, 4);       % Store indices for reference
    % Get all the spike train data
    for iTrialType=1:length(trialTypeName)
        concatenated_spike_trains=cell(size(dataNamePair));
        concatenated_permuted_spike_trains=cell(size(dataNamePair));
        for iData=1:length(dataNamePair)
            load(['Preprocessed data/' dataNamePair{iData} trialTypeName{iTrialType}])
            load(['Results\Self ReBaCCA for MDS ' dataNamePair{iData}])
            concatenated_spike_trains{iData}=success_spike_counts_concantenated;
            trial_length = timeAroundEvent * 1000;
            num_trials=size(success_spike_counts_concantenated,2)/trial_length;     
            
            for iTrial = 1:num_trials      
                spike_trains{iData, iTrialType}(iTrial, :, :) = success_spike_counts_concantenated(:,...
                    (iTrial - 1) * trial_length + 1 : iTrial * trial_length);             
            end
            selected_indices_together(iData, :) = selected_indices;
        end
    end
    
    %% Perform ReBaCCA-ss for same event same animal
    alpha = 0.5;
    tol = 1e-2;            % Tolerance for convergence
    maxIter = 64;         % Maximum iterations
    percentVar = 1-1e-2; 
    extractIdx = 500:5499;
    kernel_width_pool = logspace(0, 2, 16);
    num_selected = 8;
    
    tic;
    targetTrialType = [1 2 3 4];
    num_pairs = (4 * num_selected)^2;
    
    ReBaCCATotal_vec = zeros(num_pairs, 1);
    optimalKernelTotal_vec = zeros(num_pairs, 1);
    ReBaCCATotal = ones(4 * num_selected);
    optimalKernelTotal = zeros(4 * num_selected);
    % Generate all unique pairs (iTrial, jTrial)
    [pairs_i, pairs_j] = find(ones(4 * num_selected));
    
    parfor k = 1:length(pairs_i)
	    iTrial = pairs_i(k);
	    jTrial = pairs_j(k);
        trialType1 = ceil(iTrial / num_selected);
        trialType2 = ceil(jTrial / num_selected);
        trialIdx1 = mod(iTrial - 1, num_selected) + 1;
        trialIdx2 = mod(jTrial - 1, num_selected) + 1;
        tempIdx1 = selected_indices_together{1, trialType1}(trialIdx1);
        tempIdx2 = selected_indices_together{2, trialType2}(trialIdx2);
    
	    tempSpikeTrains1 = spike_trains{1, trialType1};
	    tempSpikeTrains2 = spike_trains{2, trialType2};
    
	    spike_trains1 = squeeze(tempSpikeTrains1(tempIdx1, :, :));
	    spike_trains2 = squeeze(tempSpikeTrains2(tempIdx2, :, :));
        
	    [optimal_kernel, optimal_value] = ReBaCCAss(...
		    spike_trains1, spike_trains2, kernel_width_pool, alpha, tol, maxIter, percentVar, extractIdx);
    
	    ReBaCCATotal_vec(k) = optimal_value;
	    optimalKernelTotal_vec(k) = optimal_kernel;
    end
    
    for k = 1:num_pairs
	    iTrial = pairs_i(k);
	    jTrial = pairs_j(k);
	    ReBaCCATotal(iTrial, jTrial) = ReBaCCATotal_vec(k);
	    optimalKernelTotal(iTrial, jTrial) = optimalKernelTotal_vec(k);
    end
    
    
    totalTime = toc;
    % Display computation time
    disp(['Current computation time: ', num2str(totalTime), ' seconds']); 
    
    save(['Results/Cross Session ReBaCCA for MDS ' dataNamePair{1} ' vs ' dataNamePair{2}], 'ReBaCCATotal', 'optimalKernelTotal',...
	    'selected_indices_together', 'num_selected', 'totalTime')

end