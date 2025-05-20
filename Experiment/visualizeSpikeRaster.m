clc,clear,close all

rng(42)

num_selected = 8;
dataNameTotal = {'1029_1', '1053_1', '1053_6'};
dataNameLabel = {'Rat A', 'Rat B session 1', 'Rat B session 2'};
trialTypeName={'_LEFT_nonmatch','_LEFT_sample','_RIGHT_nonmatch','_RIGHT_sample'};
trialTypeTitle = {'LN', 'LS', 'RN', 'RS'};
dataIdx = 0;
extractIdx = 500:5499;
for dataName = dataNameTotal
    dataIdx = dataIdx + 1;
    spike_trains=cell(length(dataName),length(trialTypeName));
    permuted_spike_trains=cell(length(dataName),length(trialTypeName));
    
    % Get all the spike train data
    for iTrialType=1:length(trialTypeName)
        concatenated_spike_trains=cell(size(dataName));
        concatenated_permuted_spike_trains=cell(size(dataName));
        for iData=1:length(dataName)
            load(['Preprocessed data/' dataName{iData} trialTypeName{iTrialType}])
    
            concatenated_spike_trains{iData}=success_spike_counts_concantenated;
            trial_length = timeAroundEvent * 1000;
            num_trials=size(success_spike_counts_concantenated,2)/trial_length;     
            
            for iTrial = 1:num_trials      
                spike_trains{iData, iTrialType}(iTrial, :, :) = success_spike_counts_concantenated(:,...
                    (iTrial - 1) * trial_length + 1 : iTrial * trial_length);
            end
        end
    end 

    % Step 1: Randomly sample trials from each target trial type
    selected_indices = cell(1, 4);       % Store indices for reference
    for iTrialType = 1:4
        num_trials = size(spike_trains{iTrialType}, 1);
        if num_trials < 8
            error(['Not enough trials (' num2str(num_trials) ') for trial type ' trialTypeName{iTrialType}]);
        end
        selected_indices{iTrialType} = randperm(num_trials, num_selected);
    end 
    
    for iTrialType = 1:4 
        tempSpikeTrains = spike_trains{iTrialType};
        tempIdx = selected_indices{iTrialType};
        for iTrial = 1:2   % Visualize 2 trials for each event
            subplot(3,8,(dataIdx-1)*8 + (iTrialType-1)*2 + iTrial)
            temp_spike_trains = squeeze(tempSpikeTrains(tempIdx(iTrial), :, :));
            temp_spike_trains = temp_spike_trains(:, extractIdx);
            plotSpikeRaster(temp_spike_trains, size(temp_spike_trains,1),...
                size(temp_spike_trains,2) / 1000, '', size(temp_spike_trains,1))
          
            if dataIdx == 1 
                title(trialTypeTitle{iTrialType})
            end
            if iTrialType == 1 && iTrial == 1
                ylabel(dataNameLabel{dataIdx}, 'Interpreter','none')
                yl = get(gca, 'YLabel');set(yl, 'Visible', 'on'); % Explicitly set visibility
            end
            axis off
            
        end
    end
    
end

export_pdf_figure(1,'Fig/Spike rasters for 3 datasets',1)