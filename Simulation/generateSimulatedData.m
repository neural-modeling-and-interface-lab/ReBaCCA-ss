clc,clear,close all

totalRepeatTime=8;

for iRepeat=1:totalRepeatTime
    
    % Parameters
    num_neurons1 = 20;
    num_neurons2 = 20;
    trial_length = 5000; % in ms
    half_trial_length = trial_length / 2;
    num_trials = 20; % Number of trials
    
    % Gaussian parameters for first dataset
    mu1_1 = [2000 2000]; % Mean for the first group of neurons
    sigma1_1 = 50; % Standard deviation for the first group of neurons
    mu2_1 = [3000, 3000]; % Mean for the second group of neurons
    sigma2_1 = 150; % Standard deviation for the second group of neurons 150
    split_index1 = 10;
    basicFiring1 = [10 10];
    backgroundFiring1 = [0.002 0.002];
    
    % Gaussian parameters for second dataset
    mu1_2 = [2000, 2000]; % Mean for the first group of neurons
    sigma1_2 = 50; % Standard deviation for the first group of neurons
    mu2_2 = [3000, 3000]; % Mean for the second group of neurons
    sigma2_2 = 150; % Standard deviation for the second group of neurons 150
    split_index2 = 16;
    basicFiring2 = [10 10];
    backgroundFiring2 = [0.002 0.002];
    
    % Generate datasets
    [spike_trains1, firing_rates1] = generate_spike_trains(basicFiring1, backgroundFiring1, num_neurons1, trial_length, mu1_1, sigma1_1, mu2_1, sigma2_1, split_index1, num_trials);
    [spike_trains2, firing_rates2] = generate_spike_trains(basicFiring2, backgroundFiring2, num_neurons2, trial_length, mu1_2, sigma1_2, mu2_2, sigma2_2, split_index2, num_trials);
    
    % Initialize permuted spike trains
    permuted_spike_trains1 = spike_trains1;
    permuted_spike_trains2 = spike_trains2;
    randIdx1 = randperm(trial_length);
    randIdx2 = randperm(trial_length);
    
    % Randomly permute spikes within each trial for permuted_spike_trains1
    for trial = 1:num_trials
        permuted_spike_trains1(:, :, trial) = spike_trains1(:, randperm(trial_length), trial);
    end
    
    % Randomly permute spikes within each trial for permuted_spike_trains2
    for trial = 1:num_trials
        permuted_spike_trains2(:, :, trial) = spike_trains2(:, randperm(trial_length), trial);
    end
    
    % Concatenate spike trains across all trials for the first dataset
    concatenated_spike_trains1 = reshape(permute(spike_trains1, [1 2 3]), num_neurons1, trial_length * num_trials);
    concatenated_permuted_spike_trains1 = reshape(permuted_spike_trains1, num_neurons1, trial_length * num_trials);
    
    % Concatenate spike trains across all trials for the second dataset
    concatenated_spike_trains2 = reshape(permute(spike_trains2, [1 2 3]), num_neurons2, trial_length * num_trials);
    concatenated_permuted_spike_trains2 = reshape(permuted_spike_trains2, num_neurons2, trial_length * num_trials);

    
    save(['Data/Simulated data_' num2str(iRepeat)]);
end

%% Plot the simulated data
figure
subplot(2,2,1);
plot((1:trial_length)/1000, firing_rates1(1:split_index1, :, 1),'b','Linewidth',2);
hold on;
plot((1:trial_length)/1000, firing_rates1(split_index1+1:end, :, 1),'r','Linewidth',2);
xlabel('');xticklabels('');
title('Firing Probability of Neurons - Dataset 1');
xlim([0 trial_length/1000])

subplot(2,2,2);
plot((1:trial_length)/1000, firing_rates2(1:split_index2, :, 1),'b','Linewidth',2);
hold on;
plot((1:trial_length)/1000, firing_rates2(split_index2+1:end, :, 1),'r','Linewidth',2);
xlabel('');xticklabels('');
title('Firing Probability of Neurons - Dataset 2');
xlim([0 trial_length/1000])

subplot(2,2,3);
plotSpikeRaster(spike_trains1(:, :, 2),num_neurons1,trial_length/1000,...
    'Simulated Spike Trains',split_index1)
xlabel('Time (s)')
xlim([0 trial_length/1000])

subplot(2, 2, 4);
plotSpikeRaster(spike_trains2(:, :, 2),num_neurons2,trial_length/1000,...
    'Simulated Spike Trains',split_index2)
xlabel('Time (s)')
xlim([0 trial_length/1000])

% Function to generate spike trains
function [spike_trains, firing_rates] = generate_spike_trains(basicFiring, backgroundFiring, num_neurons, trial_length, mu1, sigma1, mu2, sigma2, split_index, num_trials)
    time_vector = 1:trial_length;
    
    % Initialize firing rates
    firing_rates = zeros(num_neurons, trial_length, num_trials);

    % Generate firing rates
    for i = 1:num_neurons
        for trial = 1:num_trials
            if i <= split_index
                % Neurons up to split_index have two Gaussian firing rates in each trial
                firing_rate_1 = backgroundFiring(1) + normpdf(time_vector, mu1(1), sigma1) * basicFiring(1);
                firing_rate_2 = backgroundFiring(1) + normpdf(time_vector, mu1(2), sigma1) * basicFiring(1);
                firing_rates(i, :, trial) = (firing_rate_1 + firing_rate_2) / 2; % Average the two firing rates
            else
                % Neurons after split_index have two Gaussian firing rates in each trial
                firing_rate_1 = backgroundFiring(2) + normpdf(time_vector, mu2(1), sigma2) * basicFiring(2);
                firing_rate_2 = backgroundFiring(2) + normpdf(time_vector, mu2(2), sigma2) * basicFiring(2);
                firing_rates(i, :, trial) = (firing_rate_1 + firing_rate_2) / 2; % Average the two firing rates
            end
        end
    end

    % Simulate spike trains
    spike_trains = zeros(num_neurons, trial_length, num_trials);
    for trial = 1:num_trials
        for i = 1:num_neurons
            for t = 1:trial_length
                if rand < firing_rates(i, t, trial)  % Convert rate to probability per ms
                    spike_trains(i, t, trial) = 1;
                end
            end
        end
    end
end