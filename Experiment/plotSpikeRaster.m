% Define function to plot spike raster
function plotSpikeRaster(spike_trains, num_channels, duration_s, title_str, split_index)
    % spike_trains: channel x time, 0 1 matrix
    hold on;
    for channel = 1:num_channels
        % Find spike times for this channel (convert to seconds)
        spikes = find(spike_trains(channel, :)) / 1000;
        if isempty(spikes)
            continue; % Skip if no spikes in this channel
        end
        % Prepare x and y data for a single line object per channel
        x = reshape([spikes; spikes; nan(1, length(spikes))], [], 1);
        y = repmat([channel-0.4; channel+0.4; nan], length(spikes), 1);
        % Set color based on channel position relative to split_index
        if channel <= split_index
            color = 'k'; % Black for channels up to split_index
        else
            color = 'r'; % Red for channels after split_index
        end
        % Plot all spikes for this channel as one line object
        line(x, y, 'Color', color);
    end
    hold off;
    xlim([0 duration_s]);
    ylim([0.5, num_channels + 0.5]);
    xlabel('Time (s)');
    ylabel('Neuron index');
    title(title_str);
    % yticks(0:num_channels) % Uncomment if you want specific y-ticks
end