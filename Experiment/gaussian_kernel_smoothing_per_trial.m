% Function to perform Gaussian kernel smoothing
function [smoothed_data,effectiveTrialLength] = gaussian_kernel_smoothing_per_trial(data, kernel_width)

    kernel_size = 2 * kernel_width;
    kernel = normpdf(-kernel_size:kernel_size, 0, kernel_width);
    kernel = kernel / sum(kernel); % Normalize the kernel

    % Apply convolution to each neuron's spike train
    for iTrial=1:size(data,3)
	    smoothed_data(:,:,iTrial)=conv2(data(:,:,iTrial),kernel,'same');
    end
	
    effectiveTrialLength=size(smoothed_data,2);
    smoothed_data=reshape(smoothed_data,size(smoothed_data,1),size(smoothed_data,2)*size(smoothed_data,3));
end