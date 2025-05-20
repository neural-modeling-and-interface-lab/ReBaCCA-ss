% Function to perform Gaussian kernel smoothing
function smoothed_data = gaussian_kernel_smoothing(data, kernel_width)

    kernel_size = 2 * kernel_width;
    kernel = normpdf(-kernel_size:kernel_size, 0, kernel_width);
    kernel = kernel / sum(kernel); % Normalize the kernel

    % Apply convolution to each neuron's spike train
	smoothed_data=conv2(data,kernel,'same');
	
end