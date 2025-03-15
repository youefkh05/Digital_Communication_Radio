clc; 
clear;
close all;

% Parameters
A = 4;                  % Amplitude
N_realizations = 500;     % Number of waveforms (ensemble size)
num_bits = 100+1;         % Bits per waveform and one extra bit for shifting
bit_duration = 70e-3;   % 70 ms per bit
dac_interval = 10e-3;   % DAC updates every 10 ms
samples_per_bit = round(bit_duration / dac_interval);  % 7 samples per bit
total_time = num_bits * bit_duration;                  % Total waveform duration
t = 0:dac_interval:(total_time - dac_interval);        % Time vector

% Preallocate matrices for efficiency
Unipolar_All = zeros(N_realizations, length(t), 'int8');
PolarNRZ_All = zeros(N_realizations, length(t), 'int8');
PolarRZ_All = zeros(N_realizations, length(t), 'int8');

% Generate and store 500 realizations
for i = 1:N_realizations
    Data = randi([0, 1], 1, num_bits, 'int8');  % Random bit sequence
    
    %encode the data
    [Unipolar, PolarNRZ, PolarRZ] = generate_linecodes(Data, A, samples_per_bit);
    
    % Store in matrices
    Unipolar_All(i,:) = Unipolar;
    PolarNRZ_All(i,:) = PolarNRZ;
    PolarRZ_All(i,:) = PolarRZ;
    
    % Plot the first realization
    if i==1
        plot_linecodes(Data, Unipolar, PolarNRZ, PolarRZ, t, num_bits-1,'Realization 1');       
    end
end



% Plot the first 5 realizations as a sample
figure;
for i = 1:5
    subplot(5,1,i);
    stairs(t, Unipolar_All(i,:), 'b', 'LineWidth', 1.5);
    title(['Unipolar NRZ - Realization ' num2str(i)]);
    grid on;
end
xlabel('Time (s)');

% Apply random shift to all realizations
[Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted] = apply_random_shift_fixed_size(Unipolar_All, PolarNRZ_All, PolarRZ_All, samples_per_bit);

t_shifted = t(1:length(Unipolar_Shifted)); % Ensure the time vector matches

%plot after shift
plot_linecodes(Data, Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted, t_shifted, num_bits-1, 'Realization 1 shifted');

% Convert to double for accuracy
Unipolar_Shifted = double(Unipolar_Shifted);
PolarNRZ_Shifted = double(PolarNRZ_Shifted);
PolarRZ_Shifted = double(PolarRZ_Shifted);

%calculate the mean acrros time (question 1)
Unipolar_Mean = calculate_mean(Unipolar_Shifted);
PolarNRZ_Mean = calculate_mean(PolarNRZ_Shifted);
PolarRZ_Mean = calculate_mean(PolarRZ_Shifted);

%plot the mean across time
plot_mean_waveforms(t_shifted, Unipolar_Mean, PolarNRZ_Mean, PolarRZ_Mean, A);

% Compute variance for each code line
Unipolar_Var = calculate_variance(Unipolar_Shifted);
PolarNRZ_Var = calculate_variance(PolarNRZ_Shifted);
PolarRZ_Var = calculate_variance(PolarRZ_Shifted);

%plot the vatiance across time
plot_variance(t_shifted, Unipolar_Var, PolarNRZ_Var, PolarRZ_Var);

% Determine max_lag dynamically
max_lag = size(Unipolar_Shifted, 2) - 1;

%calculate the autocorrelation (question 3)
[Unipolar_AutoCorr, PolarNRZ_AutoCorr, PolarRZ_AutoCorr] = compute_stat_autocorr(Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted, max_lag);   

%plot the autocorrelation
plot_autocorrelation(Unipolar_AutoCorr, PolarNRZ_AutoCorr, PolarRZ_AutoCorr, max_lag)

% Compute time autocorrelation
[R_unipolar_nrz_t, R_polar_nrz_t, R_polar_rz_t, taw2] = compute_time_autocorr(Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted);

% Plot the time autocorrelation
plot_time_autocorrelation(R_unipolar_nrz_t, R_polar_nrz_t, R_polar_rz_t, taw2);

% Compute time mean for each line code
Unipolar_TimeMean = compute_time_mean(Unipolar_Shifted);
PolarNRZ_TimeMean = compute_time_mean(PolarNRZ_Shifted);
PolarRZ_TimeMean = compute_time_mean(PolarRZ_Shifted);

% Plot the time mean
plot_time_mean(t_shifted, Unipolar_TimeMean, PolarNRZ_TimeMean, PolarRZ_TimeMean);

% Compute time mean for all realizations
Unipolar_TimeMean_All = arrayfun(@estimate_time_mean, Unipolar_Shifted, 'UniformOutput', true);
PolarNRZ_TimeMean_All = arrayfun(@estimate_time_mean, PolarNRZ_Shifted, 'UniformOutput', true);
PolarRZ_TimeMean_All = arrayfun(@estimate_time_mean, PolarRZ_Shifted, 'UniformOutput', true);

% Compute the final estimated time mean (average over all realizations)
Final_Unipolar_TimeMean = mean(Unipolar_TimeMean_All);
Final_PolarNRZ_TimeMean = mean(PolarNRZ_TimeMean_All);
Final_PolarRZ_TimeMean = mean(PolarRZ_TimeMean_All);

% Display results
disp(['Estimated Unipolar NRZ Time Mean: ', num2str(Final_Unipolar_TimeMean)]);
disp(['Estimated Polar NRZ Time Mean: ', num2str(Final_PolarNRZ_TimeMean)]);
disp(['Estimated Polar RZ Time Mean: ', num2str(Final_PolarRZ_TimeMean)]);

fs = 100; % Sampling frequency 
[BW_unipolar, BW_polarNRZ, BW_polarRZ] = estimate_all_bandwidths(R_unipolar_nrz, R_polar_nrz, R_polar_rz, fs);


%-----------------------Functions----------------------------

function [Unipolar, PolarNRZ, PolarRZ] = generate_linecodes(Data, A, samples_per_bit)
    % Ensure input Data is of type int8
    Data = int8(Data);
    
    % Convert samples_per_bit to double for safe calculations
    samples_per_bitd = double(samples_per_bit);
    
    % Unipolar NRZ: 0 → 0V, 1 → A
    Unipolar = int8(Data * A);
    Unipolar = repelem(Unipolar, samples_per_bit); % Repeat each bit for duration
    
    % Polar NRZ: 0 → -A, 1 → +A
    PolarNRZ = int8((2 * Data - 1) * A);
    PolarNRZ = repelem(PolarNRZ, samples_per_bit);

    % Polar Return-to-Zero (RZ): Same as Polar NRZ but second half set to 0
    PolarRZ = PolarNRZ;
   
    % Apply RZ rule: second half of each bit period should be zero
    i = length(Data);  % Start from the last bit
    while i > 0
        end_idx = i * samples_per_bitd;  % Last sample of the bit
        start_idx = end_idx - floor(samples_per_bitd / 2) + 1;  % Start of the second half
        PolarRZ(start_idx:end_idx) = 0;  % Set the second half of the bit period to zero
        i = i - 1;  % Move to the previous bit
    end
end

function plot_linecodes(Data, Unipolar, PolarNRZ, PolarRZ, t, num_bits_to_show, plot_title)
    % Ensure num_bits_to_show does not exceed the actual number of bits
    num_samples_per_bit = ceil(length(t) / length(Data));
    num_samples_to_show = num_bits_to_show * num_samples_per_bit;

    % Trim the signals to display only the required number of bits
    t_show = t(1:num_samples_to_show);
    Unipolar_show = Unipolar(1, 1:num_samples_to_show); % Select row 1 explicitly
    PolarNRZ_show = PolarNRZ(1, 1:num_samples_to_show);
    PolarRZ_show = PolarRZ(1, 1:num_samples_to_show);

    % Convert Data into a sample-wise representation for accurate plotting
    Data_show = repelem(Data(1:num_bits_to_show), num_samples_per_bit);
    Data_t = t(1:length(Data_show)); % Adjust time axis

    % Plot the signals
    figure;
    sgtitle(plot_title); % Set a title for the entire figure
    
    subplot(4,1,1);
    stairs(Data_t, Data_show, 'k', 'LineWidth', 2);
    title('Orignal Data');
    ylim([-0.5, 1.5]); % Keep the binary level range
    yticks([0 1]);
    yticklabels({'0', '1'});
    grid on;

    subplot(4,1,2);
    stairs(t_show, Unipolar_show, 'b', 'LineWidth', 2);
    title('Unipolar NRZ');
    grid on;

    subplot(4,1,3);
    stairs(t_show, PolarNRZ_show, 'r', 'LineWidth', 2);
    title('Polar NRZ');
    grid on;

    subplot(4,1,4);
    stairs(t_show, PolarRZ_show, 'm', 'LineWidth', 2);
    title('Polar RZ');
    grid on;

    xlabel('Time (s)');
end

function [Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted] = apply_random_shift_fixed_size(Unipolar_All, PolarNRZ_All, PolarRZ_All, samples_per_bit)
    % Define parameters
    N_realizations = size(Unipolar_All, 1);  % 500 realizations
    extended_samples = size(Unipolar_All, 2); % 707 samples
    total_samples = 700;  % Fixed output size

    % Initialize shifted matrices
    Unipolar_Shifted = zeros(N_realizations, total_samples, 'int8');
    PolarNRZ_Shifted = zeros(N_realizations, total_samples, 'int8');
    PolarRZ_Shifted = zeros(N_realizations, total_samples, 'int8');

    % Apply random shift to each realization
    for i = 1:N_realizations
        % Generate random shift in range [0, samples_per_bit-1] samples
        random_shift_bits = randi([0, samples_per_bit-1]);

        % Extract shifted region
        Unipolar_Shifted(i, :) = Unipolar_All(i, random_shift_bits+1 : random_shift_bits+total_samples);
        PolarNRZ_Shifted(i, :) = PolarNRZ_All(i, random_shift_bits+1 : random_shift_bits+total_samples);
        PolarRZ_Shifted(i, :) = PolarRZ_All(i, random_shift_bits+1 : random_shift_bits+total_samples);
    end
end

function mean_waveform = calculate_mean(waveform_matrix)
    % Calculates the mean across all realizations without using the mean function
    % waveform_matrix: Matrix where each row is a realization
    
    [num_realizations, num_samples] = size(waveform_matrix); % Get matrix dimensions
    mean_waveform = sum(waveform_matrix, 1) / num_realizations; % Sum and divide by count
    
end

function plot_mean_waveforms(t, Unipolar_Mean, PolarNRZ_Mean, PolarRZ_Mean, A)
    % Function to plot the mean waveforms of different line codes across time
    % Inputs:
    %   t - Time vector
    %   Unipolar_Mean - Mean waveform of Unipolar NRZ
    %   PolarNRZ_Mean - Mean waveform of Polar NRZ
    %   PolarRZ_Mean - Mean waveform of Polar RZ
    %   A - Amplitude limit
    
    figure;

    subplot(3,1,1);
    plot(t, Unipolar_Mean, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Mean of Unipolar NRZ');
    grid on;
    ylim([-A, A]); % Set y-axis limits

    subplot(3,1,2);
    plot(t, PolarNRZ_Mean, 'g', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Mean of Polar NRZ');
    grid on;
    ylim([-A, A]);

    subplot(3,1,3);
    plot(t, PolarRZ_Mean, 'b', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Mean of Polar RZ');
    grid on;
    ylim([-A, A]);

    % Add a super title for the whole figure
    sgtitle('Mean of Different Line Codes');
end

function variance_waveform = calculate_variance(waveform_matrix)
     % Calculates the variance across all realizations (column-wise)
    % waveform_matrix: Matrix where each row is a realization (num_realizations x num_samples)
    % Returns: variance_waveform (1 x num_samples), representing the variance of each sample point
    
    % Compute the mean using the previously implemented function
    mean_waveform = calculate_mean(waveform_matrix);
    
    [num_realizations, num_samples] = size(waveform_matrix); % Get dimensions

    % Ensure mean_waveform is the same size for element-wise subtraction
    mean_waveform = repmat(mean_waveform, num_realizations, 1);
    
     % Compute variance manually using the variance formula
    variance_waveform = sum((waveform_matrix - mean_waveform).^2, 1) / num_realizations; % Population variance
    
end

function plot_variance(t, Unipolar_Var, PolarNRZ_Var, PolarRZ_Var)
    % Plots the variance of different line codes over time
    % t: Time vector
    % Unipolar_Var, PolarNRZ_Var, PolarRZ_Var: Variance waveforms

    figure;

    subplot(3,1,1);
    plot(t, Unipolar_Var, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Variance');
    title('Variance of Unipolar NRZ');
    grid on;

    subplot(3,1,2);
    plot(t, PolarNRZ_Var, 'g', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Variance');
    title('Variance of Polar NRZ');
    grid on;

    subplot(3,1,3);
    plot(t, PolarRZ_Var, 'b', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Variance');
    title('Variance of Polar RZ');
    grid on;

    % Add a super title for clarity
    sgtitle('Variance of Different Line Codes');
end

function [Unipolar_AutoCorr, PolarNRZ_AutoCorr, PolarRZ_AutoCorr] = compute_stat_autocorr(Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted, max_lag)
    % Compute Statistical Autocorrelation for given signals using calculate_mean
    % Inputs:
    %   Unipolar_Shifted - Shifted signal matrix for Unipolar NRZ
    %   PolarNRZ_Shifted - Shifted signal matrix for Polar NRZ
    %   PolarRZ_Shifted - Shifted signal matrix for Polar RZ
    % Outputs:
    %   Unipolar_AutoCorr - Computed autocorrelation for Unipolar NRZ
    %   PolarNRZ_AutoCorr - Computed autocorrelation for Polar NRZ
    %   PolarRZ_AutoCorr - Computed autocorrelation for Polar RZ


    % Set x-axis limits dynamically
    x_limit = max_lag / 10;

    % Initialize autocorrelation arrays
    Unipolar_AutoCorr = zeros(1, max_lag + 1);
    PolarNRZ_AutoCorr = zeros(1, max_lag + 1);
    PolarRZ_AutoCorr = zeros(1, max_lag + 1);

    % Compute mean autocorrelation using calculate_mean function
    for i = 0:max_lag
        Unipolar_AutoCorr(i+1) = calculate_mean(Unipolar_Shifted(:, 1) .* Unipolar_Shifted(:, i+1));
        PolarNRZ_AutoCorr(i+1) = calculate_mean(PolarNRZ_Shifted(:, 1) .* PolarNRZ_Shifted(:, i+1));
        PolarRZ_AutoCorr(i+1) = calculate_mean(PolarRZ_Shifted(:, 1) .* PolarRZ_Shifted(:, i+1));
    end

    % Time axis for plotting
    t = -max_lag:max_lag;

    % Compute symmetric autocorrelation values
    Unipolar_AutoCorr = [fliplr(Unipolar_AutoCorr), Unipolar_AutoCorr(2:end)];
    PolarNRZ_AutoCorr = [fliplr(PolarNRZ_AutoCorr), PolarNRZ_AutoCorr(2:end)];
    PolarRZ_AutoCorr = [fliplr(PolarRZ_AutoCorr), PolarRZ_AutoCorr(2:end)];

    
end

function plot_autocorrelation(Unipolar_AutoCorr, PolarNRZ_AutoCorr, PolarRZ_AutoCorr, max_lag)
    % Plots the statistical autocorrelation of Unipolar NRZ, Polar NRZ, and Polar RZ
    % Inputs:
    %   Unipolar_AutoCorr - Autocorrelation for Unipolar NRZ
    %   PolarNRZ_AutoCorr - Autocorrelation for Polar NRZ
    %   PolarRZ_AutoCorr - Autocorrelation for Polar RZ
    %   max_lag - Maximum lag value used for the time axis
    
    % Time axis
    t = -max_lag:max_lag;

    % Compute x-axis limit
    x_limit = max_lag / 10;

    % Plot results
    figure("name", "Statistical Autocorrelation");

    subplot(3,1,1);
    plot(t, Unipolar_AutoCorr, 'g');
    xlim([-x_limit x_limit]);
    xlabel("Time axis");
    ylabel("Autocorr axis");
    title("Unipolar NRZ Statistical Autocorrelation");

    subplot(3,1,2);
    plot(t, PolarNRZ_AutoCorr, 'b');
    xlim([-x_limit x_limit]);
    xlabel("Time axis");
    ylabel("Autocorr axis");
    title("Polar NRZ Statistical Autocorrelation");

    subplot(3,1,3);
    plot(t, PolarRZ_AutoCorr, 'r');
    xlim([-x_limit x_limit]);
    xlabel("Time axis");
    ylabel("Autocorr axis");
    title("Polar RZ Statistical Autocorrelation");
end
function [R_unipolar_nrz_t, R_polar_nrz_t, R_polar_rz_t, taw2] = compute_time_autocorr(UnipolarNRZ, PolarNRZ, PolarRZ)
    % Computes time autocorrelation for Unipolar NRZ, Polar NRZ, and Polar RZ waveforms
    % Inputs:
    %   UnipolarNRZ - Matrix containing all realizations for Unipolar NRZ
    %   PolarNRZ - Matrix containing all realizations for Polar NRZ
    %   PolarRZ - Matrix containing all realizations for Polar RZ
    % Outputs:
    %   R_unipolar_nrz_t - Time autocorrelation for Unipolar NRZ
    %   R_polar_nrz_t - Time autocorrelation for Polar NRZ
    %   R_polar_rz_t - Time autocorrelation for Polar RZ
    %   taw2 - Lag values

    % Get number of realizations and samples
    [num_realizations, num_samples] = size(UnipolarNRZ);

    % Define range of time lags
    max_lag = floor(num_samples / 2); % 150 in your case
    taw2 = -max_lag:max_lag;

    % Initialize autocorrelation vectors
    R_unipolar_nrz_t = zeros(1, length(taw2));     
    R_polar_nrz_t = zeros(1, length(taw2));     
    R_polar_rz_t = zeros(1, length(taw2));     

    % Loop through lag values
    for i = taw2
        M = i + max_lag + 1; % Shift index to fit within array bounds

        % Initialize sum for each lag
        sum_unipolar = 0;
        sum_polar_nrz = 0;
        sum_polar_rz = 0;

        % Compute sum over all realizations
        for n1 = max_lag + 1:num_samples
            if (n1 + i > 0) && (n1 + i <= num_samples) % Avoid out-of-bounds indexing
                sum_unipolar = sum_unipolar + mean(UnipolarNRZ(:, n1) .* UnipolarNRZ(:, n1 + i));
                sum_polar_nrz = sum_polar_nrz + mean(PolarNRZ(:, n1) .* PolarNRZ(:, n1 + i));
                sum_polar_rz = sum_polar_rz + mean(PolarRZ(:, n1) .* PolarRZ(:, n1 + i));
            end
        end

        % Normalize autocorrelation
        valid_samples = num_samples - abs(i);
        R_unipolar_nrz_t(M) = sum_unipolar / valid_samples;
        R_polar_nrz_t(M) = sum_polar_nrz / valid_samples;
        R_polar_rz_t(M) = sum_polar_rz / valid_samples;
    end
end
function plot_time_autocorrelation(R_unipolar, R_polarNRZ, R_polarRZ, taw2)
    figure('Name', 'Time Autocorrelation');

    subplot(3,1,1);
    plot(taw2, R_unipolar, 'b', 'LineWidth', 1.5);
    grid on;
    axis([-25 25 -4 20]);
    xlabel('Time (samples)');
    ylabel('Magnitude');
    title('Time Autocorrelation of Unipolar NRZ');

    subplot(3,1,2);
    plot(taw2, R_polarNRZ, 'r', 'LineWidth', 1.5);
    grid on;
    axis([-25 25 -4 20]);
    xlabel('Time (samples)');
    ylabel('Magnitude');
    title('Time Autocorrelation of Polar NRZ');

    subplot(3,1,3);
    plot(taw2, R_polarRZ, 'g', 'LineWidth', 1.5);
    grid on;
    axis([-25 25 -4 20]);
    xlabel('Time (samples)');
    ylabel('Magnitude');
    title('Time Autocorrelation of Polar RZ');
end

function TimeMean = compute_time_mean(waveform_matrix)
    % Computes the time mean for all realizations of a given waveform
    % Inputs:
    %   waveform_matrix - Matrix where each row represents a realization
    % Output:
    %   TimeMean - The computed time mean across all realizations

    % Get the number of realizations and samples
    [num_realizations, num_samples] = size(waveform_matrix);

    % Compute time mean (mean along columns)
    TimeMean = sum(waveform_matrix, 1) / num_realizations;
end
function plot_time_mean(t, Unipolar_Mean, PolarNRZ_Mean, PolarRZ_Mean)
    % Plots the time mean of different line codes over time
    % Inputs:
    %   t - Time vector
    %   Unipolar_Mean - Time mean of Unipolar NRZ
    %   PolarNRZ_Mean - Time mean of Polar NRZ
    %   PolarRZ_Mean - Time mean of Polar RZ

    figure('Name', 'Time Mean of Different Line Codes');

    subplot(3,1,1);
    plot(t, Unipolar_Mean, 'r', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Mean of Unipolar NRZ');
    grid on;

    subplot(3,1,2);
    plot(t, PolarNRZ_Mean, 'g', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Mean of Polar NRZ');
    grid on;

    subplot(3,1,3);
    plot(t, PolarRZ_Mean, 'b', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Time Mean of Polar RZ');
    grid on;
end

function TimeMean = estimate_time_mean(signal)
    % Estimates the time mean of a single waveform realization.
    % Input:
    %   signal - A single realization of the waveform (1D array)
    % Output:
    %   TimeMean - The estimated time mean of the waveform

    % Compute time mean as the average over all time samples
    TimeMean = sum(signal) / length(signal);
end
function [BW_unipolar, BW_polarNRZ, BW_polarRZ] = estimate_all_bandwidths(R_unipolar_nrz, R_polar_nrz, R_polar_rz, fs)
    % Function to estimate the bandwidth for Unipolar NRZ, Polar NRZ, and Polar RZ
    % Inputs:
    %   R_unipolar_nrz  - Autocorrelation of Unipolar NRZ
    %   R_polar_nrz  - Autocorrelation of Polar NRZ
    %   R_polar_rz  - Autocorrelation of Polar RZ
    %   fs  - Sampling frequency
  
    % Compute FFT of autocorrelation to get Power Spectral Density (PSD)
    fft_autocorr_unipolar = fft(R_unipolar_nrz);
    fft_autocorr_shifted_unipolar = fftshift(fft_autocorr_unipolar);
    
    fft_autocorr_polarNRZ = fft(R_polar_nrz);
    fft_autocorr_shifted_polarNRZ = fftshift(fft_autocorr_polarNRZ);
    
    fft_autocorr_polarRZ = fft(R_polar_rz);
    fft_autocorr_shifted_polarRZ = fftshift(fft_autocorr_polarRZ);

    % Frequency axis
    N1 = length(fft_autocorr_shifted_unipolar);
    frequencies1 = linspace(-fs/2, fs/2, N1);
    
    N2 = length(fft_autocorr_shifted_polarNRZ);
    frequencies2 = linspace(-fs/2, fs/2, N2);
    
    N3 = length(fft_autocorr_shifted_polarRZ);
    frequencies3 = linspace(-fs/2, fs/2, N3);

    % Compute PSD (absolute value and normalize by fs)
    PSD_unipolar = abs(fft_autocorr_shifted_unipolar) / fs;
    PSD_polarNRZ = abs(fft_autocorr_shifted_polarNRZ) / fs;
    PSD_polarRZ = abs(fft_autocorr_shifted_polarRZ) / fs;

    % Estimate Bandwidth using -3dB threshold
    BW_unipolar = find_bandwidth(PSD_unipolar, frequencies1);
    BW_polarNRZ = find_bandwidth(PSD_polarNRZ, frequencies2);
    BW_polarRZ = find_bandwidth(PSD_polarRZ, frequencies3);

    % Plot PSDs
    figure("Name", "PSD Plot");

    subplot(3,1,1);
    plot(frequencies1, PSD_unipolar, 'g');
    xlabel("Frequency (Hz)");
    ylabel("PSD");
    ylim([0 2]);
    title("Unipolar NRZ PSD");

    subplot(3,1,2);
    plot(frequencies2, PSD_polarNRZ, 'b');
    xlabel("Frequency (Hz)");
    ylabel("PSD");
    title("Polar NRZ PSD");

    subplot(3,1,3);
    plot(frequencies3, PSD_polarRZ, 'r');
    xlabel("Frequency (Hz)");
    ylabel("PSD");
    title("Polar RZ PSD");
end

function BW = find_bandwidth(PSD, freq_axis)
    % Function to find the -3dB bandwidth from PSD
    PSD_max = max(PSD);       % Find max power
    threshold = PSD_max / 2;  % -3dB threshold (half power)

    % Find indices where PSD is above threshold
    valid_indices = find(PSD >= threshold);

    if isempty(valid_indices)
        BW = 0; % If no valid range is found, return 0
        return;
    end

    % Compute bandwidth as the difference between first and last crossing points
    BW = abs(freq_axis(max(valid_indices)) - freq_axis(min(valid_indices)));
end







