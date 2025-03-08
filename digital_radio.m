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

disp(Unipolar_Shifted(1, 1:10)); % Display the first 10 samples
disp(PolarNRZ_Shifted(1, 1:10));
disp(PolarRZ_Shifted(1, 1:10));

size(Unipolar_Shifted)
size(t)
size(t_shifted)


%plot after shift
plot_linecodes(Data, Unipolar_Shifted, PolarNRZ_Shifted, PolarRZ_Shifted, t_shifted, num_bits-1, 'Realization 1 shifted');


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


