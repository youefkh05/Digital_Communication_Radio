clc; 
clear;
close all;

% Parameters
A = 4;                  % Amplitude
N_realizations = 1;     % Number of waveforms (ensemble size)
num_bits = 100;         % Bits per waveform
bit_duration = 70e-3;   % 70 ms per bit
dac_interval = 10e-3;   % DAC updates every 10 ms
samples_per_bit = int8(bit_duration / dac_interval);  % 7 samples per bit
total_time = num_bits * bit_duration;           % Total waveform duration
t = 0:dac_interval:(total_time - dac_interval); % Time vector

% Generate a random bit sequence
Data = randi([0, 1], 1, num_bits, 'int8');

[Unipolar, PolarNRZ, PolarRZ] = generate_linecodes(Data, A, samples_per_bit);

plot_linecodes(Data, Unipolar, PolarNRZ, PolarRZ, t, 6);

%-----------------------Functions----------------------------

function [Unipolar, PolarNRZ, PolarRZ] = generate_linecodes(Data, A, samples_per_bit)
    % Ensure input Data is of type int8
    Data = int8(Data);
    
    % Unipolar NRZ: 0 → 0V, 1 → A
    Unipolar = int8(Data * A);
    Unipolar = repelem(Unipolar, samples_per_bit); % Repeat each bit for duration
    
    % Polar NRZ: 0 → -A, 1 → +A
    PolarNRZ = int8((2 * Data - 1) * A);
    PolarNRZ = repelem(PolarNRZ, samples_per_bit);
    samples_per_bitd= double(samples_per_bit);
    % Polar Return-to-Zero (RZ): Same as Polar NRZ but second half set to 0
    PolarRZ = PolarNRZ;
   
    i = length(Data);  % Start from the last bit
    while i > 0
        end_idx = i * samples_per_bitd;  % Last sample of the bit
        start_idx = end_idx - double(floor(samples_per_bit / 2)) + 1;  % Start of the second half
        PolarRZ(start_idx:end_idx) = 0;  % Set the second half of the bit period to zero
        i = i - 1;  % Move to the previous bit
    end
end


function plot_linecodes(Data, Unipolar, PolarNRZ, PolarRZ, t, num_bits_to_show)
    % Ensure num_bits_to_show does not exceed the actual number of bits
    num_samples_per_bit = length(t) / length(Data);
    num_samples_to_show = num_bits_to_show * num_samples_per_bit;

    % Trim the signals to display only the required number of bits
    t_show = t(1:num_samples_to_show);
    Unipolar_show = Unipolar(1:num_samples_to_show);
    PolarNRZ_show = PolarNRZ(1:num_samples_to_show);
    PolarRZ_show = PolarRZ(1:num_samples_to_show);

    % Convert Data into a sample-wise representation for accurate plotting
    Data_show = repelem(Data(1:num_bits_to_show), num_samples_per_bit);
    Data_t = t(1:length(Data_show)); % Adjust time axis

    % Plot the signals
    figure;
    subplot(4,1,1);
    stairs(Data_t, Data_show, 'k', 'LineWidth', 2);
    title('Binary Data');
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
