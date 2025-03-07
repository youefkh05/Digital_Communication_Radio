clc; 
clear;
close all;

% Parameters
A = 4;                  % Amplitude
N_realizations = 1;     % Number of waveforms (ensemble size)
num_bits = 100;         % Bits per waveform
bit_duration = 70e-3;   % 70 ms per bit
dac_interval = 10e-3;   % DAC updates every 10 ms
samples_per_bit = bit_duration / dac_interval;  % 7 samples per bit
total_time = num_bits * bit_duration;           % Total waveform duration
t = 0:dac_interval:(total_time - dac_interval); % Time vector

% Initialize arrays
Unipolar_ensemble = zeros(N_realizations, length(t));
PolarNRZ_ensemble = zeros(N_realizations, length(t));
RZ_ensemble = zeros(N_realizations, length(t));


