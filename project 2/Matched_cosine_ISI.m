clc;
clear;
close all;

% Parameters
N = 5+1;                   % Number of samples
bits_per_symbol = 10;      % 10 bits
Ts = 1;         % Sampling interval
A = 1;          % Peak amplitude at time = 0
Numofbits = bits_per_symbol * N;

% Generate triangle pulse
p = triangle_pulse(N, A);

% Define the time vector from 0 to Ts second
t = linspace(0, Ts, N);       % Time vector

% Plot
plot_pulse_shape(t, p);

% Generate random bits
bits = randi([0 1], 1, Numofbits);

% Generate the modulated waveform
y_tx = generate_pam_waveform(bits, p, N);

% Construct time vector
total_samples = length(y_tx);
Ts_symbol = 1;                         % Duration of one symbol (in seconds)
dt = Ts_symbol / N;                   % Time between samples
t_y_tx = (0:total_samples-1) * dt;    % Time vector for y_tx

% Plot the final waveform with time axis
figure;
plot(t_y_tx, y_tx, 'LineWidth', 1.5);
title('PAM Transmitted Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid 







%-----------------------Functions----------------------------

function p = triangle_pulse(N, A)
% TRIANGLE_PULSE Generates a smooth triangle pulse
%
% Inputs:
%   N  - Number of samples
%   A  - Amplitude of the triangle at time = 0
%
% Output:
%   p  - Triangle pulse vector (linearly decreasing from A to 0, not normalized)

    % Create a linearly decreasing vector from A to 0
    p = A*linspace(1, 0, N);
    
end

function plot_pulse_shape(t, p)
% PLOT_PULSE_SHAPE Plots a given pulse shape against time
%
% Inputs:
%   t - Time vector
%   p - Pulse shape vector

    plot(t, p, 'LineWidth', 2);
    title('Smooth Triangle Pulse Shape');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;

end

function y_tx = generate_pam_waveform(bits, p, N)
    % Converts bits to PAM impulses and convolves with pulse shape p
    % bits: binary array (0, 1)
    % p: pulse shaping vector (e.g., [5 4 3 2 1]/sqrt(55))
    % N: 5 in this case
    
    % Convert bits to symbols: 0 -> -1, 1 -> +1
    symbols = 2 * bits - 1;

    % Upsample: insert (N - 1) zeros between each symbol
    impulse_train = upsample(symbols, N);

    % Convolve with pulse shaping function
    y_tx = conv(impulse_train, p);
end

