clc;
clear;
close all;

% Parameters
samples_per_bit = 5;         % Number of samples
bits_Num = 10;               % 10 bits
Ts = 1;                      % Sampling interval
A = 5;                       % Peak amplitude at time = 0
Samples_Num=samples_per_bit* bits_Num;

% Generate triangle pulse
[p,denorm_p] = triangle_pulse(samples_per_bit);

% Plot
plot_pulse_shape(p, Ts, 'Triangle Pulse Shape Normalized');

% Generate random bits
bits = randi([0 1], 1, bits_Num);

% Define the time vector from 0 to (Samples_Num-1)*Ts seconds
t = 0:Ts/samples_per_bit:Ts/samples_per_bit*(Samples_Num-1);
[~, bit_stream, ~] = generate_Impulse_linecodes(bits, 1, samples_per_bit);

%plot the bit straem
plot_linecode(t, bit_stream, 'Polar NRZ bit stream');

% Generate the PAM waveform
y_tx =  denorm_p*generate_pam_waveform(bit_stream, p);

% Plot the waveform vs time
plot_pam_waveform(Ts, y_tx, "The transmittor Output");

% Reverse the pulse shaping function p[n] to create the matched filter
p_matched = fliplr(p);





%-----------------------Functions----------------------------

function [p,denorm_p] = triangle_pulse(N)
% TRIANGLE_PULSE Generates a smooth triangle pulse
%
% Inputs:
%   N  - Number of samples
%
% Output:
%   p  - Triangle pulse vector (normalized energy)

    % Linearly decreasing integers from N to 1
    p = zeros(1, N);
    for i = 1:N
        p(i) = (N - i + 1);
    end
    
    denorm_p=norm(p);
    % Normalize to unit energy
    p = p / norm(p);
end

function plot_pulse_shape(p, Ts, plot_title)
% PLOT_PULSE_SHAPE Plots a pulse shape over time
%
% Inputs:
%   p          - Pulse shape vector
%   Ts         - Total duration (in seconds) of the pulse
%   plot_title - Optional plot title (string)

    N = length(p);                     % Number of samples
    t = 0:Ts/N:Ts*(N-1)/N;                 % Discrete time vector

    figure;
    plot(t, p, 'LineWidth', 2);
    if nargin >= 3
        title(plot_title);
    else
        title('Pulse Shape');
    end
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

function [Unipolar, PolarNRZ, PolarRZ] = generate_Impulse_linecodes(Data, A, samples_per_bit)
% GENERATE_LINECODES Generates Unipolar, Polar NRZ, and Polar RZ waveforms.
%
% Inputs:
%   Data            - Array of bits (0s and 1s)
%   A               - Amplitude value
%   samples_per_bit - Number of samples per bit (e.g., 5 for 200ms sampling)
%
% Outputs:
%   Unipolar        - Unipolar NRZ waveform
%   PolarNRZ        - Polar NRZ waveform
%   PolarRZ         - Polar RZ waveform

    % Unipolar: 0 -> 0, 1 -> +A
    unipolar_symbols = A * Data;                    % Convert bits to unipolar symbols
    Unipolar = upsample(unipolar_symbols, samples_per_bit);  % Upsample for fine resolution
    
    % Polar NRZ: 0 -> -A, 1 -> +A
    polar_symbols = A * (2 * Data - 1);             % Converts 0→-A, 1→+A
    PolarNRZ = upsample(polar_symbols, samples_per_bit);  % Upsample for fine resolution

    % Polar RZ: same amplitude as Polar NRZ but with half-bit duration pulse
    % → Non-zero followed by zero for each bit
    PolarRZ = zeros(1, length(Data) * samples_per_bit);  % Initialize a zero array for Polar RZ
    for i = 1:length(Data)
        PolarRZ((i-1)*samples_per_bit + 1) = polar_symbols(i); % Only the first sample for each bit is non-zero
    end
end

function plot_linecode(t, linecode, plot_title)
% PLOT_LINECODE Plots a single linecode signal over time.
%
% Inputs:
%   t          - Time vector
%   linecode   - Signal array representing the linecode
%   plot_title - Title string for the plot

    % Plot the signal
    figure;
    stairs(t, linecode, 'LineWidth', 2);
    title(plot_title);
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
end

function y_tx = generate_pam_waveform(upsampled_bit_stream, p)
    % Converts upsampled bit stream to PAM waveform and convolves with pulse shape p
    % upsampled_bit_stream: upsampled binary bit stream array (e.g., [1 0 0 1 ...])
    % p: pulse shaping vector (e.g., [5 4 3 2 1]/sqrt(55))

    % Convolve with pulse shaping function
    y_tx = conv(upsampled_bit_stream, p);

    % Trim the result to match the length of the impulse train (original upsampled length)
    y_tx = y_tx(1:length(upsampled_bit_stream));  % Match the length of the upsampled input
    
 
end

function plot_pam_waveform(Ts, y_tx, plot_title)
    % PLOT_PAM_WAVEFORM Plots the PAM waveform using time vector t and adds a zero amplitude line.
    %
    % Inputs:
    %   t      - Time vector (in seconds)
    %   y_tx   - PAM waveform to be plotted
    
    % Add a zero element at the end
    y_tx_z = [y_tx, 0];  % Append a zero at the end of the result
    N = length(y_tx_z)-1;                     % Number of samples
    t = 0:10*Ts/N:Ts*N/5;                 % Discrete time vector
    
    % Plot the waveform using the time vector t
    figure;
    plot(t, y_tx_z, 'LineWidth', 2);
    title(plot_title);
    xlabel('Time (s)');  % Using time as x-axis
    ylabel('Amplitude');
    grid on;

    % Add horizontal line at zero amplitude
    yline(0, 'k', 'LineWidth', 1.5);  % 'k' is for a black line
end


