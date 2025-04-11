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

% Reverse the pulse shaping function p[n] to create the matched filter
p_matched = fliplr(p);

% Generate hold filter
hold_filter = make_hold_filter(Ts, samples_per_bit, 1);

% Plot
plot_pulse_shape(p, Ts, 'Triangle Pulse Shape Normalized');
plot_pulse_shape(hold_filter, Ts, 'Hold Filter Normalized');
plot_pulse_shape(p_matched, Ts, 'Triangle Pulse Matched filter');

% Generate random bits
bits = randi([0 1], 1, bits_Num);

% Define the time vector from 0 to (Samples_Num-1)*Ts seconds
t = 0:Ts/samples_per_bit:Ts/samples_per_bit*(Samples_Num-1);
[~, bit_stream, ~] = generate_Impulse_linecodes(bits, A, samples_per_bit);

%plot the bit straem
plot_linecode(t, bit_stream, 'Polar NRZ bit stream');

%-----------------------Requiernment 1----------------------------

% Generate the PAM waveform
y_tx =  generate_pam_waveform(bit_stream, p);

% Plot the waveform vs time
plot_pam_waveform(Ts, y_tx, "The transmittor Output");

% Convolve the input signal with the matched filter
y_matched = conv(y_tx, p_matched); 

% Convolve the input signal with the hold filter
y_hold = conv(y_tx, hold_filter); 

%using correlator
y_corr = correlate_RX(y_tx, p, samples_per_bit);

%plot the matched and hold ouptut vs time
plot_matched_vs_hold(Ts, y_matched, y_hold, samples_per_bit, bits_Num)

%plot the matched and correlator ouptut vs time
plot_matched_vs_correlator(Ts, y_matched, y_tx, p, samples_per_bit, bits_Num);



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

function plot_RX_waveform(Ts, y_rx,shift_delay, samples_per_bit, bits_Num, plot_title, color)
    % PLOT_RX_WAVEFORM Plots the received waveform (after filtering) vs time vector t and adds a zero amplitude line.
    %
    % Inputs:
    %   Ts            - Symbol duration (in seconds)
    %   y_rx          - Receiver output waveform (filtered signal)
    %   samples_per_bit - Number of samples per bit
    %   bits_Num      - Total number of bits in the signal
    %   plot_title    - Title of the plot
    
    % Add a zero element at the start to ensure correct length
    y_rx_z = shift_right_zero(y_rx, shift_delay);  % Append a zero at the start of the result
    N = length(y_rx_z);                  % Number of samples
    
    % Create the time vector from 0 to Ts*(N-1) with time step Ts/samples_per_bit
    % Adjust length of t to match the length of y_rx_z
    t = 0:Ts/samples_per_bit:Ts*(N-1);
    
    % Check if the length of the time vector is greater than the signal length
    if length(t) > length(y_rx_z)
        t = t(1:length(y_rx_z));  % Trim t if necessary to match the length of y_rx_z
    end
    
    % Plot the waveform using the time vector t

    % Plot the matched filter output 
    h1 = plot(t, y_rx_z, 'Color', color, 'LineWidth', 2);   % Blue line for matched filter output
    title(plot_title);
    xlabel('Time [Ts sec]');  % Using time as x-axis
    ylabel('Amplitude');
    grid on;
    
    % Add horizontal line at zero amplitude
    yline(0, 'k', 'LineWidth', 1.5);  % Add a black horizontal line at zero
    
    % Add dots along the time axis (every Ts*samples_per_bit)
    hold on;
    % Initialize the dot_values vector with zeros
    dot_values = zeros(1, N);  % Initialize with zeros
    
    % Set the value of dot_values at each symbol's time step (every Ts*samples_per_bit)
    dot_values(1:samples_per_bit:end) = y_rx_z(1:samples_per_bit:end);  % Update only at the sampled points
    
    % Initialize dot_time for every time step
    dot_time = t;  % Time points corresponding to all samples
    

    % Plot the hollow dots (red outline, no fill)
    plot(dot_time, dot_values, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
    
    % Add vertical lines and dots manually at each sampled point
    sample_indices = 1:samples_per_bit:N;
    for i = 1:length(sample_indices)
        x = t(sample_indices(i));
        y = dot_values(sample_indices(i));

        % Draw vertical red line from 0 to sample value
        h2 = line([x x], [0 y], 'Color', 'r', 'LineWidth', 1.5);  % Save handle to one line (for legend)

        % Add a hollow red dot at the top of the line
        plot(x, y, 'o', 'MarkerSize', 6, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'w');
    end

    
    % Set the x-axis ticks and labels
    xticks(0:Ts:max(t));               % Set ticks at multiples of Ts
    xticklabels(0:length(xticks)-1);    % Label ticks as 0, 1, 2, 3, ...
    
    % Set the x-axis limits as [0, (N-samples_per_bit+1)/10]
    xlim([0 2*(N - samples_per_bit )/bits_Num]);
    
    % Add the legend with the custom line objects
    legend([h1, h2], 'Matched Filter Output', 'Sampled Output', 'Location', 'best', 'Box', 'on');
    
    hold off;

end

function y_correlated = correlate_RX(y, p, Ts)
% CORRELATE_RX Performs correlation using integrate-and-dump method
% 
% Inputs:
%   y   - Received signal (vector)
%   p   - Pulse shape (vector), length must equal Ts
%   Ts  - Number of samples per symbol (samples per bit)
%
% Output:
%   y_correlated - Correlated output per symbol

    % Ensure p is the same size as Ts
    if length(p) ~= Ts
        error('Length of pulse shape p must equal Ts (samples per symbol).');
    end

    % Pointwise multiply received signal by the pulse shape
    y_matched = y .* repmat(p, 1, ceil(length(y)/Ts));
    y_matched = y_matched(1:length(y));  % Trim any extra points from padding

    % Use built-in intdump for integrate-and-dump
    y_correlated = intdump(y_matched, Ts);  % returns one value per symbol
    
    % Upsample the correlated result to original length (for plotting, etc.)
    y_correlated = upsample(y_correlated, Ts);  % Insert Ts-1 zeros between samples
end

function y_shifted = shift_right_zero(y, shift_amt)
%SHIFT_RIGHT_ZERO Right-shifts the array by 'shift_amt' and fills with zeros at the beginning
%
% Inputs:
%   y          - Input vector
%   shift_amt  - Number of positions to shift
%
% Output:
%   y_shifted  - Shifted vector with zero-padding

    if shift_amt >= length(y)
        y_shifted = [zeros(1, shift_amt), y];  % Shift exceeds length: just prepend all zeros
    else
        y_shifted = [zeros(1, shift_amt), y(1:end)];
    end
end

function plot_matched_vs_correlator(Ts, y_filtered, y_tx, p, samples_per_bit, bits_Num)
    % PLOT_MATCHED_AND_CORRELATOR Create a figure with two subplots: 
    % one for the matched filter output and one for the correlator output.
    %
    % Inputs:
    %   Ts            - Symbol duration (in seconds)
    %   y_filtered    - Filtered receiver output signal
    %   y_tx          - Transmitted signal
    %   p             - Matched filter pulse
    %   samples_per_bit - Number of samples per bit
    %   bits_Num      - Total number of bits in the signal

    % Create a new figure for the subplots
    figure;
    
    % Create the first subplot (Matched Filter Output)
    subplot(2, 1, 1);  % Two rows, one column, first subplot
    plot_RX_waveform(Ts, y_filtered, 1, samples_per_bit, bits_Num, "The receiver output due to matched filter", 'g');

    % Create the second subplot (Correlator Output)
    subplot(2, 1, 2);  % Two rows, one column, second subplot
    y_corr = correlate_RX(y_tx, p, samples_per_bit);  % Get the correlated signal
    plot_RX_waveform(Ts, y_corr, 0, samples_per_bit, bits_Num, "The correlator Output", 'b');

end

function plot_matched_vs_hold(Ts, y_matched, y_hold, samples_per_bit, bits_Num)
    % PLOT_MATCHED_AND_CORRELATOR Create a figure with two subplots: 
    % one for the matched filter output and one for the correlator output.
    %
    % Inputs:
    %   Ts            - Symbol duration (in seconds)
    %   y_hold    - Filtered receiver output signal
    %   y_tx          - Transmitted signal
    %   p             - Matched filter pulse
    %   samples_per_bit - Number of samples per bit
    %   bits_Num      - Total number of bits in the signal

    % Create a new figure for the subplots
    figure;
    
    % Create the first subplot (Matched Filter Output)
    subplot(2, 1, 1);  % Two rows, one column, first subplot
    plot_RX_waveform(Ts, y_matched, 1, samples_per_bit, bits_Num, "The receiver output due to matched filter", 'c');

    % Create the second subplot (Correlator Output)
    subplot(2, 1, 2);  % Two rows, one column, second subplot
    plot_RX_waveform(Ts, y_hold, 0, samples_per_bit, bits_Num, "The receiver output due to hold filter", 'm');

end


function hold_filter = make_hold_filter(Ts, N, A)
    % MAKE_HOLD_FILTER Creates a flat, energy-normalized filter with scaling factor A.
    %
    % Inputs:
    %   Ts  - Duration of the filter (in seconds)
    %   N   - Number of samples
    %   A   - Amplitude scaling factor
    %
    % Outputs:
    %   hold_filter    - Filter impulse response (flat and energy normalized)

    % Time vector for the filter, from 0 to Ts, divided into N samples
    t = 0:Ts/N:Ts*(N-1)/N;  % From 0 to Ts with N samples
    
    % Filter impulse response (flat and energy normalized)
    hold_filter = A * ones(size(t)) / sqrt(Ts*N);  % Amplitude scaled 
end


