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

% theoretical calculation
Eb = denorm_p* denorm_p;     % Energy per bit (form the pulse shape)

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

% Define the time vector for the transmitted signal
t = 0:Ts/samples_per_bit:(Samples_Num-1)*Ts/samples_per_bit; % Adjusted time vector

[~, bit_stream, ~,~,bit_stream_symbols] = generate_Impulse_linecodes(bits, A, samples_per_bit);

%plot the bit straem
plot_linecode(t, bit_stream, 'Polar NRZ bit stream');

%-----------------------Requiernment 1----------------------------

% Generate the PAM waveform
y_tx =  generate_pam_waveform(bit_stream, p);

% Plot the waveform vs time
plot_pam_waveform(Ts, y_tx, 'The Transmitter Output' ,samples_per_bit);

% Convolve the input signal with the matched filter
y_matched = conv(y_tx, p_matched); 

% Convolve the input signal with the hold filter
y_hold = conv(y_tx, hold_filter); 

% using correlator
y_corr = correlate_RX(y_tx, p, samples_per_bit);

% plot the matched and hold ouptut vs time
[y_filtered_sampled ,y_corr_sampled] = plot_matched_vs_hold(Ts, y_matched, y_hold, samples_per_bit, bits_Num);

% plot the matched and correlator ouptut vs time
[~ ,y_hold_sampled] = plot_matched_vs_correlator(Ts, y_matched, y_tx, p, samples_per_bit, bits_Num);

% plot the Tx input vs All Rx output
plot_all_outputs_vs_input(Ts, bit_stream_symbols, y_filtered_sampled, y_corr_sampled, y_hold_sampled);

%-----------------------Requiernment 2----------------------------

bits_Num = 100;                % 100 bits
Samples_Num=samples_per_bit* bits_Num;

% Generate random bits
bits = randi([0 1], 1, bits_Num);

% Define the time vector for the transmitted signal
t = 0:Ts/samples_per_bit:(Samples_Num-1)*Ts/samples_per_bit; % Adjusted time vector

[~, bit_stream, ~,~,bit_stream_symbols] = generate_Impulse_linecodes(bits, A, samples_per_bit);

%plot the bit straem
plot_linecode(t, bit_stream, 'Polar NRZ bit stream');

% Generate the PAM waveform
y_tx =  generate_pam_waveform(bit_stream, p);

% Plot the waveform vs time
plot_pam_waveform(Ts, y_tx, 'The Transmitter Output' ,samples_per_bit);

% AWGN Channel
SNR_db  = -2;                % SNR given in dB from -2 dB to 5 dB in 1 dB steps
[y_tx_noise ,SNR_lin ,N0 , AWGN_scaled] = AddAWGN(y_tx, Eb, SNR_db);

% Convolve the input signal with the matched filter
y_matched = conv(y_tx_noise, p_matched); 

% Convolve the input signal with the hold filter
y_hold = conv(y_tx_noise, hold_filter); 

% using correlator
y_corr = correlate_RX(y_tx_noise, p, samples_per_bit);

% plot the matched and hold ouptut vs time
[y_filtered_sampled ,y_hold_sampled] = plot_matched_vs_hold(Ts, y_matched, y_hold, samples_per_bit, bits_Num);

% plot the matched and correlator ouptut vs time
[~ ,y_corr_sampled] = plot_matched_vs_correlator(Ts, y_matched, y_tx_noise, p, samples_per_bit, bits_Num);

% plot the Tx input vs All Rx output
plot_all_outputs_vs_input(Ts, bit_stream_symbols, y_filtered_sampled, y_corr_sampled, y_hold_sampled);

% Calculate the probability of error
polar_threshold = 0;
[BER, error_array] = calculate_error_probability(bit_stream_symbols, y_filtered_sampled, polar_threshold, A);

% Display the result
disp(['Probability of Error (BER) = ', num2str(BER)]);

% Calculate the  Theoretical probability of error
theoretical_BER = 0.5 * erfc(sqrt(SNR_lin));  % Theoretical BER

% Display the result
disp(['Theoretical BER: ', num2str(theoretical_BER)]);

SNR_db  = -2:1:5;                % SNR given in dB from -2 dB to 5 dB in 1 dB steps

% BER from matched vs hold
[BER_matched, BER_hold, theoretical_BER] = ...
    BER_vs_SNR(y_tx, Eb, p_matched, hold_filter, ...
    samples_per_bit, bits_Num, bit_stream_symbols, A, SNR_db);

% plot them
plot_BER_vs_EbN0(SNR_db, BER_matched, BER_hold, theoretical_BER);
% ===========================================================
% Project: ISI and Raised Cosine Filters - Eye Diagrams at A & B
% ===========================================================

% Parameters 
num_bits = 100;                
samples_per_symbol = 8;       
A = 1;

% Generate random bits and map: 0 -> -A, 1 -> +A
bits = randi([0 1], 1, num_bits);
symbols = A * (2 * bits - 1);  

% Upsample
tx_upsampled = upsample(symbols, samples_per_symbol);

% Configurations: [R, delay]
rolloff_values = [0, 0, 1, 1];
delay_values   = [2, 8, 2, 8];

% Figure counter
fig_num = 23;

for i = 1:4
    R = rolloff_values(i);
    delay = delay_values(i);

    % SRRC filter
    srrc_filter = rcosdesign(R, 2*delay, samples_per_symbol, 'sqrt');

    % Plot SRRC filter impulse response
    figure(fig_num);
    filter_length = length(srrc_filter); % Get the length of the filter
    t = (-(filter_length-1)/(2*samples_per_symbol)):(1/samples_per_symbol):((filter_length-1)/(2*samples_per_symbol)); % Adjusted time vector
    plot(t, srrc_filter, 'r-', 'LineWidth', 1.5); % Red line for consistency
    set(gca, 'Color', 'white');  
    set(gcf, 'Color', 'white');  
    title(['SRRC Filter - R = ' num2str(R) ', Delay = ' num2str(delay)]);
    set(get(gca, 'Title'), 'Color', 'black');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    set(get(gca, 'XLabel'), 'Color', 'black');
    set(get(gca, 'YLabel'), 'Color', 'black');
    xlabel('Time (symbols)');
    ylabel('Amplitude');
    grid on;
    fig_num = fig_num + 1;

    % Filtered signal at A
    tx_filtered = filter(srrc_filter, 1, tx_upsampled);

    % Matched filter (Rx)
    rx_filtered = filter(srrc_filter, 1, tx_filtered);

    % Remove filter transients
    filter_delay = 2 * delay * samples_per_symbol;
    valid_tx = tx_filtered(filter_delay+1:end);
    valid_rx = rx_filtered(filter_delay+1:end);

    % Plot eye diagram at A
    figure(fig_num);
    eyediagram(valid_tx, 2 * samples_per_symbol);
    set(gca, 'Color', 'white');  
    set(gcf, 'Color', 'white');  
    lines = findall(gca, 'Type', 'line');
    set(lines, 'Color', 'red'); 
    title(['Eye Diagram at Point A - R = ' num2str(R) ', Delay = ' num2str(delay)]);
    set(get(gca, 'Title'), 'Color', 'black');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    set(get(gca, 'XLabel'), 'Color', 'black');
    set(get(gca, 'YLabel'), 'Color', 'black');
    fig_num = fig_num + 1;

    % Plot eye diagram at B
    figure(fig_num);
    eyediagram(valid_rx, 2 * samples_per_symbol);
    set(gca, 'Color', 'white');
    set(gcf, 'Color', 'white');
    lines = findall(gca, 'Type', 'line');
    set(lines, 'Color', 'red');
    title(['Eye Diagram at Point B - R = ' num2str(R) ', Delay = ' num2str(delay)]);
    set(get(gca, 'Title'), 'Color', 'black');
    set(gca, 'XColor', 'black', 'YColor', 'black');
    set(get(gca, 'XLabel'), 'Color', 'black');
    set(get(gca, 'YLabel'), 'Color', 'black');
    fig_num = fig_num + 1;
end
%% 


%% 
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

function [Unipolar, PolarNRZ, PolarRZ, unipolar_symbols, polar_symbols] =...
    generate_Impulse_linecodes(Data, A, samples_per_bit)
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
    stem(t, linecode, 'LineWidth', 2);
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

function plot_pam_waveform(Ts, y_tx, plot_title,smaple_per_bit)
    % PLOT_PAM_WAVEFORM Plots the PAM waveform using time vector t and adds a zero amplitude line.
    %
    % Inputs:
    %   t      - Time vector (in seconds)
    %   y_tx   - PAM waveform to be plotted
    
    % Add a zero element at the end
    y_tx_z = [y_tx, 0];  % Append a zero at the end of the result
    N = length(y_tx_z);                     % Number of samples
    t = 0:Ts/smaple_per_bit:Ts*(N-1)/smaple_per_bit;                 % Discrete time vector
    
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

function [sampled_values] = ...
    plot_RX_waveform(Ts, y_rx,shift_delay, ...
    samples_per_bit, bits_Num, plot_title, color)
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
    xlim([0 ceil((N)/samples_per_bit)-1]);
    
    % Add the legend with the custom line objects
    legend([h1, h2], 'Matched Filter Output', 'Sampled Output', 'Location', 'best', 'Box', 'on');
    
    hold off;

    % Output: Return only the non-zero sampled values, skip first (padding), ensure row vector
    sampled_values = y_rx_z(sample_indices(2:end));
    sampled_values = sampled_values(:).';  % Force row vector
end

function y_correlated = correlate_RX(y, p, Ts)
% CORRELATE_RX Performs correlation using running integration over each symbol period.
%
% Inputs:
%   y   - Received signal (vector)
%   p   - Pulse shape (vector), must have length Ts
%   Ts  - Samples per symbol
%
% Output:
%   y_correlated - Correlated output (same length as input, running integration within each symbol)

    % Validate pulse shape length
    if length(p) ~= Ts
        error('Length of pulse shape p must equal Ts (samples per symbol).');
    end

    % Repeat the pulse shape to match the input signal length
    p_repeated = repmat(p(:), ceil(length(y)/Ts), 1);
    p_repeated = p_repeated(1:length(y));

    % Weighted input
    y_weighted = y(:) .* p_repeated;

    % Initialize the output
    y_correlated = zeros(size(y_weighted));

    % Process full symbol periods
    full_symbols = floor(length(y) / Ts);
    for i = 1:full_symbols
        idx_start = (i-1)*Ts + 1;
        idx_end = i*Ts;
        y_correlated(idx_start:idx_end) = cumsum(y_weighted(idx_start:idx_end));
    end

    % Handle remaining samples if any
    if mod(length(y), Ts) ~= 0
        idx_start = full_symbols*Ts + 1;
        y_correlated(idx_start:end) = cumsum(y_weighted(idx_start:end));
    end
end

function y_shifted = shift_right_zero(y, shift_amt)
%SHIFT_RIGHT_ZERO Right-shifts the array by 'shift_amt' and fills with zeros at the beginning
%
% Inputs:
%   y          - Input vector (row or column)
%   shift_amt  - Number of positions to shift
%
% Output:
%   y_shifted  - Shifted vector with zero-padding (same orientation as y)

    if isrow(y)
        y_shifted = [zeros(1, shift_amt), y];  % Concatenate as row
    else
        y_shifted = [zeros(shift_amt, 1); y];  % Concatenate as column
    end
end

function [y_filtered_sampled ,y_corr_smapled] = ...
    plot_matched_vs_correlator(Ts, y_filtered, y_tx, p,...
    samples_per_bit, bits_Num)
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
    
    y_filtered_sampled = plot_RX_waveform(Ts, y_filtered, 1, samples_per_bit, bits_Num, "The receiver output due to matched filter", 'g');

    % Create the second subplot (Correlator Output)
    hold on;
    y_corr = correlate_RX(y_tx, p, samples_per_bit);  % Get the correlated signal
    y_corr_smapled = plot_RX_waveform(Ts, y_corr, 1, samples_per_bit, bits_Num, "The correlator Output", 'b');
end

function [y_matched_smapled ,y_hold_sampled] = ...
    plot_matched_vs_hold(Ts, y_matched, y_hold, ...
    samples_per_bit, bits_Num)
    % PLOT_MATCHED_AND_HOLD Create a figure with two subplots: 
    % one for the matched filter output and one for the hold filter output.
    %
    % Inputs:
    %   Ts            - Symbol duration (in seconds)
    %   y_matched     - Matched filter output signal
    %   y_hold        - Hold filter output signal
    %   samples_per_bit - Number of samples per bit
    %   bits_Num      - Total number of bits in the signal

    % Create a new figure for the subplots
    figure;
    
    % Create the first subplot (Matched Filter Output)
    subplot(2, 1, 1);  % Two rows, one column, first subplot
    y_matched_smapled = plot_RX_waveform(Ts, y_matched, 1, samples_per_bit, bits_Num, "The receiver output due to matched filter", 'c');
    pause(3);
    
    % Create the second subplot (Hold Filter Output)
    subplot(2, 1, 2);  % Two rows, one column, second subplot
    y_hold_sampled = plot_RX_waveform(Ts, y_hold, 1, samples_per_bit, bits_Num, "The receiver output due to hold filter", 'm');
    pause(3);
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

function plot_all_outputs_vs_input(Ts, message, y_filtered_sampled, y_corr_sampled, y_hold_sampled)
    % PLOT_ALL_OUTPUTS_VS_INPUT
    %   Plots the original bit stream and the receiver sampled outputs in 1x4 subplots.
    %
    % Inputs:
    %   Ts                  - Symbol duration (in seconds)
    %   message             - Original input bit stream (1×N)
    %   y_filtered_sampled  - Output of matched filter sampled (1×N)
    %   y_corr_sampled      - Output of correlator sampled (1×N)
    %   y_hold_sampled      - Output of hold circuit sampled (1×N)

    % Validate input lengths
    N = length(message);
    if ~isequal(length(y_filtered_sampled), N) || ...
       ~isequal(length(y_corr_sampled), N) || ...
       ~isequal(length(y_hold_sampled), N)
        error('All input vectors must be the same length as the message.');
    end

    % Time vector
    t = 0:Ts:(N-1)*Ts;

    figure;

    % --- Subplot 1: Original message ---
    subplot(4, 1, 1);
    stem(t, message, 'k', 'filled', 'LineWidth', 1.5);
    title('Input Message');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    xlim([0 t(end)]);

    % --- Subplot 2: Matched Filter Output (Impulse representation) ---
    subplot(4, 1, 2);
    stem(t, y_filtered_sampled, 'b', 'filled', 'LineWidth', 2);
    title('Matched Filter Output');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    xlim([0 t(end)]);

    % --- Subplot 3: Correlator Output (Impulse representation) ---
    subplot(4, 1, 3);
    stem(t, y_corr_sampled, 'r', 'filled', 'LineWidth', 2);
    title('Correlator Output');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    xlim([0 t(end)]);

    % --- Subplot 4: Hold Output (Impulse representation) ---
    subplot(4, 1, 4);
    stem(t, y_hold_sampled, 'g', 'filled', 'LineWidth', 2);
    title('Hold Output');
    xlabel('Time (s)');
    ylabel('Amplitude');
    grid on;
    xlim([0 t(end)]);
end

function [Pe, error_array] = calculate_error_probability(input, output, threshold, A)
    % CALCULATE_ERROR_PROBABILITY
    %   Compares the detected symbols with the transmitted symbols using a threshold.
    %   Calculates the probability of error (BER) and provides an error array for visualization.
    %
    % Inputs:
    %   input      - Transmitted symbol stream (1xN) [-A, A]
    %   output     - Received output (1xN) after filtering and sampling
    %   threshold  - Decision threshold (e.g., for polar, it would be 0)
    %   A          - The amplitude used for the symbols (e.g., +A and -A)
    %
    % Outputs:
    %   Pe         - Probability of error (BER)
    %   error_array - Array of errors (1xN), 0 for no error, 1 for error

    % Decision rule (Apply threshold to output to decide the received symbols)
    decision = output > threshold;  % Decision rule: 1 if output > threshold, else 0
    decision = A*(2 * decision - 1);    % Map 1 to +A and 0 to -A (polar decision)

    % Compare the detected symbols (decision) with the transmitted symbols (input)
    error_array = (decision ~= input);  % 0 for matched, 1 for error
    
    % Calculate probability of error (BER)
    Pe = sum(error_array) / length(input);  % Ratio of errors to total bits
end

function [y_tx_noise ,SNR_lin ,N0 , AWGN_scaled] = AddAWGN(y_tx, Eb, SNR_db)
    % AddAWGN - Add Additive White Gaussian Noise (AWGN) to a transmitted signal
    % 
    % Inputs:
    %   y_tx       - Transmitted signal (1xN array)
    %   Eb         - Energy per bit (scalar)
    %   SNR_db     - Signal-to-Noise Ratio in dB (scalar)
    %
    % Output:
    %   y_tx_noise - Noisy received signal (1xN array)
    
    % Length of transmitted signal
    N_tx = length(y_tx);

    % Generate unity variance, zero mean white Gaussian noise
    AWGN = randn(1, N_tx);  % Generate noise with the same length as y_tx

    % Calculate the noise power spectral density based on Eb and SNR
    SNR_lin = 10^(SNR_db/10.0);  % Convert SNR from dB to linear scale
    N0 = Eb / SNR_lin;           % Calculate noise power (N0)
    
    % Scale the noise to have variance N0/2
    AWGN_scaled = sqrt(N0 / 2) * AWGN;

    % Add the noise to the transmitted signal
    y_tx_noise = y_tx + AWGN_scaled;  % Noisy received signal
end

function [BER_matched, BER_hold, theoretical_BER] = ...
    BER_vs_SNR(y_tx, Eb, p_matched, hold_filter, ...
    samples_per_bit, bits_Num, bit_stream_symbols, A, Eb_N0_dB)
    % BER_VS_SNR - Calculate and plot the Bit Error Rate (BER) vs Eb/N0 for different filters
    %
    % Inputs:
    %   y_tx             - Transmitted signal (1xN array)
    %   Eb               - Energy per bit (scalar)
    %   p_matched        - Matched filter (1xN array)
    %   hold_filter      - Hold filter (1xN array)
    %   samples_per_bit  - Number of samples per bit (scalar)
    %   bits_Num         - Number of bits (scalar)
    %   bit_stream_symbols - Transmitted bit stream symbols (1xN array)
    %   A                - The amplitude used for the symbols (+A or -A)
    %   Eb_N0_dB         - Range of Eb/N0 in dB
    
    % Initialize BER arrays for matched filter, hold filter, and theoretical BER
    BER_matched = zeros(size(Eb_N0_dB));  
    BER_hold = zeros(size(Eb_N0_dB));  
    theoretical_BER = zeros(size(Eb_N0_dB));  
    
    % Loop over all Eb/N0 values
    for i = 1:length(Eb_N0_dB)
        % Add noise to the transmitted signal
        [y_tx_noise, SNR_lin, N0, AWGN_scaled] = AddAWGN(y_tx, Eb, Eb_N0_dB(i));
        
        % Convolve the noisy signal with the matched filter
        y_matched = conv(y_tx_noise, p_matched); 
        
        % Convolve the noisy signal with the hold filter
        y_hold = conv(y_tx_noise, hold_filter); 
        
        % Sample the matched filter and hold filter outputs
        [y_filtered_sampled, y_hold_sampled] = plot_matched_vs_hold(1, y_matched, y_hold, samples_per_bit, bits_Num);
        
        % Calculate the BER for matched filter output
        polar_threshold = 0;
        [BER_matched(i), ~] = calculate_error_probability(bit_stream_symbols, y_filtered_sampled, polar_threshold, A);
        
        % Calculate the BER for hold filter output
        [BER_hold(i), ~] = calculate_error_probability(bit_stream_symbols, y_hold_sampled, polar_threshold, A);
        
        % Calculate theoretical BER
        theoretical_BER(i) = 0.5 * erfc(sqrt(SNR_lin));  % Theoretical BER formula
    end
end

function plot_BER_vs_EbN0(Eb_N0_dB, BER_matched, BER_hold, theoretical_BER)
    % plot_BER_vs_EbN0 - Plot the Bit Error Rate (BER) vs Eb/N0 for different filters
    %
    % Inputs:
    %   Eb_N0_dB       - Eb/N0 values in dB (vector)
    %   BER_matched    - BER for the matched filter (vector)
    %   BER_hold       - BER for the hold filter (vector)
    %   theoretical_BER - Theoretical BER (vector)
    
    % Create figure for the plot
    figure;
    
    % Plot BER vs Eb/N0 on a semi-logarithmic scale (log scale on y-axis)
    semilogy(Eb_N0_dB, BER_matched, 'b-o', 'LineWidth', 2); % Matched filter
    hold on;
    semilogy(Eb_N0_dB, BER_hold, 'r-x', 'LineWidth', 2); % Hold filter
    semilogy(Eb_N0_dB, theoretical_BER, 'k--', 'LineWidth', 2); % Theoretical BER
    hold off;
    
    % Set the labels and title
    xlabel('Eb/N0 (dB)');
    ylabel('Bit Error Rate (BER)');
    title('BER vs Eb/N0 with Matched Filter and Hold Filter');
    
    % Add legend to identify each curve
    legend('Matched Filter', 'Hold Filter', 'Theoretical BER');
    
    % Enable grid for better visibility
    grid on;
end