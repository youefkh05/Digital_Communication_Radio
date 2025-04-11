clc; 
clear;
close all;

%using the functions from project1
clc;
clear;
close all;

%using the functions from project1
N = 50;  % Number of samples for smooth triangle
p = linspace(1, 0, N);
p = p / norm(p);  % Normalize energy


plot(p, 'LineWidth', 2);
title('Smooth Triangle Pulse Shape');
xlabel('Sample Index');
ylabel('Amplitude');
grid on;
 


% Generate bits and waveform
bits = randi([0 1], 1, num_symbols);
y_tx = generate_pam_waveform(bits, p, samples_per_symbol);

% Time vector for the output signal
t = 0:0.2:(length(y_tx)-1)*0.2;

% Plot
stairs(t, y_tx, 'LineWidth', 1.5);
title('PAM Realization (10 symbols)');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;













%-----------------------Functions----------------------------

function y_tx = generate_pam_waveform(bits, p, samples_per_symbol)
    % Converts bits to PAM impulses and convolves with pulse shape p
    % bits: binary array (0, 1)
    % p: pulse shaping vector (e.g., [5 4 3 2 1]/sqrt(55))
    % samples_per_symbol: 5 in this case
    
    % Convert bits to symbols: 0 -> -1, 1 -> +1
    symbols = 2 * bits - 1;

    % Upsample: insert (samples_per_symbol - 1) zeros between each symbol
    impulse_train = upsample(symbols, samples_per_symbol);

    % Convolve with pulse shaping function
    y_tx = conv(impulse_train, p);
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
