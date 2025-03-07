clc; clear; close all;

% Parameters
A = 4; % Amplitude
N_realizations = 500; % Number of waveforms (ensemble size)
num_bits = 100; % Bits per waveform
bit_duration = 70e-3; % 70 ms per bit
dac_interval = 10e-3; % DAC updates every 10 ms
samples_per_bit = bit_duration / dac_interval; % 7 samples per bit
total_time = num_bits * bit_duration; % Total waveform duration
t = 0:dac_interval:(total_time - dac_interval); % Time vector

% Initialize arrays
Unipolar_ensemble = zeros(N_realizations, length(t));
PolarNRZ_ensemble = zeros(N_realizations, length(t));
RZ_ensemble = zeros(N_realizations, length(t));

for i = 1:N_realizations
    % Generate a random bit sequence
    Data = randi([0, 1], 1, num_bits);

    % Unipolar NRZ: 0 → 0V, 1 → A
    Unipolar = Data * A;
    Unipolar_expanded = repmat(Unipolar, samples_per_bit, 1);
    Unipolar_ensemble(i, :) = reshape(Unipolar_expanded, [], 1);

    % Polar NRZ: 0 → -A, 1 → +A
    PolarNRZ = ((2 * Data) - 1) * A;
    PolarNRZ_expanded = repmat(PolarNRZ, samples_per_bit, 1);
    PolarNRZ_ensemble(i, :) = reshape(PolarNRZ_expanded, [], 1);

    % Return-to-Zero (RZ): 0 → 0V, 1 → A but returns to 0 in half-bit time
    RZ = zeros(1, num_bits * 2);
    RZ(1:2:end) = Data * A; % Set amplitude for first half of each bit period
    RZ_expanded = repmat(RZ, samples_per_bit / 2, 1);
    RZ_ensemble(i, :) = reshape(RZ_expanded, [], 1);
end

%% Compute the Statistical Mean
mean_unipolar = mean(Unipolar_ensemble, 1);
mean_polarNRZ = mean(PolarNRZ_ensemble, 1);
mean_RZ = mean(RZ_ensemble, 1);

%% Check Stationarity (Plot Mean Over Time)
figure;
subplot(3,1,1); plot(t, mean_unipolar, 'r'); title('Mean of Unipolar NRZ');
subplot(3,1,2); plot(t, mean_polarNRZ, 'g'); title('Mean of Polar NRZ');
subplot(3,1,3); plot(t, mean_RZ, 'b'); title('Mean of RZ');
xlabel('Time (s)'); grid on;

%% Compute Ensemble Autocorrelation R_x(τ)
tau_values = -100:100;
Rx_unipolar = zeros(size(tau_values));
Rx_polarNRZ = zeros(size(tau_values));
Rx_RZ = zeros(size(tau_values));

for tau = tau_values
    shift = abs(tau); % Shift in samples
    for i = 1:N_realizations
        Rx_unipolar(tau + 101) = Rx_unipolar(tau + 101) + ...
            mean(Unipolar_ensemble(i, 1:end-shift) .* Unipolar_ensemble(i, shift+1:end));
        Rx_polarNRZ(tau + 101) = Rx_polarNRZ(tau + 101) + ...
            mean(PolarNRZ_ensemble(i, 1:end-shift) .* PolarNRZ_ensemble(i, shift+1:end));
        Rx_RZ(tau + 101) = Rx_RZ(tau + 101) + ...
            mean(RZ_ensemble(i, 1:end-shift) .* RZ_ensemble(i, shift+1:end));
    end
end

Rx_unipolar = Rx_unipolar / N_realizations;
Rx_polarNRZ = Rx_polarNRZ / N_realizations;
Rx_RZ = Rx_RZ / N_realizations;

% Plot Autocorrelation
figure;
subplot(3,1,1); plot(tau_values, Rx_unipolar, 'r'); title('Autocorrelation of Unipolar NRZ');
subplot(3,1,2); plot(tau_values, Rx_polarNRZ, 'g'); title('Autocorrelation of Polar NRZ');
subplot(3,1,3); plot(tau_values, Rx_RZ, 'b'); title('Autocorrelation of RZ');
xlabel('τ (Lag Samples)'); grid on;

%% Compute Time Mean and Autocorrelation for One Waveform
one_waveform = PolarNRZ_ensemble(1, :);
time_mean = mean(one_waveform);
time_autocorr = xcorr(one_waveform, 'unbiased');

% Plot Time Mean & Autocorrelation
figure;
subplot(2,1,1); plot(t, one_waveform, 'k'); title('One Realization of Polar NRZ');
subplot(2,1,2); plot(time_autocorr, 'm'); title('Time Autocorrelation of One Waveform');
grid on;

%% Check for Ergodicity
% A process is ergodic if time averages and ensemble averages are approximately equal
ergodicity_check = abs(time_mean - mean(mean_polarNRZ)) < 0.05;
disp(['Is the process Ergodic? ', num2str(ergodicity_check)]);

%% Compute Bandwidth using FFT
Fs = 1/dac_interval; % Sampling frequency
f_axis = linspace(-Fs/2, Fs/2, length(t));

% Compute FFT for each line code
fft_unipolar = abs(fftshift(fft(mean_unipolar)));
fft_polarNRZ = abs(fftshift(fft(mean_polarNRZ)));
fft_RZ = abs(fftshift(fft(mean_RZ)));

% Plot Frequency Domain Analysis
figure;
subplot(3,1,1); plot(f_axis, fft_unipolar, 'r'); title('Spectrum of Unipolar NRZ');
subplot(3,1,2); plot(f_axis, fft_polarNRZ, 'g'); title('Spectrum of Polar NRZ');
subplot(3,1,3); plot(f_axis, fft_RZ, 'b'); title('Spectrum of RZ');
xlabel('Frequency (Hz)'); grid on;

