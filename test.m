clear;
clc;
close all;

%% Data Generation
A = 4; % amplitude of line code
rows_num = 500; % number of realizations
real_length = 100; % length of realization (in bits)
samples_num = 7; % number of samples per bit
samples_per_real = samples_num * real_length; % total samples per realization (700 samples)
fs = 100; % sampling frequency

Data = randi([0,1], rows_num, real_length+1);

%% Generation of polar_NRZ realizations
polar_NRZ_amplitude = (2 * Data - 1) * A;
polar_NRZ_reals = repelem(polar_NRZ_amplitude, 1, samples_num);

%% Generation of uni_polar_NRZ realizations
uni_polar_NRZ_amplitude = Data * A;
uni_polar_NRZ_reals = repelem(uni_polar_NRZ_amplitude, 1, samples_num);

%% Generation of polar_RZ realizations
polar_RZ_reals = polar_NRZ_reals;
for i = 1:rows_num
    for j = 5:7:705
        polar_RZ_reals(i, j:j+2) = 0;
    end
end

%% Random initial time shift
start_indices = randi([0, 6], rows_num, 1);
Rndm_shift_polar_NRZ = zeros(rows_num, samples_per_real);
Rndm_shift_uni_polar_NRZ = zeros(rows_num, samples_per_real);
Rndm_shift_polar_RZ = zeros(rows_num, samples_per_real);

%% Generation of random initial time shift versions of line codes
for i = 1:rows_num
    start_index = start_indices(i);
    end_index = samples_per_real + start_index;
    Rndm_shift_polar_NRZ(i, :) = polar_NRZ_reals(i, start_index+1:end_index);
    Rndm_shift_uni_polar_NRZ(i, :) = uni_polar_NRZ_reals(i, start_index+1:end_index);
    Rndm_shift_polar_RZ(i, :) = polar_RZ_reals(i, start_index+1:end_index);
end

%% Display the final line codes
figure("Name", "Line Codes");

% Uni-polar NRZ
subplot(3,1,1);
stairs(Rndm_shift_uni_polar_NRZ(1, :), 'g');
axis([0 700 -(A+1) (A+1)]);
xlabel("Time axis");
ylabel("Amplitude axis");
title("Uni-polar NRZ");

% Polar NRZ
subplot(3,1,2);
stairs(Rndm_shift_polar_NRZ(1, :), 'b');
axis([0 700 -(A+1) (A+1)]);
xlabel("Time axis");
ylabel("Amplitude axis");
title("Polar NRZ");

% Polar RZ
subplot(3,1,3);
stairs(Rndm_shift_polar_RZ(1, :), 'r');
axis([0 700 -(A+1) (A+1)]);
xlabel("Time axis");
ylabel("Amplitude axis");
title("Polar RZ");

%% Question (1) >> Statistical or ensemble mean
Col_Sum_pol_NRZ = zeros(1, samples_per_real); % initialization of sum matrix
Col_Sum_uni_pol_NRZ = zeros(1, samples_per_real); % initialization of sum matrix
Col_Sum_pol_RZ = zeros(1, samples_per_real); % initialization of sum matrix

for k = 1 : samples_per_real
    for j = 1 : rows_num
        Col_Sum_pol_NRZ(k) = Col_Sum_pol_NRZ(k) + Rndm_shift_polar_NRZ(j, k);
        Col_Sum_uni_pol_NRZ(k) = Col_Sum_uni_pol_NRZ(k) + Rndm_shift_uni_polar_NRZ(j, k);
        Col_Sum_pol_RZ(k) = Col_Sum_pol_RZ(k) + Rndm_shift_polar_RZ(j, k);
    end
end

%Averaging the vector values over the number of realizations
pol_NRZ_Stat_mean = Col_Sum_pol_NRZ .* (1 / rows_num);
uni_pol_NRZ_Stat_mean = Col_Sum_uni_pol_NRZ .* (1 / rows_num);
pol_RZ_Stat_mean = Col_Sum_pol_RZ .* (1 / rows_num);

% plotting the Statistical mean
t = 1 : samples_per_real; % samples_per_real = 7*100 = 700
figure("name", "Statistical mean plot");

% stat_mean of uni polar NRZ
subplot(3,1,1);
plot(t, uni_pol_NRZ_Stat_mean,'g');
axis([0 700 -(A+1) (A+1)]); % x-axis = 0:700 & y-axis = -5:5
xlabel("Time axis");
ylabel("Mean axis");
title("uni polar NRZ mean");

% stat_mean of polar NRZ
subplot(3,1,2);
plot(t, pol_NRZ_Stat_mean,'b');
axis([0 700 -(A+1) (A+1)]); % x-axis = 0:700 & y-axis = -5:5
xlabel("Time axis");
ylabel("Mean axis");
title("polar NRZ mean");

% stat_mean of polar RZ
subplot(3,1,3);
plot(t, pol_RZ_Stat_mean,'r');
axis([0 700 -(A+1) (A+1)]); % x-axis = 0:700 & y-axis = -5:5
xlabel("Time axis");
ylabel("Mean axis");
title("polar RZ mean");
%% Question (2) >> Statistical Autocorrelation
averages_pol_NRZ = zeros(1, 700);
averages_uni_pol_NRZ = zeros(1, 700);
averages_pol_RZ = zeros(1, 700);

for i = 0:699
    result_column_pol_NRZ = Rndm_shift_polar_NRZ(:, 1) .* Rndm_shift_polar_NRZ(:, i+1);
    averages_pol_NRZ(i+1) = mean(result_column_pol_NRZ);
end

for i = 0:699
    result_column_uni_pol_NRZ = Rndm_shift_uni_polar_NRZ(:, 1) .* Rndm_shift_uni_polar_NRZ(:, i+1);
    averages_uni_pol_NRZ(i+1) = mean(result_column_uni_pol_NRZ);
end

for i = 0:699
    result_column_pol_RZ = Rndm_shift_polar_RZ(:, 1) .* Rndm_shift_polar_RZ(:, i+1);
    averages_pol_RZ(i+1) = mean(result_column_pol_RZ);
end

t = -699 : 699;
figure("name", "Statistical Autocorrelation");

% stat_autocorr of uni polar NRZ
stat_autocorr_uni_pol_NRZ = fliplr(averages_uni_pol_NRZ);
averages_combined_uni_pol_NRZ = [stat_autocorr_uni_pol_NRZ, averages_uni_pol_NRZ(2:700)];
subplot(3,1,1);
plot(t, averages_combined_uni_pol_NRZ,'g');
xlim([-50 50]);
xlabel("Time axis");
ylabel("Autocorr axis");
title("uni polar NRZ Stat Autocorr");

% time_autocorr of polar NRZ
stat_autocorr_pol_NRZ = fliplr(averages_pol_NRZ);
averages_combined_pol_NRZ = [stat_autocorr_pol_NRZ, averages_pol_NRZ(2:700)];
subplot(3,1,2);
plot(t, averages_combined_pol_NRZ,'b');
xlim([-50 50]);
xlabel("Time axis");
ylabel("Autocorr axis");
title("polar NRZ Stat Autocorr");

% time_autocorr of polar RZ
stat_autocorr_pol_RZ = fliplr(averages_pol_RZ);
averages_combined_pol_RZ = [stat_autocorr_pol_RZ, averages_pol_RZ(2:700)];
subplot(3,1,3);
plot(t, averages_combined_pol_RZ,'r');
xlim([-50 50]);
xlabel("Time axis");
ylabel("Autocorr axis");
title("polar RZ Stat Autocorr");

%% Question (4) >> average time mean
% initilize mean vectors
polar_NRZ_time_mean = zeros(500, 1);
uni_polar_NRZ_time_mean = zeros(500, 1);
polar_RZ_time_mean = zeros(500, 1);

for i = 1: 700
    polar_NRZ_time_mean = polar_NRZ_time_mean + Rndm_shift_polar_NRZ (: , i);
    uni_polar_NRZ_time_mean = uni_polar_NRZ_time_mean + Rndm_shift_uni_polar_NRZ (: , i);
    polar_RZ_time_mean = polar_RZ_time_mean + Rndm_shift_polar_RZ (: , i);
end

%Averaging the vector values over the time
polar_NRZ_time_mean = polar_NRZ_time_mean /700;
uni_polar_NRZ_time_mean = uni_polar_NRZ_time_mean /700;
polar_RZ_time_mean = polar_RZ_time_mean /700;

t = 1 : 500;
figure("name", " average time mean plot");

% time_mean of uni polar NRZ
subplot(3,1,1);
plot(t,uni_polar_NRZ_time_mean,'g');
axis([0 500 -(A+1) (A+1)]); % x-axis = 0:700 & y-axis = -5:5
xlabel("Time axis");
ylabel("average time mean axis");
title("uni polar NRZ time mean");

% time_mean of polar NRZ
subplot(3,1,2);
plot(t,polar_NRZ_time_mean,'b');
axis([0 500 -(A+1) (A+1)]); % x-axis = 0:700 & y-axis = -5:5
xlabel("Time axis");
ylabel("average time mean axis");
title("polar NRZ time mean");

% time_mean of polar RZ
subplot(3,1,3);
plot(t,polar_RZ_time_mean,'r');
axis([0 500 -(A+1) (A+1)]); % x-axis = 0:700 & y-axis = -5:5
xlabel("Time axis");
ylabel("average time mean axis");
title("polar RZ time mean");
