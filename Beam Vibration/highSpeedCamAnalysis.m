%% Position-Based Oscillation Analysis (With Velocity & Acceleration)
% Uses maxima to estimate b, omega, k and plots all motion states

clear; close all; clc;

%% Parameters
m = 0.02;         % Effective mass (kg)

% === Load position data ===
data = readmatrix('position_120Hz.csv');  % Format: [time (s), position (m)]
t_raw = data(:,1);
pos_raw = data(:,2);

% === Manual cropping ===
figure;
plot(t_raw, pos_raw, 'b');
title('Click to Crop Start/End of Oscillations');
xlabel('Time (s)'); ylabel('Position (m)');
grid on;
[crop_times, ~] = ginput(2);
crop_start = min(crop_times);
crop_end = max(crop_times);
close;

% === Apply crop ===
crop_idx = (t_raw >= crop_start) & (t_raw <= crop_end);
t = t_raw(crop_idx);
pos = pos_raw(crop_idx);

% === Derivatives: velocity and acceleration ===
v = gradient(pos, mean(diff(t)));        % Velocity (m/s)
a = gradient(v, mean(diff(t)));          % Acceleration (m/s²)

% === FFT of position to estimate frequency ===
Fs = 1 / mean(diff(t));
Y = fft(pos);
P2 = abs(Y / length(Y));
P1 = P2(1:floor(length(Y)/2)+1);
P1(2:end-1) = 2 * P1(2:end-1);
f_fft = Fs * (0:(length(P1)-1)) / length(pos);
[~, idx_fft] = max(P1);
fft_freq = f_fft(idx_fft);

% === Find maxima ===
[pks, locs] = findpeaks(pos, t, ...
    'MinPeakHeight', 0.001, ...
    'MinPeakDistance', 0.2);

if length(pks) < 8
    warning("Not enough maxima found to estimate damping.");
    return;
end

A1 = pks(1); t1 = locs(1);
A2 = pks(8); t2 = locs(8);
delta_t = t2 - t1;

if A1 <= A2
    warning("Maxima not decaying. Cannot compute damping.");
    return;
end

% === Estimate system parameters ===
b_est = (-2 * m * log(A2 / A1)) / delta_t;
omega = 7 / delta_t * 2 * pi;
k_est = m * omega^2;

% === Plot Position with Peaks ===
figure;
plot(t, pos, 'b', 'DisplayName', 'Position (m)'); hold on;
plot(locs, pks, 'k.', 'DisplayName', 'Local Maxima', 'MarkerSize', 10);
plot([t1 t2], [A1 A2], 'ro', 'DisplayName', 'Used Maxima');
title('Position with Local Maxima');
xlabel('Time (s)'); ylabel('Position (m)');
legend show; grid on;

% === Plot Velocity and Acceleration ===
figure;
subplot(2,1,1);
plot(t, v, 'm');
title('Velocity'); ylabel('m/s'); grid on;

subplot(2,1,2);
plot(t, a, 'k');
title('Acceleration'); ylabel('m/s²'); xlabel('Time (s)');
grid on;

% === Output Summary ===
fprintf('\n=== Parameter Estimates from Position (Maxima) ===\n');
fprintf('Natural Frequency (FFT): %.3f Hz\n', fft_freq);
fprintf('Damping coefficient (b): %.4f Ns/m\n', b_est);
fprintf('Angular frequency (omega): %.4f rad/s\n', omega);
fprintf('Spring constant (k): %.2f N/m\n', k_est);
