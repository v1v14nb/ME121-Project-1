%% Ruler Oscillation Analysis with RK4 + Parameter Uncertainty Output
% Uses AcY (col 3), estimates b/k from local minima, tail centering,
% and outputs uncertainty as ± std across trials.

clear; close all; clc;

%% Parameters
m = 0.02;         % Effective mass (kg)
dt = 0.001;       % RK4 time step
tail_window = 0.5; % seconds for mean centering

files = {'beam1.csv', 'beam2.csv', 'beam3.csv'};
num_trials = length(files);

% Storage
b_vals = zeros(1, num_trials);
omega_vals = zeros(1, num_trials);
k_vals = zeros(1, num_trials);
fft_freqs = zeros(1, num_trials);
rmse_vals = zeros(1, num_trials);

acc_all = cell(1, num_trials);
acc_sim_all = cell(1, num_trials);
t_all = cell(1, num_trials);
crop_windows = zeros(num_trials, 2);

%% === Step 1: Crop all trials ===
for i = 1:num_trials
    data = readmatrix(files{i});
    t = data(:,1) / 1000;
    acc_y = data(:,3);

    figure;
    plot(t, acc_y, 'b');
    title(sprintf('Trial %d: Click to crop start/end of oscillation', i));
    xlabel('Time (s)'); ylabel('AcY (m/s^2)'); grid on;
    [crop_times, ~] = ginput(2);
    crop_start = min(crop_times);
    crop_end = max(crop_times);
    crop_windows(i, :) = [crop_start, crop_end];
    close;
end

%% === Step 2: Analyze Trials ===
t_ref = [];

for i = 1:num_trials
    data = readmatrix(files{i});
    t = data(:,1) / 1000;
    acc_y = data(:,3);

    crop_start = crop_windows(i, 1);
    crop_end = crop_windows(i, 2);
    start_idx = find(t >= crop_start, 1);
    end_idx = find(t <= crop_end, 1, 'last');

    if isempty(start_idx) || isempty(end_idx) || end_idx <= start_idx
        warning("Trial %d: Invalid crop. Skipping.", i);
        continue;
    end

    t = t(start_idx:end_idx);
    acc_y = acc_y(start_idx:end_idx);

    % Tail-based mean subtraction
    tail_idx = t > (t(end) - tail_window);
    acc_y = acc_y - mean(acc_y(tail_idx));

    if isempty(t_ref)
        t_ref = t(1):dt:t(end);
    end

    % FFT for natural frequency
    Fs = 1 / mean(diff(t));
    Y = fft(acc_y);
    P2 = abs(Y / length(Y));
    P1 = P2(1:floor(length(Y)/2)+1);
    P1(2:end-1) = 2 * P1(2:end-1);
    f_fft = Fs * (0:(length(P1)-1)) / length(acc_y);
    [~, idx_fft] = max(P1);
    fft_freqs(i) = f_fft(idx_fft);

    % Find local minima
    [pks, locs] = findpeaks(-acc_y, t, ...
        'MinPeakHeight', 0.01, ...
        'MinPeakDistance', 0.2);
    pks = abs(pks);

    if length(pks) < 8
        warning("Trial %d: Not enough minima.", i);
        continue;
    end

    A1 = pks(1); t1 = locs(1);
    A2 = pks(8); t2 = locs(8);
    delta_t = t2 - t1;

    if A1 <= A2
        warning("Trial %d: Minima not decaying. Skipping.", i);
        continue;
    end

    % Estimate system parameters
    b_est = (-2 * m * log(A2 / A1)) / delta_t;
    omega = 7 / delta_t * 2 * pi;
    k_est = m * omega^2;

    if b_est <= 0
        warning("Trial %d: Non-positive damping. Skipping.", i);
        continue;
    end

    b_vals(i) = b_est;
    omega_vals(i) = omega;
    k_vals(i) = k_est;

    % Estimate initial displacement from A1
    y0 = [-A1 * m / k_est; 0];

    % RK4 simulation
    f = @(t, y) [y(2); (1/m)*(-k_est*y(1) - b_est*y(2))];
    y_sim = rk4_solver(f, t_ref, y0, dt);
    acc_sim = gradient(gradient(y_sim(:,1), dt), dt);

    % Store for overlay and error
    acc_all{i} = acc_y;
    acc_sim_all{i} = acc_sim;
    t_all{i} = t;

    common_len = min(length(t_ref), length(acc_y));
    rmse_vals(i) = sqrt(mean((acc_sim(1:common_len) - acc_y(1:common_len)).^2));

    % Trial figure
    figure;
    plot(t, acc_y, 'b', 'DisplayName', 'Measured AcY'); hold on;
    plot(t_ref, acc_sim, 'r--', 'DisplayName', 'RK4 Simulated');
    plot(locs, pks, 'k.', 'DisplayName', 'Local Minima', 'MarkerSize', 10);
    plot([t1 t2], [A1 A2], 'ro', 'DisplayName', 'Used Minima');
    title(sprintf('Trial %d: Measured vs Simulated (AcY)', i));
    xlabel('Time (s)'); ylabel('Acceleration (m/s^2)');
    legend show; grid on;
end

%% === Step 3: Summary Plot ===
figure;
hold on;
valid_trials = find(~cellfun(@isempty, acc_sim_all));
for i = valid_trials
    plot(t_all{i}, acc_all{i}, ...
        'Color', [0.2 0.2 1 0.4], ...
        'DisplayName', sprintf('Trial %d - Measured', i));
    plot(t_ref, acc_sim_all{i}, ...
        'Color', [1 0 0 0.4], ...
        'LineStyle', '--', ...
        'DisplayName', sprintf('Trial %d - Simulated', i));
end
xlabel('Time (s)');
ylabel('Acceleration (m/s^2)');
title('All Trials Overlay: Measured and Simulated');
legend show; grid on;

%% === Step 4: Output Parameter Estimates with Uncertainty ===
valid_b = b_vals(b_vals > 0);
valid_omega = omega_vals(omega_vals > 0);
valid_k = k_vals(k_vals > 0);
valid_rmse = rmse_vals(rmse_vals > 0);
valid_fft = fft_freqs(fft_freqs > 0);

fprintf('\n=== Parameter Estimates with Uncertainty (n = %d trials) ===\n', length(valid_b));

fprintf('Natural Frequency (FFT): %.3f Hz\n', mean(valid_fft));
fprintf('  ± %.3f Hz (standard deviation)\n', std(valid_fft));

fprintf('Damping coefficient (b): %.4f Ns/m\n', mean(valid_b));
fprintf('  ± %.4f Ns/m (standard deviation)\n', std(valid_b));

fprintf('Angular frequency (omega): %.4f rad/s\n', mean(valid_omega));
fprintf('  ± %.4f rad/s (standard deviation)\n', std(valid_omega));

fprintf('Spring constant (k): %.2f N/m\n', mean(valid_k));
fprintf('  ± %.2f N/m (standard deviation)\n', std(valid_k));

fprintf('Model RMSE: %.3f m/s²\n', mean(valid_rmse));
fprintf('  ± %.3f m/s² (standard deviation)\n', std(valid_rmse));

%% === RK4 Solver ===
function y = rk4_solver(f, tspan, y0, dt)
    N = length(tspan);
    y = zeros(N, length(y0));
    y(1,:) = y0';
    for i = 1:N-1
        t = tspan(i);
        y_i = y(i,:)';
        k1 = f(t, y_i);
        k2 = f(t + dt/2, y_i + dt/2 * k1);
        k3 = f(t + dt/2, y_i + dt/2 * k2);
        k4 = f(t + dt, y_i + dt * k3);
        y(i+1,:) = (y_i + dt/6 * (k1 + 2*k2 + 2*k3 + k4))';
    end
end
