clear
close all

% --- USER PARAMETERS ---
mass = 0.30;           % Mass in kg
numTrials = 5;          % Number of trials per orientation

% Store results
results = struct('Z', [], 'Y', []);

% --- ANALYZE HORIZONTAL (Z) TRIALS ---
fprintf('\n=== HORIZONTAL (Z axis) TRIALS ===\n');

for trial = 1:numTrials
    filename = sprintf('Z%d.csv', trial);
    fprintf('\nProcessing %s...\n', filename);
    data = readtable(filename);
    time = data{:, 1};
    accZ = data{:, end-1};  % assuming Z is last column

    % Plot and select bounds by clicking
    figure;
    plot(time, accZ); title(sprintf('trialZ%d - Click to bound', trial));
    xlabel('Time (s)'); ylabel('Z Acceleration (m/s^2)'); grid on;
    [xBounds, ~] = ginput(2);
    t_start = min(xBounds); t_end = max(xBounds);

    idx = (time >= t_start) & (time <= t_end);
    t = time(idx); a = accZ(idx);

    % Peak Detection
    [peakVals, peakTimes] = findpeaks(a, t, 'MinPeakProminence', 0.2);
    if numel(peakTimes) < 2, fprintf('Not enough peaks. Skipping.\n'); continue; end

    % Compute Oscillation Parameters
    
    % Damped Natural Frequency w_d: intervals between peaks
    T = mean(diff(peakTimes));
    wd = 2 * pi / T;

    % Logarithmic Decrement delta: ln of ratio of successive peak amplitudes
    delta = mean(log(peakVals(1:end-1) ./ peakVals(2:end)));
    
    % Damping Ratio zeta:
    zeta = delta / sqrt(4 * pi^2 + delta^2);

    % Undamped natural frequency w_n: 
    wn = wd / sqrt(1 - zeta^2);

    % Spring constant k:
    k = mass * wn^2;

    % Damping coeff c:
    c = 2 * mass * zeta * wn;

    results.Z(end+1,:) = [wd, wn, zeta, delta, k, c];

    % Plot
    figure; plot(t, a); hold on;
    plot(peakTimes, peakVals, 'ro', 'MarkerFaceColor', 'r');
    title(sprintf('trialZ%d - Peaks', trial)); grid on;
end

% --- ANALYZE VERTICAL (Y) TRIALS ---
fprintf('\n=== VERTICAL (Y axis) TRIALS ===\n');

for trial = 1:numTrials
    filename = sprintf('Y%d.csv', trial);
    fprintf('\nProcessing %s...\n', filename);
    data = readtable(filename);
    time = data{:, 1};
    accY = data{:, end-2};  % assuming Y is second-to-last column

    % Plot and select bounds
    figure;
    plot(time, accY); title(sprintf('trialY%d - Click to bound', trial));
    xlabel('Time (s)'); ylabel('Y Acceleration (m/s^2)'); grid on;
    [xBounds, ~] = ginput(2);
    t_start = min(xBounds); t_end = max(xBounds);

    idx = (time >= t_start) & (time <= t_end);
    t = time(idx); a = accY(idx);

    % Peak Detection
    [peakVals, peakTimes] = findpeaks(a, t, 'MinPeakProminence', 0.2);
    if numel(peakTimes) < 2, fprintf('Not enough peaks. Skipping.\n'); continue; end

    T = mean(diff(peakTimes));
    wd = 2 * pi / T;
    delta = mean(log(peakVals(1:end-1) ./ peakVals(2:end)));
    zeta = delta / sqrt(4 * pi^2 + delta^2);
    wn = wd / sqrt(1 - zeta^2);
    k = mass * wn^2;
    c = 2 * mass * zeta * wn;

    results.Y(end+1,:) = [wd, wn, zeta, delta, k, c];

    % Plot
    figure; plot(t, a); hold on;
    plot(peakTimes, peakVals, 'ro', 'MarkerFaceColor', 'r');
    title(sprintf('trialY%d - Peaks', trial)); grid on;
end

% --- Display Summary Stats ---
labels = {'wd (rad/s)', 'wn (rad/s)', 'zeta', 'delta', 'k (N/m)', 'c (Ns/m)'};

fprintf('\n=== AVERAGE RESULTS ===\n');
for ax = ["Z", "Y"]
    group = results.(ax);
    if isempty(group), continue; end
    means = mean(group, 1);
    stds = std(group, 0, 1);
    fprintf('\n%s-axis trials:\n', ax);
    for i = 1:length(labels)
        fprintf('  %-10s: %.4f ± %.4f\n', labels{i}, means(i), stds(i));
    end
end
