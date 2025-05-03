% --- Full FRF Analysis: Unmodified vs Modified Plate Comparison ---
clear; clc; close all;

% --- PARAMETERS ---
Npoints = 5;    % 5 tap points + 1 background noise
Ntrials = 3;    % 3 trials per point
Npeaks = 6;     % Number of peaks to find
freq_low_limit = 400;
freq_high_limit = 1900;

% --- STORAGE ---
all_f_common = cell(Npoints,1); all_Y_matrix = cell(Npoints,1); all_peakFreqs_matrix = cell(Npoints,1);
all_f_common_mass = cell(Npoints,1); all_Y_matrix_mass = cell(Npoints,1); all_peakFreqs_matrix_mass = cell(Npoints,1);

all_avg_peaks = zeros(Npoints, Npeaks); all_std_peaks = zeros(Npoints, Npeaks);
all_new_avg_peaks = zeros(Npoints, Npeaks); all_new_std_peaks = zeros(Npoints, Npeaks);
all_avg_peaks_mass = zeros(Npoints, Npeaks); all_std_peaks_mass = zeros(Npoints, Npeaks);
all_new_avg_peaks_mass = zeros(Npoints, Npeaks); all_new_std_peaks_mass = zeros(Npoints, Npeaks);

final_unmod_FFT = {}; final_mod_FFT = {}; final_unmod_f = {}; final_mod_f = {};

%% --- UNMODIFIED DATA ---
figure; tiledlayout(Npoints,1,'TileSpacing','compact');
for pt = 1:Npoints
    filenames = arrayfun(@(tr) sprintf('%d%d.wav', pt, tr), 1:Ntrials, 'UniformOutput', false);
    [f_common, Y_matrix, peakFreqs_matrix, ~, ~] = processTapFiles(filenames, Npeaks, freq_low_limit, freq_high_limit);

    all_f_common{pt} = f_common; all_Y_matrix{pt} = Y_matrix; all_peakFreqs_matrix{pt} = peakFreqs_matrix;

    for k = 1:Ntrials
        final_unmod_FFT{end+1} = Y_matrix(k,:);
        final_unmod_f{end+1} = f_common;
    end

    nexttile; hold on;
    colors = lines(Ntrials);
    for k = 1:Ntrials
plot(f_common, Y_matrix(k,:), 'Color', colors(k,:), 'DisplayName', sprintf('Trial %d',k));
    end
    title(sprintf('Unmodified Point %d', pt));
    xlabel('Frequency (Hz)'); ylabel('Normalized Mag'); xlim([freq_low_limit freq_high_limit]); grid on; hold off;

    avg_peaks = mean(peakFreqs_matrix,1,'omitnan');
    std_peaks = std(peakFreqs_matrix,0,1,'omitnan');
    all_avg_peaks(pt,:) = avg_peaks;
    all_std_peaks(pt,:) = std_peaks;
end

% Pre-manual Uncertainty Plot (Unmod)
figure; hold on;
for pt = 1:Npoints-1
    errorbar(1:Npeaks, all_avg_peaks(pt,:), all_std_peaks(pt,:), '-o', 'DisplayName', sprintf('Point %d',pt));
end
xlabel('Peak Number'); ylabel('Frequency (Hz)');
title('Unmodified Pre-Manual Peak Uncertainty');
legend('show'); grid on; hold off;

%% --- MODIFIED DATA ---
figure; tiledlayout(Npoints,1,'TileSpacing','compact');
for pt = 1:Npoints
    filenames = arrayfun(@(tr) sprintf('mass%d%d.wav', pt, tr), 1:Ntrials, 'UniformOutput', false);
    [f_common, Y_matrix, peakFreqs_matrix, ~, ~] = processTapFiles(filenames, Npeaks, freq_low_limit, freq_high_limit);

    all_f_common_mass{pt} = f_common; all_Y_matrix_mass{pt} = Y_matrix; all_peakFreqs_matrix_mass{pt} = peakFreqs_matrix;

    for k = 1:Ntrials
        final_mod_FFT{end+1} = Y_matrix(k,:);
        final_mod_f{end+1} = f_common;
    end

    nexttile; hold on;
    colors = lines(Ntrials);
    for k = 1:Ntrials
        plot(f_common, Y_matrix(k,:), 'Color', colors(k,:), 'DisplayName', sprintf('Trial %d',k));

    end
    title(sprintf('Modified Point %d', pt));
    xlabel('Frequency (Hz)'); ylabel('Normalized Mag'); xlim([freq_low_limit freq_high_limit]); grid on; hold off;

    avg_peaks = mean(peakFreqs_matrix,1,'omitnan');
    std_peaks = std(peakFreqs_matrix,0,1,'omitnan');
    all_avg_peaks_mass(pt,:) = avg_peaks;
    all_std_peaks_mass(pt,:) = std_peaks;
end

% Pre-manual Uncertainty Plot (Mod)
figure; hold on;
for pt = 1:Npoints-1
    errorbar(1:Npeaks, all_avg_peaks_mass(pt,:), all_std_peaks_mass(pt,:), '-o', 'DisplayName', sprintf('Point %d',pt));
end
xlabel('Peak Number'); ylabel('Frequency (Hz)');
title('Modified Pre-Manual Peak Uncertainty');
legend('show'); grid on; hold off;
%% --- Manual Peak Selection: Unmodified (Flexible) ---
fprintf('\n>>> UNMODIFIED: Click LOWER and UPPER bounds for each peak, in pairs.\n');
fprintf('    Press ENTER when done selecting for a point to move to the next one.\n');

region_bounds = cell(Npoints, 1);

figure(1);
for pt = 1:Npoints
    nexttile(pt);
    title(sprintf('Unmodified Point %d — Select peak bounds (press ENTER to finish)', pt));
    hold on;
    regions = [];
    done = false;

    while ~done
        [x,~,button] = ginput(1);
        if isempty(x)  % ENTER key ends input
            done = true;
            break;
        end

        [x2,~,button2] = ginput(1);
        if isempty(x2)
            warning('Second point was not selected — skipping last range.');
            break;
        end

        pair = sort([x, x2]);
        regions = [regions; pair];
        plot(pair, [1 1], 'k--', 'LineWidth', 1.2); % Visual confirmation
    end

    region_bounds{pt} = regions;
end

% --- Peak Picking (Unmodified) ---
for pt = 1:Npoints
    f_common = all_f_common{pt};
    Y_matrix = all_Y_matrix{pt};
    regions = region_bounds{pt};
    Npeaks_actual = size(regions,1);
    peakFreqs_matrix_new = NaN(Ntrials, Npeaks_actual);

    for k = 1:Ntrials
        for p = 1:Npeaks_actual
            lowB = regions(p,1);
            highB = regions(p,2);
            idx_region = (f_common >= lowB) & (f_common <= highB);
            if any(idx_region)
                [~, idx_local] = max(Y_matrix(k,idx_region));
                idx_global = find(idx_region);
                true_idx = idx_global(idx_local);
                peakFreqs_matrix_new(k,p) = f_common(true_idx);
            end
        end
    end

    all_new_avg_peaks(pt,1:Npeaks_actual) = mean(peakFreqs_matrix_new,1,'omitnan');
    all_new_std_peaks(pt,1:Npeaks_actual) = std(peakFreqs_matrix_new,0,1,'omitnan');
end


% Post-manual Uncertainty Plot (Unmod)
figure; hold on;


for pt = 1:Npoints-1
    nActualPeaks = sum(~isnan(all_new_avg_peaks(pt,:)));
    errorbar(1:nActualPeaks, ...
             all_new_avg_peaks(pt,1:nActualPeaks), ...
             all_new_std_peaks(pt,1:nActualPeaks), ...
             '-o', 'DisplayName', sprintf('Point %d',pt));
end



xlabel('Peak Number'); ylabel('Frequency (Hz)');
title('Unmodified Post-Manual Peak Uncertainty');
legend('show'); grid on; hold off;
%% --- Manual Peak Selection: Modified (Flexible) ---
fprintf('\n>>> MODIFIED: Click LOWER and UPPER bounds for each peak, in pairs.\n');
fprintf('    Press ENTER when done selecting for a point to move to the next one.\n');

region_bounds_mass = cell(Npoints, 1);

figure(3);
for pt = 1:Npoints
    nexttile(pt);
    title(sprintf('Modified Point %d — Select peak bounds (press ENTER to finish)', pt));
    hold on;
    regions = [];
    done = false;

    while ~done
        [x,~,button] = ginput(1);
        if isempty(x)  % ENTER key ends input
            done = true;
            break;
        end

        [x2,~,button2] = ginput(1);
        if isempty(x2)
            warning('Second point was not selected — skipping last range.');
            break;
        end

        pair = sort([x, x2]);
        regions = [regions; pair];
        plot(pair, [1 1], 'k--', 'LineWidth', 1.2); % Visual confirmation
    end

    region_bounds_mass{pt} = regions;
end

% --- Peak Picking (Modified) ---
for pt = 1:Npoints
    f_common = all_f_common_mass{pt};
    Y_matrix = all_Y_matrix_mass{pt};
    regions = region_bounds_mass{pt};
    Npeaks_actual = size(regions,1);
    peakFreqs_matrix_new = NaN(Ntrials, Npeaks_actual);

    for k = 1:Ntrials
        for p = 1:Npeaks_actual
            lowB = regions(p,1);
            highB = regions(p,2);
            idx_region = (f_common >= lowB) & (f_common <= highB);
            if any(idx_region)
                [~, idx_local] = max(Y_matrix(k,idx_region));
                idx_global = find(idx_region);
                true_idx = idx_global(idx_local);
                peakFreqs_matrix_new(k,p) = f_common(true_idx);
            end
        end
    end

    all_new_avg_peaks_mass(pt,1:Npeaks_actual) = mean(peakFreqs_matrix_new,1,'omitnan');
    all_new_std_peaks_mass(pt,1:Npeaks_actual) = std(peakFreqs_matrix_new,0,1,'omitnan');
end

% Post-manual Uncertainty Plot (Mod)
figure; hold on;


for pt = 1:Npoints-1
    nActualPeaks = sum(~isnan(all_new_avg_peaks_mass(pt,:)));
    errorbar(1:nActualPeaks, ...
             all_new_avg_peaks_mass(pt,1:nActualPeaks), ...
             all_new_std_peaks_mass(pt,1:nActualPeaks), ...
             '-o', 'DisplayName', sprintf('Point %d',pt));
end


xlabel('Peak Number'); ylabel('Frequency (Hz)');
title('Modified Post-Manual Peak Uncertainty');
legend('show'); grid on; hold off;

%% --- Final FFT Overlay: Frequency vs Normalized Magnitude ---
figure; hold on;
for idx = 1:length(final_unmod_FFT)
    plot(final_unmod_f{idx}, final_unmod_FFT{idx}, 'Color', [0 0 1 0.3], 'LineWidth', 1.2);
end
for idx = 1:length(final_mod_FFT)
    plot(final_mod_f{idx}, final_mod_FFT{idx}, 'Color', [1 0 0 0.3], 'LineWidth', 1.2);
end
xlabel('Frequency (Hz)'); ylabel('Normalized Magnitude');
title('Overlay of All FFTs (Blue = Unmodified, Red = Modified)');
xlim([freq_low_limit freq_high_limit]); grid on; hold off;


figure;
hold on;
for pt = 1:Npoints
    n_unmod = sum(~isnan(all_new_avg_peaks(pt,:)));
    n_mod   = sum(~isnan(all_new_avg_peaks_mass(pt,:)));

    errorbar(1:n_unmod, ...
             all_new_avg_peaks(pt,1:n_unmod), ...
             all_new_std_peaks(pt,1:n_unmod), ...
             '-o', 'Color', [0 0 1], 'DisplayName', sprintf('Point %d (Unmod)', pt));  % Blue

    errorbar(1:n_mod, ...
             all_new_avg_peaks_mass(pt,1:n_mod), ...
             all_new_std_peaks_mass(pt,1:n_mod), ...
             '-s', 'Color', [1 0 0], 'DisplayName', sprintf('Point %d (Mod)', pt));   % Red
end
xlabel('Peak Number');
ylabel('Frequency (Hz)');
title('Tap Point Resonant Frequencies (Post-Manual Peaks)');
grid on;
legend('show', 'Location', 'eastoutside');
hold off;


%% --- Function for Tap Processing ---
function [f_common, Y_matrix, peakFreqs_matrix, localPeakFreqs, localPeakMags] = processTapFiles(filenames, Npeaks, freq_low_limit, freq_high_limit)
    numFiles = length(filenames);
    signals = {}; lengths = zeros(1,numFiles);

    for k = 1:numFiles
        [signal, fs] = audioread(filenames{k});
        if size(signal,2) > 1
            signal = signal(:,1); % take first channel if stereo
        end
        % --- Normalize raw time-domain signal by RMS energy ---
        rms_val = rms(signal);
        if rms_val > 0
            signal = signal / rms_val;
        end
        signals{k} = signal;
        lengths(k) = length(signal);
    end

    maxLength = max(lengths);
    Y_matrix = [];
    peakFreqs_matrix = zeros(numFiles, Npeaks);
    localPeakFreqs = cell(numFiles,1);
    localPeakMags = cell(numFiles,1);

    for k = 1:numFiles
        % --- Zero-pad to match max length ---
        signal_padded = zeros(maxLength,1);
        signal_padded(1:length(signals{k})) = signals{k};

        % --- FFT ---
        N = length(signal_padded);
        Y = fft(signal_padded);
        f = (0:N-1)*(fs/N);
        halfIdx = 1:floor(N/2);
        Y_half = abs(Y(halfIdx));
        f_half = f(halfIdx);

        % --- Frequency window ---
        freqIdx = (f_half >= freq_low_limit) & (f_half <= freq_high_limit);
        f_half = f_half(freqIdx);
        Y_half = Y_half(freqIdx);

        % --- Save FFT row ---
        Y_matrix(k,:) = Y_half;

        % --- Peak detection ---
        peak_height_thresh = 0.002;
        if max(Y_half) < peak_height_thresh
            peak_height_thresh = 0;  % avoid invalid threshold
        end
        [peakValsAll, peakLocsAll] = findpeaks(Y_half, f_half, ...
            'MinPeakProminence', 0.005, ...
            'MinPeakDistance', 30, ...
            'MinPeakHeight', peak_height_thresh);

        localPeakFreqs{k} = peakLocsAll;
        localPeakMags{k} = peakValsAll;

        if length(peakValsAll) >= Npeaks
            [~, idx] = maxk(peakValsAll, Npeaks);
            peakLocs = peakLocsAll(idx);
        else
            peakLocs = nan(1,Npeaks);
            peakLocs(1:length(peakLocsAll)) = peakLocsAll;
        end
        peakFreqs_matrix(k,:) = sort(peakLocs);
    end

    f_common = f_half;
end
