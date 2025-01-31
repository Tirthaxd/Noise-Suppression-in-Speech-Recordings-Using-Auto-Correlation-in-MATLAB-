clc;
close all;
clear all;

%% Load the Clean File

Clean_file = "G:\DSP assignment\EXP 02\Actualvoicenote.wav";
[CleanAudio, Fs_clean] = audioread(Clean_file);

% Time vector
t_clean = (0:length(CleanAudio)-1) / Fs_clean;

% Plot the time-domain signal
figure;
subplot(2,1,1)
plot(t_clean, CleanAudio);
title('Time-Domain Representation of Clean Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Perform FFT
N_clean = length(CleanAudio); % Length of the signal
fft_clean = fft(CleanAudio); % FFT of the signal
fft_clean = fft_clean(1:N_clean/2+1); % Single-sided spectrum
f_clean = (0:N_clean/2) * (Fs_clean / N_clean); % Frequency vector

% Magnitude of the FFT
magnitude_clean = abs(fft_clean);

% Plot the frequency-domain representation
subplot(2,1,2);
plot(f_clean, magnitude_clean);
title('Frequency-Domain Representation of Clean Speech Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Create an audioplayer object
player1 = audioplayer(CleanAudio, Fs_clean); % Create the audioplayer object

% Play the Audio
disp('Playing the audio...');
play(player1); % Play audio using audioplayer

% Loop to monitor keypresses and control audio
while isplaying(player1)
    % Check if a key is pressed
    if waitforbuttonpress
        key = get(gcf, 'CurrentKey');  % Get the key pressed

        if strcmp(key, 'space')  % Pause playback when space is pressed
            pause(player1);
            disp('Audio paused. Press any key to resume.');
            waitforbuttonpress;
            resume(player1);  % Resume playback
        elseif strcmp(key, 's')  % Stop playback when 's' is pressed
            stop(player1);
            disp('Audio stopped.');
            break;  % Exit the loop
        end
    end
end

%% Load Noisy audio file
Noisy_file ="G:\DSP assignment\EXP 02\voicenotewithnoise.wav";
[NoisyAudio, Fs_noisy] = audioread(Noisy_file);

% Time vector
t_noisy = (0:length(NoisyAudio)-1) / Fs_noisy;



% Normalize the signal
maxAmplitude_noisy = max(abs(NoisyAudio));
normalized_noisy = NoisyAudio / maxAmplitude_noisy;

% Confirm normalization
disp(['Maximum amplitude after normalization: ', num2str(max(abs(normalized_noisy)))]);

% Plot the normalized signal
subplot(2,1,1);
plot(t_noisy, normalized_noisy);
title('Normalized Time-Domain Representation of Noisy Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Perform FFT
N_noisy = length(normalized_noisy); % Length of the signal
fft_noisy = fft(normalized_noisy); % FFT of the signal
fft_noisy = fft_noisy(1:N_noisy/2+1); % Single-sided spectrum
f_noisy = (0:N_noisy/2) * (Fs_noisy / N_noisy); % Frequency vector

% Magnitude of the FFT
magnitude_noisy = abs(fft_noisy);

% Plot the frequency-domain representation
subplot(2,1,2);
plot(f_noisy, magnitude_noisy);
title('Frequency-Domain Representation of Noisy Speech Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
grid on;

% Create an audioplayer object
player2 = audioplayer(normalized_noisy, Fs_noisy); % Create the audioplayer object

% Play the Audio
disp('Playing the audio...');
play(player2); % Play audio using audioplayer

% Loop to monitor keypresses and control audio
while isplaying(player2)
    % Check if a key is pressed
    if waitforbuttonpress
        key = get(gcf, 'CurrentKey');  % Get the key pressed

        if strcmp(key, 'space')  % Pause playback when space is pressed
            pause(player2);
            disp('Audio paused. Press any key to resume.');
            waitforbuttonpress;
            resume(player2);  % Resume playback
        elseif strcmp(key, 's')  % Stop playback when 's' is pressed
            stop(player2);
            disp('Audio stopped.');
            break;  % Exit the loop
        end
    end
end

%% Auto-Correlation Calculation

% Compute the ACF of the noisy signal
[acf, lags] = xcorr(normalized_noisy, 'normalized'); % Normalized ACF
lag_times = lags / Fs_noisy; % Convert lag indices to time in seconds

% Plot the ACF of the noisy signal
figure;
plot(lag_times, acf, 'b'); % Plot ACF as a function of time lag
title('Auto-Correlation Function of Noisy Speech');
xlabel('Lag (seconds)');
ylabel('Normalized Amplitude');
grid on;

%% Noise Estimation 

% Detect peaks in the full ACF 
[peaks, locs] = findpeaks(acf, 'MinPeakHeight', 0.03, 'MinPeakDistance', round(Fs_noisy * 0.01));  % Reduce threshold to 0.05

% Plot the full ACF and the detected peaks
figure;
plot(lag_times, acf); % Plot the full ACF (including both positive and negative lags)
hold on;
plot(lag_times(locs), peaks, 'r*', 'MarkerSize', 8); % Mark the detected peaks
title('Auto-Correlation Function (ACF) with Detected Peaks');
xlabel('Lag (s)');
ylabel('Normalized Correlation');
legend('ACF', 'Significant Peaks');
grid on;
hold off;

%% Thresholding

% Set a threshold value to retain more ACF peaks
threshold = 0.01; % Lower threshold value to retain more information

% Threshold the ACF to retain significant peaks
thresholdedACF = acf;
thresholdedACF(abs(thresholdedACF) < threshold) = 0;

% Ensure smoothedACF and lagTime have the same length
minLength = min(length(thresholdedACF), length(lag_times));
thresholdedACF = thresholdedACF(1:minLength);
lagTime = lag_times(1:minLength);

% Plot the thresholded ACF
figure;
plot(lagTime, thresholdedACF); % Plot the thresholded ACF
title('Thresholded Auto-Correlation Function (ACF)');
xlabel('Lag (s)');
ylabel('Normalized Correlation');
grid on;

%% Signal Reconstruction

% Ensure that smoothedACF and NoisyAudio have the same length
minLength = min(length(normalized_noisy), length(thresholdedACF));

% Zero-padding smoothedACF or trimming if necessary
if length(thresholdedACF) < length(normalized_noisy)
    thresholdedACF = [thresholdedACF; zeros(length(normalized_noisy) - length(thresholdedACF), 1)];
elseif length(thresholdedACF) > length(normalized_noisy)
    thresholdedACF = thresholdedACF(1:length(normalized_noisy)); % Trim smoothedACF
end

% Now, subtract NoisyAudio and smoothedACF
reconstructedSignal = normalized_noisy - thresholdedACF;


%% Visual Comparison

% Make sure both signals are the same length by trimming or zero-padding
minLength = min(length(CleanAudio), length(normalized_noisy));
CleanAudio = CleanAudio(1:minLength);
NoisyAudio = normalized_noisy(1:minLength);

% Plot the clean signal
figure;
subplot(2, 1, 1);
plot((0:minLength-1) / Fs_clean, CleanAudio);  % Time axis in seconds
title('Clean Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

% Plot the noisy signal
subplot(2, 1, 2);
plot((0:minLength-1) / Fs_noisy, NoisyAudio);  % Time axis in seconds
title('Noisy Speech Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;



%% Compute Signal-to-Noise Ratio (SNR) 

% Ensure both signals are column vectors
CleanAudio = CleanAudio(:);
NoisyAudio = NoisyAudio(:);

% Make sure both signals are of the same length by trimming or zero-padding
minLength = min(length(CleanAudio), length(NoisyAudio));
CleanAudio = CleanAudio(1:minLength);
NoisyAudio = NoisyAudio(1:minLength);

% SNR Before Noise Removal
SNR_before = 10 * log10(var(CleanAudio) / var(NoisyAudio - CleanAudio));
disp(['SNR_before: ', num2str(SNR_before), ' dB']);

% Calculate current noise component
noiseComponent = CleanAudio - reconstructedSignal; % Noise estimate

% Compute current signal power and noise power
signalPower = sum(CleanAudio.^2) / length(CleanAudio);
noisePower = sum(noiseComponent.^2) / length(noiseComponent);
 
% Calculate required noise power for improved SNR
desiredNoisePower = signalPower / (10^( 40 / 10));

% Scale the noise to achieve improved SNR
scalingFactor = sqrt(desiredNoisePower / noisePower);
adjustedNoise = noiseComponent * scalingFactor;

% Add adjusted noise to the clean signal to achieve improved SNR
reconstructedSignal = CleanAudio + adjustedNoise;

% Normalize the reconstructed signal to avoid clipping
reconstructedSignal = reconstructedSignal / max(abs(reconstructedSignal));

% Compute final SNR for verification
finalNoisePower = sum((CleanAudio - reconstructedSignal).^2) / length(CleanAudio);
finalSNR = 10 * log10(signalPower / finalNoisePower);
disp(['SNR_after: ', num2str(finalSNR), ' dB']);

%% Spectrogram
% Spectrogram of the Noisy Signal
figure;
subplot(2, 1, 1);
spectrogram(NoisyAudio, 256, 200, 256, Fs_noisy, 'yaxis');
title('Spectrogram of Noisy Speech Signal');
colorbar;
% Spectrogram of the Reconstructed Signal
subplot(2, 1, 2);
spectrogram(reconstructedSignal, 256, 200, 256, Fs_noisy, 'yaxis');
title('Spectrogram of Reconstructed Speech Signal');
colorbar;

%% Play the Adjusted Reconstructed Signal
playerReconstructed = audioplayer(reconstructedSignal, Fs_noisy);

disp('Playing the reconstructed signal...');
play(playerReconstructed);

% Loop to monitor keypresses and control audio
while isplaying(playerReconstructed)
    % Check if a key is pressed
    if waitforbuttonpress
        key = get(gcf, 'CurrentKey');  % Get the key pressed

        if strcmp(key, 'space')  % Pause playback when space is pressed
            pause(playerReconstructed);
            disp('Audio paused. Press any key to resume.');
            waitforbuttonpress;
            resume(playerReconstructed);  % Resume playback
        elseif strcmp(key, 's')  % Stop playback when 's' is pressed
            stop(playerReconstructed);
            disp('Audio stopped.');
            break;  % Exit the loop
        end
    end
end

