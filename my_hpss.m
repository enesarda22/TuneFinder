% this function calculates the harmonics and percussion of a segment
function [h, p] = my_hpss(segment, Fs, threshold, N)
    
    win = sqrt(hann(N,"periodic")); % window length
    overlap_length = floor(numel(win)/2); % overlap length
    n_fft = 2^nextpow2(numel(win) + 1); % fft length

    % spectrogram is found
    y = stft(segment, ...
        "Window",win, ...
        "OverlapLength",overlap_length, ...
        "FFTLength",n_fft, ...
        "Centered",true);

    half_idx = 1:ceil(size(y,1)/2); % last positive index
    y_half = y(half_idx,:); % positive frequency components
    y_mag = abs(y_half); % magnitudes of positive frequency components
    
    % median smoothing along the time axis to enhance harmonic part
    time_filter_length = 0.2; 
    n_time_filter = time_filter_length/((numel(win) - overlap_length)/Fs); % time filter length in samples
    ymag_harm = movmedian(y_mag,n_time_filter,2);
    
    % median smoothing along the frequency axis to enhance percussive part
    freq_filter_length = 500;
    n_freq_filter = freq_filter_length/(Fs/n_fft);
    ymag_perc = movmedian(y_mag,n_freq_filter,1);
    
    tot_mag_per_bin = ymag_harm + ymag_perc; % total magnitude per bin
    
    harmonic_mask = ymag_harm > (tot_mag_per_bin*threshold); % binary harmonic mask is formed
    percussive_mask = ymag_perc > (tot_mag_per_bin*threshold); % binary percussive mask is formed
     
    y_harm = harmonic_mask.*y_half; % harmonic part is found
    y_perc = percussive_mask.*y_half; % percussive part is found
    
    % half sided spectrum is mirrored to get two sided conj symmetric
    % spectrum
    y_harm = cat(1,y_harm,flipud(conj(y_harm)));
    y_perc = cat(1,y_perc,flipud(conj(y_perc)));
    
    % inverse stft is formed to return the signals to time domain
    h = istft(y_harm, ...
        "Window",win, ...
        "OverlapLength",overlap_length, ...
        "FFTLength",n_fft, ...
        "ConjugateSymmetric",true);
    p = istft(y_perc, ...
        "Window",win, ...
        "OverlapLength",overlap_length, ...
        "FFTLength",n_fft, ...
        "ConjugateSymmetric",true);


end