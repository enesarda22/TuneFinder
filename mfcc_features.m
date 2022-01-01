function [mfcc_coeffs] = mfcc_features(segment, Fs)
    % a function for getting Mel Frequency Coefficients with the following
    % fixed hyperparameters
   
    FFT_size = 2048; % which defines the length of the FFT
    hop_size = 10; % hop size in ms
    frame_size = 25; % the frame size in ms
    num_of_filters = 64; % number of filters 
    dct_coefficients = 20; % number of desired mfc coefficients

    [coeffs, ~]=utilities.mfcc_calc(segment, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true, false); % this function gets the mfccs
    coeffs = coeffs';
    coeffs = coeffs(all(~isnan(coeffs),2),:); % get only number vals

    mfcc_mean = mean(coeffs); % get mean
    mfcc_var = var(coeffs); % get variance

    mfcc_coeffs = [mfcc_mean, mfcc_var];

end