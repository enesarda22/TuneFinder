function [mfcc_coeffs] = mfcc_features(segment, Fs)

    FFT_size = 2048; % which defines the length of the FFT
    hop_size = 10; % hop size in ms
    frame_size = 25; % the frame size in ms
    num_of_filters = 64;
    dct_coefficients = 20;

    [coeffs, ~]=utilities.mfcc_calc(segment, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
    coeffs = coeffs';
    coeffs = coeffs(all(~isnan(coeffs),2),:);

    mfcc_mean = mean(coeffs);
    mfcc_var = var(coeffs);

    mfcc_coeffs = [mfcc_mean, mfcc_var];

end