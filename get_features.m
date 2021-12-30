function [features] = get_features(segment, Fs)

    features = zeros(1, 46);

    features(1, 1:4) = hpss_features(segment, Fs); % hpss features are added
    features(1, 5:6) = rms_features(segment, Fs); % rms features are added
    features(1, 7:46) = mfcc_features(segment, Fs); % mfcc features are added

end