function [rms_coeffs] = rms_features(segment, Fs)
    
    hop_length = 0.010;
    window_length = 0.025;
    
    frame_len = floor(window_length*Fs);
    frame_shift = floor(hop_length*Fs);
    
    n_frames = floor((length(segment)-frame_len) / frame_shift)+1;
    RMS = zeros(n_frames, 1);
    
    for i = 1:n_frames
        windowed = segment((i-1)*frame_shift+1:(i-1)*frame_shift+frame_len) .* hamming(frame_len);
        RMS(i) = rms(windowed);
    end
    
    RMS = RMS/rms(segment);
    RMS_mean = mean(RMS);
    RMS_var = var(RMS);

    rms_coeffs = [RMS_mean, RMS_var];

end