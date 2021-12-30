function [rms_coeffs] = rms_features(segment, Fs)
    
    hop_length = 0.010; % 10msec
    window_length = 0.025; % 25msec
    
    frame_len = floor(window_length*Fs); % calculate the frame width
    frame_shift = floor(hop_length*Fs); % calculate the shift size of frames
    
    n_frames = floor((length(segment)-frame_len) / frame_shift)+1; % calculate the number of frames
    RMS = zeros(n_frames, 1); % initialize RMS matrix
    
    for i = 1:n_frames % each frame
        windowed = segment((i-1)*frame_shift+1:(i-1)*frame_shift+frame_len) .* hamming(frame_len); % apply hamming filter
        RMS(i) = rms(windowed); % get rms value
    end
    
    RMS = RMS/rms(segment); % normalize the matrix with the rms of the overall data
    RMS_mean = mean(RMS);
    RMS_var = var(RMS);

    rms_coeffs = [RMS_mean, RMS_var];

end