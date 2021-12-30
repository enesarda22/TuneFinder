function [spectral_coeffs] = spectral_features(segment, Fs)

    hop_length = 0.010;
    window_length = 0.025;
    n_fft = 2048;
    
    frame_len = floor(window_length*Fs);
    frame_shift = floor(hop_length*Fs);
    
    n_frames = floor((length(segment)-frame_len) / frame_shift)+1;
    
    frames = zeros(frame_len, n_frames);
    for i = 1:n_frames
        frame = segment((i-1)*frame_shift+1:(i-1)*frame_shift+frame_len);
        frames(:, i) = frame;
    end
    
    frames(:,all(frames == 0))=[]; % removes zero columns
    n_frames = size(frames, 2);

    repeated_win = repmat(hamming(frame_len), [1, n_frames]);
    windowed = frames .* repeated_win;

    W = fft(windowed, n_fft);

    half_idx = n_fft/2;
    nngt_indices = (1:half_idx)'; % non-negative indices are defined
    repeated_indices = repmat(nngt_indices, [1, n_frames]);

    normalized_W = abs(W) ./ sum(abs(W(nngt_indices, :)));
    C = sum(repeated_indices .* normalized_W(nngt_indices, :)); % spectral centroid is defined

    second_moment = sum((repeated_indices-C).^2 .* normalized_W(nngt_indices, :));
    third_moment = sum((repeated_indices-C).^3 .* normalized_W(nngt_indices, :)) ./ (second_moment.^(3/2));
    fourth_moment = sum((repeated_indices-C).^4 .* normalized_W(nngt_indices, :)) ./ (second_moment.^2);

    spectral_coeffs = [mean(second_moment), var(second_moment), mean(third_moment), var(third_moment), mean(fourth_moment), var(fourth_moment)];

end