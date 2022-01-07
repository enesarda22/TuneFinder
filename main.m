%% reading audio

[segment, Fs] = audioread("data/genre_data/blues/blues.00000.wav");
segment = segment(1:4*Fs);

spectrogram(segment,1024,512,1024,Fs,"yaxis");
title("Auido");
sound(segment, Fs);

%% Harmonic Percussive Source Separation
th1 = 0.7;
N1 = 4096;
[h1, p1] = my_hpss(segment, Fs, th1, N1);

th2 = 0.6;
N2 = 256;
[h2, p2] = my_hpss(segment, Fs, th2, N2);

%% harmonic
sound(h1, Fs);

spectrogram(h1,1024,512,1024,Fs,"yaxis");
title("Harmonic Part");

%% percussive
sound(p2, Fs);

spectrogram(p2,1024,512,1024,Fs,"yaxis")
title("Percussive Part")

%% spectral features

hop_length = 0.010; % hop length is set in s
window_length = 0.025; % window length is set in s
n_fft = 2048; % fft length is set
    
frame_len = floor(window_length*Fs); % frame length in samples
frame_shift = floor(hop_length*Fs); % frame shift in samples
    
n_frames = floor((length(segment)-frame_len) / frame_shift)+1; % number of frames
    
% frames is initialized
frames = zeros(frame_len, n_frames); 
for i = 1:n_frames
    frame = segment((i-1)*frame_shift+1:(i-1)*frame_shift+frame_len);
    frames(:, i) = frame;
end
    
frames(:,all(frames == 0))=[]; % removes zero columns
n_frames = size(frames, 2); % number of frames updated

repeated_win = repmat(hamming(frame_len), [1, n_frames]); % window is repeated for each column
windowed = frames .* repeated_win;

W = fft(windowed, n_fft); % fft of windowed frames is taken

half_idx = n_fft/2; % last non-negative index
nngt_indices = (1:half_idx)'; % non-negative indices are defined
repeated_indices = repmat(nngt_indices, [1, n_frames]);

normalized_W = abs(W) ./ sum(abs(W(nngt_indices, :))); % W is normalized
C = sum(repeated_indices .* normalized_W(nngt_indices, :)); % spectral centroid is defined

% moments are calculated
second_moment = sum((repeated_indices-C).^2 .* normalized_W(nngt_indices, :)); % spread
third_moment = sum((repeated_indices-C).^3 .* normalized_W(nngt_indices, :)) ./ (second_moment.^(3/2)); % skewness
fourth_moment = sum((repeated_indices-C).^4 .* normalized_W(nngt_indices, :)) ./ (second_moment.^2); % kurtosis

% plot two frames
subplot(211);
plot(normalized_W(1:half_idx, 1));
xlim([1, half_idx]);
ylabel("Amplitude");
xlabel("Non-negative Indices");
title("center=" + num2str(C(1)) + " spread=" + num2str(second_moment(1) + " skewness=" + num2str(third_moment(1) + " kurtosis=" + num2str(fourth_moment(1)))));
grid on;
subplot(212);
plot(normalized_W(1:half_idx, 7));
xlim([1, half_idx]);
ylabel("Amplitude");
xlabel("Non-negative Indices");
title("center=" + num2str(C(7)) + " spread=" + num2str(second_moment(7) + " skewness=" + num2str(third_moment(7) + " kurtosis=" + num2str(fourth_moment(7)))));
grid on;

%% feature extraction
classes = ["blues", "classical", "country", "disco", "hiphop", "jazz", "metal", "pop", "reggae", "rock"]; % classes are defined
Fs = 22050; % sampling rate
dur = 3; % duration of segments in sec
n_samples = Fs*dur; % number of samples
n_segments = floor(30/dur); % number of segments each song

rootdir = 'data';
filelist = dir(fullfile(rootdir, '**/*.wav'));
data = zeros(10000, 52); % data is initialized
y = string([10000, 1]); % target labels are initialized

idx = 1; % sample index is initialized

for i = 1:length(filelist) % iterates through files

    % genre name is extracted as the label
    file = filelist(i);
    str = strsplit(file.name, '.');
    label = str{1};

    % tries to read the audio file
    try
        track = audioread(strcat(file.folder, '/', file.name));
    catch
        fprintf("%d is broken\n", i);
    end

    for j = 1:n_segments-1 % iterates through all segments
        segment = track((j-1)*n_samples+1:j*n_samples); % segment is formed
        data(idx, :) = get_features(segment, Fs); % features are extracted

        y(idx) = label; % label is formed
        idx = idx+1; % index is incremented
    end

    segment = track((n_segments-1)*n_samples+1:end); % last segment
    data(idx, :) = get_features(segment, Fs); % last data
    
    y(idx) = label; % label is formed
    idx = idx+1; % index is incremented
    fprintf("%s is added\n", file.name);
end

% matrices are stored
writematrix(data, 'data.txt');
writematrix(y, 'y.txt');

%% shazam feature extraction
classes = ["all_i_need", "fade_to_black", "sarkilarin_gozu", "yalan"];

for i = 1:length(classes)
    class = classes(i);

    [track, Fs] = audioread("data/shazam_data/" + class + ".mp3");
    track = track(:, 1);
    track = track(30*Fs:90*Fs);
    
    window_length = 5; % window in sec
    hop_length = 3; % hop in sec

    frame_len = floor(window_length*Fs);
    frame_shift = floor(hop_length*Fs);

    n_frames = floor((length(track)-frame_len) / frame_shift)+1;
    
    data = zeros(n_frames, 52);
    
    for j = 1:n_frames
        segment = track((j-1)*frame_shift+1:(j-1)*frame_shift+frame_len);
        data(j, :) = get_features(segment, Fs);
    end
    
%     segment = track((n_segments-1)*n_samples+1:end);
%     data(n_segments, :) = get_features(segment, Fs);
    
    writematrix(data, 'extracted_features/shazam/recording/5sec_new' + class + '.txt');
end

%% read shazam features
rootdir = 'extracted_features/shazam/recording';
filelist = dir(fullfile(rootdir, '*.txt'));

file = filelist(1).name;
A = readmatrix(strcat(rootdir, '/', file));
labels = ones(size(A, 1), 1);
fprintf("%s is added \n", file);

for i = 2:length(filelist)
    file = filelist(i).name;
    data = readmatrix(strcat(rootdir, '/', file));
    A = [A; data];
    labels = [labels; ones(size(data, 1), 1)*i];

    fprintf("%s is added \n", file);
end

%% k-nearest neigbor
Fs = 44100;
classes = ["all i need", "fade to black", "sarkilarin gozu", "yalan"];

recObj = audiorecorder(Fs, 16, 1);
disp('Listening....');

recordblocking(recObj, 5);
track = getaudiodata(recObj);
track = track(1500:end);

test = get_features(track, Fs);

k = 9;

distances = vecnorm(A' - test'); % finds the distance to training samples
[~, idx] = sort(distances, 'ascend'); % sorts the distances in ascending order

% sims = cosineSimilarity(normalized_A, test);
% [~, idx] = sort(sims, 'ascend');

y = mode(labels(idx(1:k))); % assigned class is the one that appears the most among the k nearest neighbor
sum(labels(idx(1:k)) == y)/k
fprintf("%s\n", classes(y));


%% fisher linear discriminant
A = readmatrix("extracted_features/after_spec/data.txt");
classes = ["blues", "classical", "country", "disco", "hiphop", "jazz", "metal", "pop", "reggae", "rock"]; % classes are defined
d_prime = 52; % # of features
n = 1000; % sample size

% M is the mean matrix with each column corresponding to mean of that class
M = zeros(d_prime, size(classes, 2));
% B is created to find the S_W
B = zeros(d_prime, size(A, 1));
for i = 1:size(classes, 2)

    reduced_A = A(n*(i-1)+1:n*i, :);
    M(:, i) = mean(reduced_A)';
    B(:, n*(i-1)+1:n*i) = reduced_A' - M(:, i);

end

% S_W is calculated to hold the seperation within each class from the mean
% of the class
S_W = B*B';

% S_B is calculated to hold the seperation between the mean of each class
% from the mean of all the dataset
S_B = zeros(d_prime, d_prime);
m_all = mean(A)';
for col = M
    S_B = S_B + (col-m_all)*(col-m_all)';
end
S_B = S_B * n;

% W matrix is a (d_prime * c-1) matrix whose columns are the c-1
% eigenvectors of S_W\S_B corresponding to c-1 greatest eigenvalues
[W, lambdas] = eig(S_W\S_B, 'vector');
[~, ind] = sort(abs(lambdas), 'descend');
W = W(:, 1:size(classes, 2)-1);

% J is calculated using frobenius norm
J = norm(W'*S_B*W, "fro") / norm(W'*S_W*W, "fro");

fprintf("J for the 10 class FLD = %.2f\n", J);



