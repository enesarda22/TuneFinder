% main
%%
close all
clearvars *
clc
%%
[segment, Fs] = audioread("Data/genres_original/classical/classical.00055.wav");
%data = data / max(abs(data)); % normalize
%% hyperparameters
Fs = 22050; % sampling rate

FFT_size = 2048; % which defines the length of the FFT
sample_rate = Fs;
hop_size = 10; % hop size in ms
frame_size = 25; % the frame size in ms
num_of_filters = 64;
dct_coefficients = 20;
%% windowing the input sound
[segment, Fs] = audioread("Data/genres_original/blues/blues.00000.wav");
segment = segment(1:3*Fs);

% [coeffs, filtered]=utilities.mfcc_calc(segment, sample_rate, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
% coeffs = coeffs';
% coeffs = coeffs(all(~isnan(coeffs),2),:);
% 
% features = [mean(coeffs), var(coeffs)];
sound(segment, Fs);

%% HPSS
th1 = 0.7;
N1 = 4096;
[h1, p1] = my_hpss(segment, Fs, th1, N1);

th2 = 0.6;
N2 = 256;
[h2, p2] = my_hpss(segment, Fs, th2, N2);

%% harmonic
sound(h1, Fs);

spectrogram(h1, 1024,512,1024,Fs,"yaxis")
title("Recovered Harmonic Audio")

%% percussive
sound(p2, Fs);

spectrogram(p2,1024,512,1024,Fs,"yaxis")
title("Recovered Percussive Audio")

%% feature extraction
classes = ["blues", "classical", "country", "disco", "hiphop", "jazz", "metal", "pop", "reggae", "rock"]; % classes are defined
Fs = 22050; % sampling rate
dur = 30; % duration of segments in sec
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
classes = ["all_i_need"];

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
    
    writematrix(data, 'extracted_features/shazam/recording/5sec_' + class + '.txt');
end

%% fischer
d = 52;
yalan_coeffs = readmatrix("extracted_features/shazam/yalan.txt");
impromtu_coeffs = readmatrix("extracted_features/shazam/impromptu.txt");

% S_W is initialized as zero matrix
S_W = zeros(d, d);

% means of the partitioned matrices are found and S_W matrix is formed from
% the definition
m1 = mean(impromtu_coeffs)';
for x = impromtu_coeffs'
    S_W = S_W + (x-m1)*(x-m1)';
end

m2 = mean(yalan_coeffs)';
for x = yalan_coeffs'
    S_W = S_W + (x-m2)*(x-m2)';
end

% S_B matrix and w are found from the definitions
S_B = (m1-m2) * (m1-m2)';
w = S_W \ (m1-m2);

% J that should be maximised is found
J = (w'*S_B*w) / (w'*S_W*w)

% projections of the samples on one dimension are calculated
y_impromptu = w'*impromtu_coeffs';
y_yalan = w'*yalan_coeffs';

% projections are plotted
x_impromptu = zeros(1, size(y_impromptu, 2));
x_yalan = zeros(1, size(y_yalan, 2));
scatter(y_impromptu, x_impromptu, 'filled', 'g');
hold on;
scatter(y_yalan, x_yalan, 'filled', 'b');
grid on;
legend('impromptu','yalan');
title(["training data J =", num2str(J)]);

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

%% test
recObj = audiorecorder(Fs, 16, 1);
disp('Listening....');

recordblocking(recObj, 3);
track = getaudiodata(recObj);

test = get_features(track, Fs);

y_old = test*w;
is_found = false;

while ~is_found

    recordblocking(recObj, 3);
    track = getaudiodata(recObj);
    
    test = get_features(track, Fs);
    
    y_new = test*w;
    if y_new<8
        if y_old<8
            is_found = true;
            c = "yalan\n";
        else
            y_old = y_new;
        end
    else
        if y_old>=8
            is_found = true;
            c = "impromptu\n";
        else
            y_old = y_new;
        end
    end

end

fprintf(c);


%% k-nearest neigbor

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

%%
classes = ["fade to black", "impromptu"];

recObj = audiorecorder(Fs, 16, 1);
disp('Listening....');

recordblocking(recObj, 5);
track = getaudiodata(recObj);
track = track ./ max(track);

window_length = 3; % window in sec
hop_length = 1; % hop in sec

frame_len = floor(window_length*Fs);
frame_shift = floor(hop_length*Fs);

n_frames = floor((length(track)-frame_len) / frame_shift)+1;

for j = 1:n_frames
    segment = track((j-1)*frame_shift+1:(j-1)*frame_shift+frame_len);
    test = get_features(segment, Fs);

    k = 5;
    
    distances = vecnorm(A' - test'); % finds the distance to training samples
    [~, idx] = sort(distances, 'ascend'); % sorts the distances in ascending order
    
    % sims = cosineSimilarity(normalized_A, test);
    % [~, idx] = sort(sims, 'ascend');
    
    y = mode(labels(idx(1:k))); % assigned class is the one that appears the most among the k nearest neighbor
    sum(labels(idx(1:k)) == y);
    fprintf("%s\n", classes(y));

end



%% random forest

% model = fitctree(A, labels);
% t = templateTree('NumVariablesToSample','all',...
%     'PredictorSelection','interaction-curvature','Surrogate','on');
% model = fitrensemble(A, labels, 'Method','Bag','NumLearningCycles',125, ...
%     'Learners',t);

% classes = ["bir kadin cizeceksin", "bu son olsun", "ihtimal", "impromptu", "yalan"];
classes = ["bit ihtimal", "bir kadin", "bu son"];
% model = fitcknn(A, labels);

%%

recObj = audiorecorder(Fs, 16, 1);
disp('Listening....');

recordblocking(recObj, 3);
track = getaudiodata(recObj);
test = get_features(track, Fs);
test*W


% y = model.predict(test);
% fprintf("%s\n", classes(round(y)));

%% fisher
classes = ["fade to black", "in da club", "sarkilarin gozu", "yalan"];
d_prime = 52;
n = 19;

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

fprintf("J for the 3 class FLD = %.2f\n", J);

