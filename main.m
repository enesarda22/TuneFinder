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
segment = segment(1:4*Fs);

[coeffs, filtered]=utilities.mfcc_calc(segment, sample_rate, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
coeffs = coeffs';
coeffs = coeffs(all(~isnan(coeffs),2),:);

features = [mean(coeffs), var(coeffs)];
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

% %% residual
% sound(r, Fs)
% 
% spectrogram(r,1024,512,1024,Fs,"yaxis")
% title("Recovered Residual Audio")


%% feature extraction
classes = ["blues", "classical", "country", "disco", "hiphop", "jazz", "metal", "pop", "reggae", "rock"];
Fs = 22050; % sampling rate
dur = 3; % duration in sec
n_samples = Fs*dur; % number of samples

rootdir = 'Data';
filelist = dir(fullfile(rootdir, '**/*.wav'));
data = zeros(10000, 46);
y = string([10000, 1]);

idx = 1;

for i = 1:length(filelist)

    file = filelist(i);
    str = strsplit(file.name, '.');
    label = str{1};

    try
        track = audioread(strcat(file.folder, '/', file.name));
    catch
        fprintf("%d is broken\n", i);
    end

    for j = 1:9
        segment = track((j-1)*n_samples+1:j*n_samples);
        data(idx, :) = get_features(segment, Fs);

        y(idx) = label;
        idx = idx+1;
    end

    segment = track(9*n_samples+1:end);
    data(idx, :) = get_features(segment, Fs);
    
    y(idx) = label;
    idx = idx+1;
    fprintf("%s is added\n", file.name);
end

writematrix(data, 'data.txt');
writematrix(y, 'y.txt');

%% shazam

[track, Fs] = audioread("yenidena.mp3");
track = track(:, 1);

dur = 3; % duration in sec
n_samples = Fs*dur; % number of samples
n_segments = floor(length(track)/n_samples);

data = zeros(n_segments, 44);

for j = 1:n_segments-1
    segment = track((j-1)*n_samples+1:j*n_samples);
    [coeffs, ~]=utilities.mfcc_calc(segment, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
    coeffs = coeffs';
    coeffs = coeffs(all(~isnan(coeffs),2),:);

    data(j, 1:4) = hpss_features(segment, Fs); % hpss featuers are added
    data(j, 5:24) = mean(coeffs); % mean of the mfcc are added
    data(j, 25:44) = var(coeffs); % variance of the mfcc are added

end

segment = track((n_segments-1)*n_samples+1:end);
[coeffs, ~]=utilities.mfcc_calc(segment, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
coeffs = coeffs';
coeffs = coeffs(all(~isnan(coeffs),2),:);

data(n_segments, 1:4) = hpss_features(segment, Fs); % hpss featuers are added
data(n_segments, 5:24) = mean(coeffs); % mean of the mfcc are added
data(n_segments, 25:44) = var(coeffs); % variance of the mfcc are added

writematrix(data, 'yenidena.txt');

%% fischer
d = 44;
yalan_coeffs = readmatrix("yalan.txt");
yenidena_coeffs = readmatrix("yenidena.txt");

% S_W is initialized as zero matrix
S_W = zeros(d, d);

% means of the partitioned matrices are found and S_W matrix is formed from
% the definition
m1 = mean(yenidena_coeffs)';
for x = yenidena_coeffs'
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
y_yenidena = w'*yenidena_coeffs';
y_yalan = w'*yalan_coeffs';

% projections are plotted
x_yenidena = zeros(1, size(y_yenidena, 2));
x_yalan = zeros(1, size(y_yalan, 2));
scatter(y_yenidena, x_yenidena, 'filled', 'g');
hold on;
scatter(y_yalan, x_yalan, 'filled', 'b');
grid on;
legend('yenidena','yalan');
title(["training data J =", num2str(J)]);

%% test
recObj = audiorecorder(Fs, 16, 1);
disp('Listening....');

recordblocking(recObj, 3);
track = getaudiodata(recObj);

[coeffs, ~]=utilities.mfcc_calc(track, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
coeffs = coeffs';
coeffs = coeffs(all(~isnan(coeffs),2),:);

test = zeros(1, 44);
test(1:4) = hpss_features(track, Fs); % hpss featuers are added
test(5:24) = mean(coeffs); % mean of the mfcc are added
test(25:44) = var(coeffs); % variance of the mfcc are added

y_old = test*w;
is_found = false;

while ~is_found

    recordblocking(recObj, 3);
    track = getaudiodata(recObj);
    
    [coeffs, ~]=utilities.mfcc_calc(track, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
    coeffs = coeffs';
    coeffs = coeffs(all(~isnan(coeffs),2),:);
    
    test = zeros(1, 44);
    test(1:4) = hpss_features(track, Fs); % hpss featuers are added
    test(5:24) = mean(coeffs); % mean of the mfcc are added
    test(25:44) = var(coeffs); % variance of the mfcc are added
    
    y_new = test*w;
    if y_new<0.9
        if y_old<0.9
            is_found = true;
            c = "yalan\n";
        else
            y_old = y_new;
        end
    else
        if y_old>=0.9
            is_found = true;
            c = "yenidena\n";
        else
            y_old = y_new;
        end
    end

end

fprintf(c);

%% k-nearest neigbor

recObj = audiorecorder(Fs, 16, 1);
disp('Listening....');

recordblocking(recObj, 3);
track = getaudiodata(recObj);

[coeffs, ~]=utilities.mfcc_calc(track, Fs, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, true);
coeffs = coeffs';
coeffs = coeffs(all(~isnan(coeffs),2),:);

test = zeros(1, 44);
test(1:4) = hpss_features(track, Fs); % hpss featuers are added
test(5:24) = mean(coeffs); % mean of the mfcc are added
test(25:44) = var(coeffs); % variance of the mfcc are added

labels = repelem([1, 2], [size(yenidena_coeffs, 1), size(yalan_coeffs, 1)]); % 1 is yenidena 2 is yalan
k = 21;

distances = vecnorm(A' - test'); % finds the distance to training samples
[~, idx] = sort(distances, 'ascend'); % sorts the distances in ascending order
% similarities = cosineSimilarity(A, test);
% [~, idx] = sort(similarities, 'descend'); % sorts the similarities in descending order
y = mode(labels(idx(1:k))); % assigned class is the one that appears the most among the k nearest neighbor

if y == 1
    fprintf("yenidena\n");
else
    fprintf("yalan\n");
end

