function [hpss_coeffs] = hpss_features(track, Fs)

th1 = 0.7;
N1 = 4096;
[h, ~] = my_hpss(track, Fs, th1, N1);

th2 = 0.6;
N2 = 256;
[~, p] = my_hpss(track, Fs, th2, N2);

hpss_coeffs = zeros(1, 4);
hpss_coeffs(1:2) = [mean(h), var(h)];
hpss_coeffs(3:4) = [mean(p), var(p)];

end