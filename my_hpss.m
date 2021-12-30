function [h, p] = my_hpss(segment, Fs, threshold, N)
    % this function calculates the harmonics and percussion of a segment
    win = sqrt(hann(N,"periodic"));
    overlapLength = floor(numel(win)/2);
    fftLength = 2^nextpow2(numel(win) + 1);
    y = stft(segment, ...
        "Window",win, ...
        "OverlapLength",overlapLength, ...
        "FFTLength",fftLength, ...
        "Centered",true);
    halfIdx = 1:ceil(size(y,1)/2);
    yhalf = y(halfIdx,:);
    ymag = abs(yhalf);
    
    timeFilterLength = 0.2;
    timeFilterLengthInSamples = timeFilterLength/((numel(win) - overlapLength)/Fs);
    ymagharm = movmedian(ymag,timeFilterLengthInSamples,2);
    
    frequencyFilterLength = 500;
    frequencyFilterLengthInSamples = frequencyFilterLength/(Fs/fftLength);
    ymagperc = movmedian(ymag,frequencyFilterLengthInSamples,1);
    
    totalMagnitudePerBin = ymagharm + ymagperc;
    
    harmonicMask = ymagharm > (totalMagnitudePerBin*threshold);
    percussiveMask = ymagperc > (totalMagnitudePerBin*threshold);
    
    yharm = harmonicMask.*yhalf;
    yperc = percussiveMask.*yhalf;
    
    yharm = cat(1,yharm,flipud(conj(yharm)));
    yperc = cat(1,yperc,flipud(conj(yperc)));
    
    h = istft(yharm, ...
        "Window",win, ...
        "OverlapLength",overlapLength, ...
        "FFTLength",fftLength, ...
        "ConjugateSymmetric",true);
    p = istft(yperc, ...
        "Window",win, ...
        "OverlapLength",overlapLength, ...
        "FFTLength",fftLength, ...
        "ConjugateSymmetric",true);


end