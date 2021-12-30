classdef utilities
    methods(Static)
        function [mel_filters] = get_filters(lowest_freq, highest_freq, num_filters, num_fft, sampling_rate, is_normalized)
            % the function calculates the mel filter banks by taking following inputs
            % lowest_freq: lowest frequency value desired to be considered in energy calculations 
            % highest_freq: highest frequencies desired to be considered in energy calculations
            % num_filters: the number of mel filters used in filtering operations
            % num_fft : the number of points of which fast fourier taken
            % sampling_rate : the sampling rate of the audio
            % is_normalized : a boolean to indicate the filters should be normalized according to their bandwidth in frequency domain
            
            center_pts = utilities.calculate_centers(lowest_freq, highest_freq, num_filters, num_fft, sampling_rate); % get center frequencies of triangle filters
            center_pts = center_pts+1; % to get 1 index as starting point
            frame_size = num_fft / 2 + 1; % calculate the window length for fft storage
            mel_filters = zeros(num_filters,frame_size); % initialize the vector that stores all filters
            for i=2:num_filters+1 % for each filter
                mel_filter = zeros(1,frame_size); % initialize the individual filter vectors
                previous_center = center_pts(i-1); % a variable to store previous filter's center point
                center = center_pts(i); % get current filter's center point
                next_center = center_pts(i+1); % get next filter's center point
                bandwidth = next_center-previous_center; % get the bandwidth of the filter
                amplitude = 1; % initialize the amplitude of the filter
                if is_normalized % if filters desired to be normalized 
                    amplitude= 1/bandwidth;
                end
                delta_rise = amplitude/(center - previous_center); % slope of the rising edge
                delta_fall = amplitude/(center - next_center); % slope of the falling edge
                mel_filter(previous_center:center)=(0:delta_rise:amplitude); % assign rising edge
                mel_filter(center:next_center) = (amplitude:delta_fall:0); % assign falling edge
                mel_filters(i-1,:) = mel_filter; % assign finalized filter
            end
            % following lines are coded to plot the generated mel filter bank
%             widths = center_pts(2:end) - center_pts(1:end-1); 
%             figure
%             hold on
%             freq = (0:1:num_fft/2);
%             for i =1:num_filters
%                 plot(freq, mel_filters(i,:))
%             end
%             title("Mel Weights")
%             xlabel('freq')
%             ylabel('weight')
%             grid on
%             xlim([0 num_fft/2])
        end

        function [mel_centers] = calculate_centers(lowest_freq, highest_freq, num_filters, num_fft, sampling_rate)
            % a function to calculate center points of the mel filters by the following inputs
            % lowest_freq: lowest frequency value desired to be considered in energy calculations 
            % highest_freq: highest frequencies desired to be considered in energy calculations
            % num_filters: the number of mel filters used in filtering operations
            % num_fft : the number of points of which fast fourier taken
            % sampling_rate : the sampling rate of the audio
            mel_low = utilities.freq2mel(lowest_freq); % get the lowest point in the mel domain
            mel_high = utilities.freq2mel(highest_freq); % get the highest point in the mel domain
            mel_centers = (mel_low:(mel_high-mel_low)/(num_filters+1):mel_high); % divide mel domain into equally spaced parts
            mel_centers = utilities.mel2freq(mel_centers); % convert division points into the frequency domain
            mel_centers = floor((num_fft+1).*mel_centers./sampling_rate); % convert center values into integers to be used in fft indexes
        end
        % the conversion is done by the formulae here https://en.wikipedia.org/wiki/Mel_scale
        
        function [freq_data]=mel2freq(mel_data_points)
            freq_data=700.*(10.^(mel_data_points./2595)-1);
        end

        function [mel_data] = freq2mel(freq_data_points)
            mel_data = 2595.*log10(1+freq_data_points./700);
        end

        function [coeff,filtered]=mfcc_calc(data, sample_rate, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, is_normalized)
            % the function calculates the Mel Frequency Cepstral Coefficients by getting following inputs
            % data: one channel audio signal 
            % sampling_rate : the sampling rate of the audio
            % frame_size : the window size for each time frame 
            % hop_size : the shifting size of the time windows
            % num_of_filters: the number of mel filters used in filtering operations
            % FFT_size : the number of points of which fast fourier taken
            % dct_coefficients : the number of coefficients that will be calculated
            % is_normalized : a boolean to indicate the filters should be normalized according to their bandwidth in frequency domain
            
            frame_len = floor(sample_rate * frame_size / 1000); % calculate the length of each sample frame
            frame_shift = sample_rate * hop_size / 1000; % calculate the length of each shifting
            frame_num = floor((length(data)-frame_len) / frame_shift)+1; % calculate the resulting number of frames
            frames = zeros(FFT_size/2+1,frame_num); % initialize the frames
            hamming_win = hamming(frame_len); % get the hamming filter
            for i = 1:frame_num % for each frame
                windowed = data((i-1)*frame_shift+1:(i-1)*frame_shift+frame_len).*hamming_win; % get hamming filtered windowed data
                fft_result = fft(windowed,FFT_size); % calculate the fft of the current window
                frames(:,i) = fft_result(1:end/2+1,1); % assign the fft into previously initialized frames
            end
            
            periodogram_estimate=(abs(frames).^2)./((FFT_size/2)+1); % get periodogram estimate by calculating the energy and normalizing the vector with the length
            filters = utilities.get_filters(1, sample_rate/2, num_of_filters, FFT_size, sample_rate, is_normalized); % get the mel filter bank
            filtered = filters*periodogram_estimate; % filter the FFT data matrix
            filtered = 10.*log10(filtered); % get dB values
            dct_matrix=utilities.get_dct(dct_coefficients, num_of_filters); % get the dct matrix to turn data from cepstrum to coefficients
            coeff = dct_matrix*filtered; % calculate the coefficients for each time frame
        end
        
        function [dct_matrix]=get_dct(num_dct_coeff, num_mel_filters)
            % the function calculates the dct matrix which turns data from cesptrum to coefficients if it is multiplied with mel cepstrum
            % num_dct_coeff: the number of desired mfc coefficients
            % sampling_rate : the sampling rate of the audio
            % frame_size : the window size for each time frame 
            % The DCT-II unitary algorithm is defined in "Discrete Time Signal
            % Processing", equations (8.158), (8.159), (8.160).
            dct_matrix  = zeros(num_dct_coeff,num_mel_filters);
            A       = sqrt(1/num_mel_filters);
            B       = sqrt(2/num_mel_filters);
            C       = 2*num_mel_filters;
            fixed_coeff = 2*pi/C;

            dct_matrix(1,:) = A;
            for k = 1:num_mel_filters
                for n = 2:num_dct_coeff
                    dct_matrix(n,k) = B*cos(fixed_coeff*(n-1)*(k-0.5));
                end
            end
        end
    end
end