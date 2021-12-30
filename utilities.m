classdef utilities
    methods(Static)
        function [mel_filters] = get_filters(lowest_freq, highest_freq, num_filters, num_fft, sampling_rate, is_normalized)
            center_pts = utilities.calculate_centers(lowest_freq, highest_freq, num_filters, num_fft, sampling_rate);
            center_pts = center_pts+1; % to get 1 index as starting point
            frame_size = num_fft/2+1;
            mel_filters = zeros(num_filters,frame_size);
            for i=2:num_filters+1
                mel_filter = zeros(1,frame_size);
                previous_center = center_pts(i-1);
                center = center_pts(i);
                next_center = center_pts(i+1);
                bandwidth = next_center-previous_center;
                amplitude = 1;
                if is_normalized
                    amplitude= 1/bandwidth;
                end
                delta_rise = amplitude/(center - previous_center);
                delta_fall = amplitude/(center - next_center);
                mel_filter(previous_center:center)=(0:delta_rise:amplitude);
                mel_filter(center:next_center) = (amplitude:delta_fall:0);
                mel_filters(i-1,:) = mel_filter;
            end
            % For normalization divide amplitude of the triangles with their widths
            % widths = center_pts(2:end) - center_pts(1:end-1); 
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
            mel_low = utilities.freq2mel(lowest_freq);
            mel_high = utilities.freq2mel(highest_freq);
            mel_centers = (mel_low:(mel_high-mel_low)/(num_filters+1):mel_high);
            mel_centers = utilities.mel2freq(mel_centers);
            mel_centers = floor((num_fft+1).*mel_centers./sampling_rate);
        end
        function [freq_data]=mel2freq(mel_data_points)
            freq_data=700.*(10.^(mel_data_points./2595)-1);
        end

        function [mel_data] = freq2mel(freq_data_points)
            mel_data = 2595.*log10(1+freq_data_points./700);
        end

        function [coeff,filtered]=mfcc_calc(data, sample_rate, frame_size, hop_size, num_of_filters, FFT_size, dct_coefficients, is_normalized)
            frame_len = floor(sample_rate * frame_size / 1000); % calculate the length of each sample frame
            frame_shift = sample_rate * hop_size / 1000; % calculate the length of each shifting
            frame_num = floor((length(data)-frame_len) / frame_shift)+1;
            frames = zeros(FFT_size/2+1,frame_num);
            hamming_win = hamming(frame_len);
            for i = 1:frame_num
                windowed = data((i-1)*frame_shift+1:(i-1)*frame_shift+frame_len).*hamming_win;
                fft_result = fft(windowed,FFT_size);
                frames(:,i) = fft_result(1:end/2+1,1);
            end
            periodogram_estimate=(abs(frames).^2)./((FFT_size/2)+1);
            filters = utilities.get_filters(1, sample_rate/2, num_of_filters, FFT_size, sample_rate, is_normalized);
            filtered = filters*periodogram_estimate;
            filtered = 10.*log10(filtered);
            dct_matrix=utilities.get_dct(dct_coefficients, num_of_filters);
            coeff = dct_matrix*filtered;
        end
        
        function [dct_matrix]=get_dct(num_dct_coeff, num_mel_filters)
            dct_matrix  = zeros(num_dct_coeff,num_mel_filters);
            A       = sqrt(1/num_mel_filters);
            B       = sqrt(2/num_mel_filters);
            C       = 2*num_mel_filters;
            fixed_coeff = 2*pi/C;

            coder.gpu.kernel;
            dct_matrix(1,:) = A;
            for k = 1:num_mel_filters
                for n = 2:num_dct_coeff
                    dct_matrix(n,k) = B*cos(fixed_coeff*(n-1)*(k-0.5));
                end
            end
        end
    end
end