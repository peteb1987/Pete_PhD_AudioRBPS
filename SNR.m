function [ snr ] = SNR( true_audio, noisy_audio )
%SNR Calculate signal to noise ratio

snr = 10*log10( sum(true_audio.^2)/sum((noisy_audio - true_audio).^2) );

end

