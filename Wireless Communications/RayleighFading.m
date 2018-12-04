% Flat fading channel
% Generate Rayleigh fading with Gans' Doppler spectrum 
clear
fs = 1000; %Hz, sampling rate
fd = 1; %Hz
Nint = 2^20; % 2^14 = 16k points of the genertaed sequence

% Doppler Spectrum
sqrtpsd=1./(1-([-fd/fs*Nint+.5:1:fd/fs*Nint-.5]/(fd/fs*Nint)).^2).^.25;

% Generate complex Gaussian
ampl=randn(1,floor(2*Nint*fd/fs))+j*randn(1,floor(2*Nint*fd/fs));

h =  ifft(ampl.*sqrtpsd,2*Nint); % row vector    
% Normalization of the channel impulse response
h = sqrt(2*Nint)*h./sqrt(sum(sqrtpsd.^2)/Nint);

% The envelope is Rayleigh distributed
Rayleigh=(real(h).^2+imag(h).^2).^.5;

plot(Rayleigh)
