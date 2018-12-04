
% Simulation of a fading multipath (mobile) channel.
% Based on Shanmugan, et. al. ('Simulation of Comm. Systems')

%function [h,L,ds] = mobilechan(tau0,T,fd,fs,Nint)

% Sample parameters:  tau0 = 5 microsec, T = 1 microsec, fd = 30 Hz,
%      fs = 1000 Hz, Nint = 2^15 (i.e. 2^15/1000 seconds)
% fd = doppler rate. e.g. 30 Hz
% fs = channel sampling rate (to track the channel variations). e.g. 1000 Hz
% Decay const. of exponentially decaying ISI  = tau0 . e.g. 5 microsec.
% Length of channel response (decay to 1%) = Tm = 4.6*tau0 
% Sampling rate of signal (=bandwidth of signal) = W = 1/T . e.g. 1 MHz 
% Length of channel = L = Tm*W+1
% Using quasi-stationary approx, number of channel snapshots = 2*Nint. e.g. 1024
% h = channel impulse response. Nint snapshots of it, one in each column
% ds = delay-spread of the channel
clear
tau0 = 5; % microsec
T = 1; % microsec
fd = 30; %Hz
fs = 1000; %Hz
Nint = 2^15;

W = 1/T;
Tm = 4.6*tau0;
L = ceil(Tm*W+1);  %L>1 means frequency selective fading
ds = Tm;

tempy=[];
sqrtpsd=[1./(1-([-fd/fs*Nint+.5:1:fd/fs*Nint-.5]/(fd/fs*Nint)).^2).^.25];
for (ii=1:L)
    ampl=randn(1,floor(2*Nint*fd/fs))+j*randn(1,floor(2*Nint*fd/fs));
    fade =  ifft(ampl.*sqrtpsd,2*Nint);
    tempy(:,ii) = fade';
end
% Normalization of the channel impulse response
tempy = sqrt(2*Nint)*tempy./sqrt(sum(sqrtpsd.^2)/Nint);

% scale by exponential PDP (Power delay profile)
alpha = 1/sqrt(sum(exp(-[0:L-1]/(tau0*W))));
h = conj(tempy*diag(alpha^0.5*exp(-0.5*[0:L-1]/(tau0*W))))';

% The envelope is Rayleigh distributed
ii=15;  % path index
Rayleigh=(real(h(ii,:)).^2+imag(h(ii,:)).^2).^.5;

plot(Rayleigh)