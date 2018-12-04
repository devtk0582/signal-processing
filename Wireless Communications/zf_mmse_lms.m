
% Equalizer

T = 5000;  % total number of bits (including both training phase and data phase
Tlen = 2500;  % number of bits in training phase
Ex = 1;  % bit energy
SNR_mfb_dB = 20; % Matched filter bound in dB. Try 20 and 40 dB
mmse_len = 40;  % number of taps in MMSE-LE
mu = 0.005;  % step size in LMS algorithm
delay = floor(mmse_len/2);  % delay the bitstream to achieve causality

%hh= poly([0.8 0.6 -0.7]) % Non-singular channel
hh= poly([1 0.9 -0.95]) % Singular channel with zero=1, which is on unit circle, unstable pole for 1/H(z^-1)
%1.0000   -0.9500   -0.9050    0.8550
%z_k = x_k -0.9500*x_{k-1} -0.9050*x_{k-2}+ 0.8550*x_{k-3} + n_k

SNR_mfb = 10^(SNR_mfb_dB/10);
sigma_n = sqrt(norm(hh)^2*Ex/SNR_mfb);  % standard deviation of noise

% Channel
bits = sign(randn(T,1)); % bipolar signal
xx = sqrt(Ex)*bits;  % transmit signal
zz = conv(xx,hh); % output of linear channel
zz = zz+sigma_n*(randn(size(zz))+j*randn(size(zz)));  %AWGN channel

% ZF-LE
v_zfle = filter([1],[hh],zz);  % output of ZF-LE

% MMSE-LE
ww = zeros(mmse_len,1); % initialize tap coefficients
if (mmse_len-delay>2)
   vv = zeros(mmse_len-delay-1,1); % initialize output of MMSE-LE (for causality)
   ee = zeros(mmse_len-delay-1,1); % initialize error (for causality)

end
%Training phase
for (ii=mmse_len-delay:Tlen)  % t= mmse_len-delay is the starting epoch (the delay is to achieve causality)
   vv(ii) = ww'*zz(ii+delay:-1:ii+delay-mmse_len+1);  % output of MMSE-LE
   ee(ii) = vv(ii) - xx(ii); % error ek = vk - xk
   ww = ww-mu*conj(ee(ii))*zz(ii+delay:-1:ii+delay-mmse_len+1);  % update tap coefficients
end
% Data phase
for (ii=Tlen+1:T-mmse_len) 
   %  ii should not be great than T-mmse_len, since we need to buffer mmse_len number of signals in the tapped delay line.
   vv(ii) = ww'*zz(ii+delay:-1:ii+delay-mmse_len+1);
   ee(ii) = vv(ii) - xx(ii); % ek = vk - xk
end


% Plot performance
figure(1)
subplot(2,2,1)
stem(xx,'k')
subplot(2,2,2)
stem(real(zz),'b');
subplot(2,2,3)
stem(real(v_zfle),'r')
subplot(2,2,4)
stem(real(vv),'m')
figure(2)
subplot(3,1,1)
plot(real(zz(1:length(xx))-xx),'b')
title('unequalized error')
subplot(3,1,2)
plot(real(v_zfle(1:length(xx)))-xx,'r')
title('zf-le error')
subplot(3,1,3)
plot(real(ee),'m')
title('mmse-le error')

