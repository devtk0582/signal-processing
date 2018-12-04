
% GPS interference cancellation

% Parameters
n = 2; % Number of antennas
N_bits = 100; % Number of bits transmitted
N_chip = 10; % Number of chips per bit
T = N_bits*N_chip; % Measurement interval (number of time samples of received signal)
lambda = 0.2; % Wavelength in meter of 1.5 GHz GPS signal
d = 0.5*lambda; % Distance between antennas (typical)
sigma_noise = 0.1; % Noise voltage relative to signal voltage
sigma_alien = 5.0; % Alien noise voltage relative to signal voltage
theta_gps = 30*pi/180; % GPS location in radians
theta_alien = 32*pi/180; % Alien location in radians


% Array responses
array_gps = exp(j*2*pi*d/lambda*sin(theta_gps)*[0:n-1]');
array_alien = exp(j*2*pi*d/lambda*sin(theta_alien)*[0:n-1]');

% Channel gains
h = array_gps;
g = array_alien;

% Signals transmitted
chip_seq = ones(1,N_chip);
x_bit = sign(randn(1,N_bits)); % GPS bits (binary data)
x = filter(chip_seq,[1],reshape([x_bit ; zeros(N_chip-1,N_bits)],1,T)); % GPS signal transmitted
a = sigma_alien*randn(1,T); % Alien interference (random Gaussian, e.g.)
noise = sigma_noise*randn(n,T); % Noise in the 'n' receiver antennas

% Received signal
y = h*x + g*a + noise;


%*********    GPS receiver processing  **************

w_match = h/norm(h)^2; % Matched filter
Proj = (eye(n) - g*g'/norm(g)^2); % Projection matrix
w_cancel = Proj*h/(h'*Proj*h);% Interference canceller


% Result of GPS receiver processing
z_match = w_match'*y;
z_cancel = w_cancel'*y;

% Plot percent-error between x (GPS signal) and z (processed data)
figure(1)
clf
subplot(2,1,1)
plot([0:T-1]*1e-6,x,'b');
hold on;
plot([0:T-1]*1e-6,y(1,:),'r')
xlabel('Time in seconds')
legend('GPS signal','signal received at antenna 1')
subplot(2,1,2)
plot([0:T-1]*1e-6,real(z_match),'g'); % Matched filter result
hold on;
plot([0:T-1]*1e-6,real(z_cancel),'r'); % Interference canceller result
xlabel('Time in seconds')
legend('GPS signal filtered with w\_match','GPS signal filtered with w\_cancel')
