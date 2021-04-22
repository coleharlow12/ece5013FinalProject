clc
close all;
clear all;

%% Physical Constants %%%
c = 3e8;            % speed of light

%% Parameters for MISO Radar System %%%
fc = 10.5e9;        % center frequency = 10.5 GHz
lambda = c/fc;      % wavelength of radar system
BW = 125e6;         % total system bandwidth = 125 MHz
fp = 1e3;           % pulse repetition frequency = 1 kHz
Tp = 1/fp;          % pulse repetition interval = 1 ms
Np = 64;            % number of pulses


beta = 120e6;       % sweep bandwidth = 120 MHz
tau = 80e-6;        % pulse width = 80 usec

Ga = 10^(3/10); %Antenna Gain for both RX/TX (linear)
F = 10^(8/10); %Noise figure of system (linear)
Pt = 10^((20-30)/10); %Transmit Power
rcs = 2;

% Here I assume SNR=0 dB, you have to simulate correct Amplitude and Noise
% Levels, using radar range equation, etc


%% Parameters for Sampled System %%%
fs = 250e6;         % sample rate = 250 Msamples/second
Ts = 1/fs;          % sample period

%% Time Vector %%%

t = 0:Ts:Np*Tp-Ts;     % time vector (Np+1) Tp long, 

%% Parameters for Target %%%
R0 = 30;            % R0 is initially 30 meters
theta = 10;         % theta = azimuth angle, initially at 10 degrees
v = 10;             % vertical  velocity is 10 m/s

%% Received signals
za1 = zeros(size(t));   % received signal from transmitter 1
za2 = zeros(size(t));   % received signal from transmitter 2
za = zeros(size(t));    % Combined RX signal

%% Simulate noise
kTo = 4*10^-21;
Pnoise = kTo*BW*F;
sigma_n=sqrt(Pnoise/2); % You should calculate using k T0 BW F
noise = sigma_n*(randn(1,length(t))) + 1i*(sigma_n*(randn(1,length(t))));

tx1=[ 0, lambda/4]; tx2=[0, -lambda/4]; rx=[0,0];

% calculate azimuth angle, range, and received signals at each location
for k=0:Np-1
    target=[R0*cosd(theta) R0*sind(theta)-k*Tp*v];
    
    Rup1 = norm(tx1-target); %tx1 to target distance
    Rup2= norm(tx2-target); %tx2 to target distance
    Rdown = norm(rx-target); %rx to target distance
    
   Ac= sqrt(Pt*Ga^2*lambda^2*rcs/((4*pi)^3*Rup1^4)); % You should  compute for each pulse using radar 
         % range eqn, OR ccompute once for some nominal range 
         % independent of k
    
    %Signal is 0 everywhere except the current pulse so addition is across
    %entire vector
    %up ramp
%     za1 = za1 + Ac*(rpulse(t-k*Tp-(Rup1+Rdown)/c,tau)) .* ...
%         (exp(-1i*(2*pi/lambda)*(Rup1+Rdown))) .* ...
%         (exp(1i*pi*(beta/tau).*(t-(tau/2)-(k*Tp+(Rup1+Rdown)/c )) .^2));
    
    za1 = za1 + Ac*(rpulse(t-k*Tp-(Rup1+Rdown)/c,tau) .* ...
        (exp(1i*2*pi*fc*(t-k*Tp-(Rup1+Rdown)/c))) .* ...
        (exp(1i*pi*(beta/tau).*(t-k*Tp-(Rup1+Rdown)/c).^2)));
    
    %down ramp
%     za2 = za2 + Ac*(rpulse(t-k*Tp-(Rup2+Rdown)/c,tau)) .* ...
%         (exp(-1i*(2*pi/lambda)*(Rup2+Rdown))) .* ...
%         (exp(-1i*pi*(beta/tau).*(t-(tau/2)-(k*Tp+(Rup2+Rdown)/c )) .^2));
    
   za2 = za2 + Ac*(rpulse(t-k*Tp-(Rup2+Rdown)/c,tau) .* ...
        (exp(1i*2*pi*fc*(t-k*Tp-(Rup2+Rdown)/c))) .* ...
        (exp(-1i*pi*(beta/tau).*(t-k*Tp-(Rup2+Rdown)/c).^2)));

end
za=za1+za2+noise;

figure(5)
plot(t,za)

%%

%Prepare match filter for Tx1
tm=0:Ts:tau-Ts;                             %time range
%hu=exp(1i*pi*(beta/tau).*(tm-(tau/2)) .^2); %TX pulse up
hu = exp(1i*(2*pi*fc*tm+pi*beta/tau*tm.^2)); %Up Ramp signal
hu=conj(fliplr(hu));                        %conjugate and flip the time.

%Prepare match filter for TX2
%hd = exp(-1i*pi*(beta/tau).*(tm-(tau/2)) .^2); %TX pulse up
hd = exp(1i*(2*pi*fc*tm-pi*beta/tau*tm.^2)); %Down ramp signal
hd = conj(fliplr(hd));


receivearray=zeros(Np,tau*fs+500); %20,000 is tau*fs
for k=0:Np-1
    receivearray(k+1,:)=za(k*uint64(Tp*fs)+1:round(k*uint64(Tp*fs)+tau*fs+500));%Essentially takes each chirp of data?
end

UpArray=zeros(Np,501);
DownArray=zeros(Np,501);
for k=0:Np-1
    UpArray(k+1,:)=conv(receivearray(k+1,:),hu,'valid'); %Convolution flips then scans
    DownArray(k+1,:)=conv(receivearray(k+1,:),hd,'valid'); %Convolution flips then scans
end

taugrid=(0:500)*Ts;
figure(1);imagesc(taugrid,1:64,abs(UpArray));
figure(2);imagesc(taugrid,1:64,abs(DownArray));
xlabel('Delay'); ylabel('Pulse No');

rangedopplerUp=fftshift( fft(UpArray,[],1),1);
rangedopplerDown=fftshift( fft(DownArray,[],1),1);

nugrid=1/Np*(-Np/2:1:(Np/2)-1);
taugrid=(0:500)*Ts;
figure(3);imagesc(taugrid,nugrid, abs(rangedopplerUp))
figure(4);imagesc(taugrid,nugrid, abs(rangedopplerDown))

xlabel('Delay (sec)'); ylabel('Normalized Frequency (sec)');

%Things to try, window the match filter, window the pulses, add ground
%return, use MTI cancelling, zeropad FFT, % label the range doppler map in 
% velocity(m/sec) and Range( meters)


%%% Creates a rectangular pulse of width tau from 0<t<tau %%%
function p = rpulse(t,tau)
    p = (t<=tau)&(t>=0);
end

