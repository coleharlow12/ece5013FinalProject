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
% za1 = zeros(size(t));   % received signal from transmitter 1
% za2 = zeros(size(t));   % received signal from transmitter 2
% za = zeros(size(t));    % Combined RX signal

%% Simulate noise
kTo = 4*10^-21;
Pnoise = kTo*BW*F;
sigma_n=sqrt(Pnoise/2); % You should calculate using k T0 BW F
noise = sigma_n*(randn(1,length(t))) + 1i*(sigma_n*(randn(1,length(t))));

%% Simulate the Signal
tx1=[ 0, lambda/4]; tx2=[0, -lambda/4]; rx=[0,0];
trueAz= zeros(Np,1); %Variable to hold the true azimuth angle of the target
% 
% %calculate azimuth angle, range, and received signals at each location
% for k=0:Np-1
%     target=[R0*cosd(theta) R0*sind(theta)-k*Tp*v] %Calculates current target location
%     
%     Rup1 = norm(tx1-target); %tx1 to target distance
%     Rup2= norm(tx2-target); %tx2 to target distance
%     Rdown = norm(rx-target); %rx to target distance
%     trueAz(k+1) = atan(target(2)/target(1))*180/pi; %add data to true azimuth location
% 
%    %Scaling due to Friis Transmission Equation
%    Ac= sqrt(Pt*Ga^2*lambda^2*rcs/((4*pi)^3*Rup1^4)); 
%     
%     %Signal is 0 everywhere except the current pulse so addition is across
%     %entire vector
%     %up ramp
%     za1 = za1 + Ac*(rpulse(t-k*Tp-(Rup1+Rdown)/c,tau)) .* ...
%         (exp(-1i*(2*pi/lambda)*(Rup1+Rdown))) .* ...
%         (exp(1i*pi*(beta/tau).*(t-(tau/2)-(k*Tp+(Rup1+Rdown)/c )) .^2));
%     
%     %down ramp
%     za2 = za2 + Ac*(rpulse(t-k*Tp-(Rup2+Rdown)/c,tau)) .* ...
%         (exp(-1i*(2*pi/lambda)*(Rup2+Rdown))) .* ...
%         (exp(-1i*pi*(beta/tau).*(t-(tau/2)-(k*Tp+(Rup2+Rdown)/c )) .^2));
% end

%Loops across pulses
t = ((1/(fs)):(1/(fs)):Tp); %Timing for a single pulse
rxSigU = zeros(size(t)); %Stores signal for each pulse
rxSigD = zeros(size(t)); %Stores signal for each pulse

for iP = 1:Np   
    target=[R0*cosd(theta) R0*sind(theta)-(iP-1)*Tp*v] %Calculates current target location
    Rup1 = norm(tx1-target); %tx1 to target distance
    Rup2= norm(tx2-target); %tx2 to target distance
    Rdown = norm(rx-target); %rx to target distance
    trueAz(iP) = atan(target(2)/target(1))*180/pi; %add data to true azimuth location
    tof = (Rup1+Rdown)/c;
    
    Ac= sqrt(Pt*Ga^2*lambda^2*rcs/((4*pi)^3*Rup1^4)); 
    for it = 1:length(t)
        %signal is zero before the pulse starts/after the pulse ends
        if (t(it)< (tof) || t(it)>(tof+tau) )
            rxSigU(it) = 0;
            rxSigD(it) = 0;
        else
            % Accounts for constant range related delay
            % then accounts for ramping based frequency
            rxSigU(it) = Ac*(exp(-1i*(2*pi/lambda)*(Rup1+Rdown)) * ...
                         exp(1i*pi*(beta/tau)*(t(it)-tau/2-(Rup1+Rdown)/c )^2));          
            rxSigD(it) = Ac*(exp(-1i*(2*pi/lambda)*(Rup2+Rdown)) * ...
                         exp(-1i*pi*(beta/tau)*(t(it)-tau/2-(Rup2+Rdown)/c )^2));
        end
    end
    
    if iP == 1
        za1=rxSigU;
        za2=rxSigD;
    else
        za1 = [za1,rxSigU];
        za2 = [za2,rxSigD];
    end
end

%%
t = 0:Ts:Np*Tp-Ts;
%Receiver receives both the up and down ramps
za = za1+za2;

%Create a plot of the signal in time
figure(4)
plot(t,abs(za))
xlabel('time (s)') 
ylabel('abs(rx)')
set(gca,'fontsize',18)

za=za+noise; %Combines the signal and the noise

%% Match Filter/Frequency
%Prepare match filter for Tx1
tm=0:Ts:tau-Ts;                             %time range
hu=exp(1i*pi*(beta/tau).*(tm-(tau/2)) .^2); %TX pulse up
hu=conj(fliplr(hu));                        %conjugate and flip the time.

%Prepare match filter for TX2
hd = exp(-1i*pi*(beta/tau).*(tm-(tau/2)) .^2); %TX pulse down
hd = conj(fliplr(hd));

Ntau = 500; %The number of delays(ranges) to test out
receivearray=zeros(Np,tau*fs+Ntau); %20,000 is tau*fs the Ntau padding is so that we can do convolution without samples not overlapping

%Prepares rx data for received filtering
%Receive array is Np x Ntau
for k=0:Np-1
    receivearray(k+1,:)=za(k*uint64(Tp*fs)+1:round(k*uint64(Tp*fs)+tau*fs+Ntau)); %Pads each pulse with Ntau points so that data
end

%Holds results of matched filtering. Np x Ntau+1.
%the 
UpArray=zeros(Np,Ntau+1); 
DownArray=zeros(Np,Ntau+1);
for k=0:Np-1
    UpArray(k+1,:)=conv(receivearray(k+1,:),hu,'valid'); %Convolution flips then scans
    DownArray(k+1,:)=conv(receivearray(k+1,:),hd,'valid'); %Convolution flips then scans
end

%% Plot Matched Filter Results
% Tau Grid has Ntau samples on it
taugrid=(0:Ntau)*Ts; %Different delays which are tested each delay corresponds to a range of 2*R/c
figure(1);imagesc(taugrid,1:64,abs(UpArray));
xlabel('Delay'); ylabel('Pulse No');

figure(2);imagesc(taugrid,1:64,abs(DownArray));
xlabel('Delay'); ylabel('Pulse No');

%% Plots the Range Doppler Maps
%Takes the fft of each column (i.e across dim1 which is rows)
%So this gives us an fft of length Np for each of the Ntau delays/ranges
%i.e we are taking fft across the number of pulses. As the object moves
%towards us there is a constant phase shift between each measurement. The
%frequency of this phase shift corresponds to the doppler frequency
rangedopplerUp=fftshift( fft(UpArray,[],1),1);
rangedopplerDown=fftshift( fft(DownArray,[],1),1);


nugrid=1/Np*(-Np/2:1:(Np/2)-1); %The nugrid contains all the doppler frequencies
taugrid=(0:Ntau)*Ts;            %The taugrid contains all the delay frequencies
figure(3);imagesc(taugrid,nugrid, abs(rangedopplerUp))
xlabel('Delay (sec)'); ylabel('Normalized Frequency (sec)');title('Up Ramp')
set(gca,'fontsize',18)

figure(4);imagesc(taugrid,nugrid, abs(rangedopplerDown))
xlabel('Delay (sec)'); ylabel('Normalized Frequency (sec)');title('Down Ramp')
set(gca,'fontsize',18)

%% Find Range / Velocity / Angle For One CPI
% Find Range
[~,d1] = max(abs(sum(rangedopplerUp,1))); %Sum integrates coherently across pulses in dimension 1
[~,d2] = max(abs(sum(rangedopplerDown,1))); 

range_d1 = taugrid(d1)*c/2;
range_d2 = taugrid(d2)*c/2;

% Find Velocity
[~,s1] = max(abs(rangedopplerUp(:,d1)));
[~,s2] = max(abs(rangedopplerDown(:,d2)));

radial_speed1 = nugrid(s1)*fp; %converts to real frequency
radial_speed1 = radial_speed1*c/(2*fc); %frequency -> velocity

radial_speed2 = nugrid(s2)*fp; %converts to real frequency
radial_speed2 = radial_speed2*c/(2*fc); %frequency -> velocity

% Find Angle
angU = angle(sum(rangedopplerUp(:,d1)));
angD = angle(sum(rangedopplerDown(:,d2)));
difAng = (angU-angD);
if(difAng<-1)
    difAng=difAng+2*pi;
end
d = lambda/2;
ang = asin(difAng*lambda/(2*pi*d))*180/pi;
