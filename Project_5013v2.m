clc
clear

%% Set Physical Constants
c = 3e8;

%% Set Radar Parameters
fc = 10.5e9; %Radar center frequency (Hz)
Bo = 120e6; %Radar Sweep Bandwidth (Hz)
Tau = 80e-6; %Ramp period (s)
prf = 1e3; %Ramp Rate (Hz)
Pt = 10^((20-30)/10);

Np = 64; %Number of pulses in CPI
fs = 250e6; %Sampling Rate (samples/sec)

beamW = 60*pi/180; %Antenna Beamwidth (radians)
Ga = 10^(3/10); %Antenna Gain for both RX/TX (linear)

NF = 10^(8/10); %Noise figure of system (linear)

%% Target Parameters
tarR = [30]; %target range (m)
velX = [10]; %Velocity Perpendicular to Boresight
velY = [0]; %Velocity Parellel to Boresight
angAz = [-10]*pi/180; %Azimuth angle ([degrees] -> radians)
rcs = [3]; %Target RCS

%             y
%             ^   O target
%             | °/
%             | /
%             |/
%      RADAR V V--------> X
%            U D
%% Create The TX Data
Tp = 1/prf; %Pulse Period
t = (1/(fs)):(1/(fs)):Tp; %Sampled Data times

stU = zeros(size(t));
stD = zeros(size(t));

for it=1:length(t)
    if (t(it) <= Tau)
        %Amplitude is not considered here because it is considered in friis
        %scaling
        stU(it) = Pt*exp(1i*(2*pi*t(it)+pi*Bo/Tau*t(it)^2)); %Up Ramp signal
        stD(it) = Pt*exp(1i*(2*pi*t(it)-pi*Bo/Tau*t(it)^2)); %Down ramp signal
    end
end

% Plots the transmitted signal for 1 pulse
figure(1)
plot(t,real(stU),'r-',t,real(stD),'b-','LineWidth',2)
title('Single Pulse')
xlabel('time (s)')
ylabel('Tx Signal')
legend('Up Ramp','Down Ramp')
set(gca,'fontsize',18)


stU = repmat(stU,1,Np); %Time Domain Data for N pulses
stD = repmat(stD,1,Np); %Time Domain Data for N pulses
t = (1:length(stD))*(1/fs);

%% Create Rx Data Using Stop and Hop Approximation
%Calculate initial x and y coordinates
Tp = 1/prf; %Pulse Period
lambda = c/fc; %Calculates the wavelength
d = lambda/2; %The spacing between the two tx antennas

x0 = tarR*sin(angAz); %Calculates the starting X location
y0 = tarR*cos(angAz); %Calculates the starting Y location

rxComb = [];

for iP = 1:Np    
    t=Tp*(iP-1);%The time the pulse starts is equal to the pulse period times the pulse number
    x = x0+velX*t; %X location at the start of the pulse
    y = y0+velY*t; %Y location at the start of the pulse
    r = sqrt(x^2+y^2); % Radius at the start of the pulse
    Az = asin(x/r); %Azimuth location. Needed for antenna phase difference
    td = d*sin(Az);  %Delay between the up and down trasnsmitter
    
    %For positive azimuth the down (right) transmitter will have a positive
    %phase shift
    
    trnd = 2*r/c;
    t = (1/(fs)):(1/(fs)):Tp; %Sampled Data Times for one pulse repition
    for it = 1:length(t)
        %signal is zero before the pulse starts/after the pulse ends
        if (t(it)< trnd || t(it)>(trnd+Tau) )
            rxsigU(it) = 0;
            rxsigD(it) = 0;
        else
            rxsigU(it) = exp(1i*(2*pi*(t(it)-trnd)+pi*Bo/Tau*(t(it)-trnd)^2)); %Up Ramp signal
            rxsigD(it) = exp(1i*(2*pi*(t(it)-trnd+td)+pi*Bo/Tau*(t(it)-trnd+td)^2)); %Down Ramp signal
        end
    end
    FriisScale = Pt*Ga^2*lambda^2*rcs/((4*pi)^3*r^4);
    rxSigU = FriisScale*rxsigU;
    rxSigD = FriisScale*rxsigD;
    
    rxComb = [rxComb, (rxSigU+rxSigD)]; %Concatenate new data
end

%% Add Noise

