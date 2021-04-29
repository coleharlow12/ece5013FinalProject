clc
close all;
clear all;
clear
%% Parameters for MISO Radar System %%%
fc = 10.5e9;        % center frequency = 10.5 GHz
BW = 125e6;         % total system bandwidth = 125 MHz
fp = 1e3;           % pulse repetition frequency = 1 kHz
Tp = 1/fp;          % pulse repetition interval = 1 ms
Np = 64;            % number of pulses
c = 3e8;
lambda = c/fc;

beta = 120e6;       % sweep bandwidth = 120 MHz
tau = 80e-6;        % pulse width = 80 usec


% Here I assume SNR=0 dB, you have to simulate correct Amplitude and Noise
% Levels, using radar range equation, etc


%% Parameters for Sampled System %%%
fs = 250e6;         % sample rate = 250 Msamples/second
Ts = 1/fs;          % sample period

%% Load Data
fid=fopen('measured_data.bin','rb');
M=256; % number of pulses in one CPI
N_cpi= 32; % number of CPIs in data

% Loops over each CPI
for iCPI=1:N_cpi
    phase_history1=zeros(301,M);
    phase_history2=zeros(301,M);
    
    for jj=1:M
          %Assume 1 is up
          phase_history1(:,jj)=fread(fid,301,'double');
          phase_history1(:,jj)=phase_history1(:,jj)+1i*fread(fid,301,'double');
          
          %Assume 2 is down
          phase_history2(:,jj)=fread(fid,301,'double');
          phase_history2(:,jj)=phase_history2(:,jj)+1i*fread(fid,301,'double');
    end
  
   
   %only Range bins >=131 have target information, earlier range bins are
   %due to antenna coupling, TX leaking into RX, delays in the radar.
   % row range-bins 0:170, columns :pulse_no 1:M
   
   % Potentially Output of Matched filters
   phase_history1=phase_history1(131:end,:);
   phase_history2=phase_history2(131:end,:);
   
   % phase_history matrix: row range-bins 0:170, columns :pulse_no 1:M
   % Construct Range Doppler map from the two antennas and display.
   
%    %Two Pulse Canceling
%    for iP=2:M
%        phase_history1_tp(:,iP)= phase_history1(:,iP) - phase_history1(:,iP-1);
%        phase_history2_tp(:,iP)= phase_history2(:,iP) - phase_history2(:,iP-1);       
%    end
%    
%    % Three pulse canceller
%     for iP=3:M
%        phase_history1_thrp(:,iP)= phase_history1_tp(:,iP) - phase_history1_tp(:,iP-1) - phase_history1_tp(:,iP-2);
%        phase_history2_thrp(:,iP)= phase_history2_tp(:,iP) - phase_history2_tp(:,iP-1) - phase_history2_tp(:,iP-2);
%     end
    
    % Average Across Slow Time
    meanPhase1 = mean(phase_history1,2);
    meanPhase2 = mean(phase_history2,2);
    
    phase_history1 = phase_history1-repmat(meanPhase1,1,M);
    phase_history2 = phase_history2-repmat(meanPhase2,1,M);
    
    % Create Window
    win = tukeywin(size(phase_history1,2),0);
    winR = repmat(win',size(phase_history1,1),1);
    
    phase_history1_win = winR.*phase_history1;
    phase_history2_win = winR.*phase_history2;
   
   % FFT Along 
   rangedopplerUp=fftshift(fft(phase_history1_win,[],2),2);
   rangedopplerDown=fftshift(fft(phase_history2_win,[],2),2);
   
   Ntau = size(rangedopplerUp,1);
   Np = size(rangedopplerUp,2);
   
   % Plot Range Doppler
   nugrid=1/Np*(-Np/2:1:(Np/2)-1); %The nugrid contains all the doppler frequencies
   taugrid=((0:Ntau)+130)*Ts*c/2;               %The taugrid contains all the delay frequencies
   
   figure(1)
   imagesc(taugrid,nugrid,abs(rangedopplerUp.'));
   xlabel('Range (m)'); ylabel('Normalized Frequency (sec)');title('Up Ramp')
   set(gca,'fontsize',18)
   % Consider canceling stationary clutter (two/three pulse canceller, or 
   % subtracting average (across slow time) range-profile from each return
   
   %% Find Range / Velocity / Angle For One CPI
    % Find Range
%     [~,d1] = max(abs(sum(rangedopplerUp,2))); %Sum integrates coherently across pulses in dimension 1
%     [~,d2] = max(abs(sum(rangedopplerDown,2))); 

    maximum = max(max(rangedopplerUp));
    [x,y] = find(rangedopplerUp == maximum);
    d1=x;
    maximum = max(max(rangedopplerDown));
    [x,y] = find(rangedopplerDown == maximum);
    d2=x;
    
    range_d1(iCPI) = taugrid(d1);
    range_d2(iCPI) = taugrid(d2);
    
    % Find Velocity
    [~,s1] = max(abs(rangedopplerUp(d1,:)));
    [~,s2] = max(abs(rangedopplerDown(d2,:)));
    
    radial_speed1(iCPI) = nugrid(s1)*fp*c/(2*fc); %converts to real frequency
    radial_speed2(iCPI) = nugrid(s2)*fp*c/(2*fc); %converts to real frequency
    
    % Find Angle
    angU = angle(sum(rangedopplerUp(d1,:)));
    angD = angle(sum(rangedopplerDown(d2,:)));
    difAng = (angU-angD);
    if(difAng<-1)
        difAng=difAng+2*pi;
    end
    d = lambda/2;
    ang(iCPI) = asin(difAng*lambda/(2*pi*d))*180/pi;
   
end
figure(2)
hold on
% plot(range_d1)
% plot(range_d2)
plot((range_d1+range_d2)/2)
hold off
title('range')

figure(3)
hold on 
% plot(radial_speed1)
% plot(radial_speed2)
plot((radial_speed1+radial_speed2)/2)
hold off
title('speed')

figure(4)
plot(abs(ang),'r-')
title('angle')

fclose(fid);