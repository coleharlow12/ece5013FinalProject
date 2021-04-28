fid=fopen('measured_data.bin','rb');

M=256; % number of pulses in one CPI
N_cpi= 32; % number of CPIs in data


for ii=1:N_cpi
    phase_history1=zeros(301,M);
    phase_history2=zeros(301,M);
    
    for jj=1:M
          phase_history1(:,jj)=fread(fid,301,'double');
          phase_history1(:,jj)=phase_history1(:,jj)+1i*fread(fid,301,'double');

          phase_history2(:,jj)=fread(fid,301,'double');
          phase_history2(:,jj)=phase_history2(:,jj)+1i*fread(fid,301,'double');
    end
  
   
   %only Range bins >=131 have target information, earlier range bins are
   %due to antenna coupling, TX leaking into RX, delays in the radar.
   % row range-bins 0:170, columns :pulse_no 1:M
   
   phase_history1=phase_history1(131:end,:);
   phase_history2=phase_history2(131:end,:);
    
   % phase_history matrix: row range-bins 0:170, columns :pulse_no 1:M
   % Construct Range Doppler map from the two antennas and display.
   
   % Consider canceling stationary clutter (two/three pulse canceller, or 
   % subtracting average (across slow time) range-profile from each return 
end

% Two pulse canceller
for ii=2:M
   phase_history1_tp(:,ii)= phase_history1(:,ii) - phase_history1(:,ii-1);
   phase_history2_tp(:,ii)= phase_history2(:,ii) - phase_history2(:,ii-1);
end

% Three pulse canceller
for ii=2:M
   phase_history1_thrp(:,ii)= phase_history1_tp(:,ii) - phase_history1_tp(:,ii-1);
   phase_history2_thrp(:,ii)= phase_history2_tp(:,ii) - phase_history2_tp(:,ii-1);
end

figure(1)
subplot(1,3,1);
imagesc(abs(phase_history1))
xlabel('Pulse No');
ylabel('Range bin');
title('Phase History 1');

subplot(1,3,2);
imagesc(abs(phase_history1_tp))
xlabel('Pulse No');
ylabel('Range bin');
title('Phase History 1: 2-Pulse Cancelled');

subplot(1,3,3);
imagesc(abs(phase_history1_thrp))
xlabel('Pulse No');
ylabel('Range bin');
title('Phase History 1: 3-Pulse Cancelled');

figure(2)
subplot(1,3,1);
imagesc(abs(phase_history2))
xlabel('Pulse No');
ylabel('Range bin');
title('Phase History 2');

subplot(1,3,2);
imagesc(abs(phase_history2_tp))
xlabel('Pulse No');
ylabel('Range bin');
title('Phase History 2: 2-Pulse Cancelled');

subplot(1,3,3);
imagesc(abs(phase_history2_thrp))
xlabel('Pulse No');
ylabel('Range bin');
title('Phase History 2: 3-Pulse Cancelled');

fclose(fid);