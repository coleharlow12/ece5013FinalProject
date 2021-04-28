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

fclose(fid);
