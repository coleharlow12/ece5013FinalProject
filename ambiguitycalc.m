fs=500e6; %Sampling Frequency
Ts=1/fs; %Sampling period
tau=80e-6; %Chirp length
b=120e6; %Bandwidth

fp=1e3; %Pulse repetition frequency
Tp=1/fp; %Pulse repetition period

Np=64; %Number of pulses

t=0:Ts:(tau-Ts); %Time data

za=exp(1i*pi*(b/tau).*(t-(tau/2)) .^2); %Pulse data

nugridsize=64; %Number of different doppler freqeuncies to try in 1 direction
% 
% M=length(t)-1;
% Mi=floor(M/taugridsize);
% Mi=0.1*Tp*fs;
% 
% taugridsize=floor(M/Mi);
% 
% Mind=-taugridsize*Mi:Mi:taugridsize*Mi; %centered at 0
% Mind=Mind + M+1; %centered at M+1 <=zero lag

maxlag=floor(10/b*fs); %Maximum time lag? 1/b is pulse width? fs is the sampling frequency. So this is max lag in samples
ambiguity=zeros(2*nugridsize+1,2*maxlag+1); %Grid of samples to test for ambiguity

for i=-nugridsize:1:nugridsize
    
    fd=  i*2*1/tau/nugridsize; % fd =-fp/2  to fp/2
    za_fd=za.*exp(j*2*pi*fd*t); % This is the delayed doppler shifted pulse
    
[acor,lag]=xcorr(za_fd,za, maxlag); %Does cross-correlation across all possible lags

%2M+1 points in acor, where M=lenght(t)-1, sample taugridsize*2+1 points

ambiguity(i+nugridsize+1,:)=acor; %stores results for the given doppler frequency and all lags
end

%% Plot Results
nugrid=linspace(-2/tau,2/tau,2*nugridsize+1)/1000;
taugrid=(-maxlag:maxlag)*Ts*1e6;

figure(1); mesh(taugrid, nugrid, abs(ambiguity))
 
 xlabel('tau (\mu sec)')
 ylabel('\nu (kHz)')
 
figure(2); contour(taugrid, nugrid, abs(ambiguity))
 
 xlabel('tau (\mu sec)')
 ylabel('\nu (kHz)')
 
%%% Creates a rectangular pulse of width tau from 0<t<tau %%%
function p = rpulse(t,tau)
    p = (t<=tau)&(t>=0);
end

