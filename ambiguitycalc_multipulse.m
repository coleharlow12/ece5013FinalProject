clc
clear

%% Set Parameters for Pulses
fs=250e6; %Sampling Frequency
Ts=1/fs; %Sampling Period
tau=80e-6; %PulseWidth
b=120e6; %Bandwidth

fp=1e3; %Pulse Repition frequnecy 
Tp=1/fp; %Pulse and quiet time length

Np=16; %Number of Pulses

t=0:Ts:Np*Tp-Ts; %Time Samples

za=zeros(size(t)); %TX signal

%% Create Full TX Signal
for k=0:Np-1,
    za=za+exp(1i*pi*(b/tau).*(t-(tau/2)-k*Tp) .^2).*rpulse(t-k*Tp,tau);
end;

figure(1)
plot(t,abs(za))
ylabel('Absolute Value of Signal')
xlabel('time (s)')

%% Locations to Sample Ambiguity Function
M=length(t)-1; %1 less than total number of samples?

Mi=0.1*Tp*fs; %want to make sure when we subsample we hit the multiples of Tp

taugridsize=floor(M/Mi); %Number of delays we test
nugridsize=32; %Number of doppler shifts to test

Mind=-taugridsize*Mi:Mi:taugridsize*Mi; %centered at 0
Mind=Mind + M+1; %centered at M+1 <=zero lag

%%

ambiguity=zeros(2*nugridsize+1,2*taugridsize+1);


for i=-nugridsize:1:nugridsize;
    
    fd=  i* 2* fp/nugridsize; %fd =-fp/2  to fp/2
    za_fd=za.*exp(j*2*pi*fd*t); %the shifted frequency data

    tic;
[acor,lag]=xcorr(za_fd,za);
toc;
%2M+1 points in acor, where M=lenght(t)-1, sample taugridsize*2+1 points

ambiguity(i+nugridsize+1,:)=acor(Mind); % store after subsampling in delay

end



nugrid=linspace(-2*fp,2*fp,2*nugridsize+1)/1000;
taugrid=(Mind-M-1)*Ts;

figure(5); mesh(taugrid, nugrid, abs(ambiguity))
 
 xlabel('tau (\mu sec)')
 ylabel('\nu (kHz)')
 

figure(6); contour(taugrid, nugrid, abs(ambiguity))
 
 xlabel('tau (\mu sec)')
 ylabel('\nu (kHz)')


%%% Creates a rectangular pulse of width tau from 0<t<tau %%%
function p = rpulse(t,tau)
    p = (t<=tau)&(t>=0);
end

