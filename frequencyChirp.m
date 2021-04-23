fc = 10.5e9;        % center frequency = 10.5 GHz
lambda = c/fc;      % wavelength of radar system
BW = 125e6;         % total system bandwidth = 125 MHz
fp = 1e3;           % pulse repetition frequency = 1 kHz
Tp = 1/fp;          % pulse repetition interval = 1 ms
Np = 64;            % number of pulses


beta = 120e6;       % sweep bandwidth = 120 MHz
tau = 80e-6;        % pulse width = 80 usec

stfu = fc-beta/2;
edfu = fc+beta/2;

stfd = fc+beta/2;
edfd = fc-beta/2;

ru = linspace(stfu,edfu,1000);
rd = linspace(stfd,edfd,1000);

fc = linspace(fc,fc,1000)
t = linspace(0,tau,1000);

figure(1)
plot(t/1e-6,ru/1e9,'r-',t/1e-6,rd/1e9,'b-',t/1e-6,fc/1e9,'--k','LineWidth',2)
xlabel('time (us)')
ylabel('frequency (GHz)')
legend('Up Ramp','Down Ramp','Carrier Frequency')
set(gca,'fontsize',18)