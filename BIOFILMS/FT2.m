ts=.001;
tau=1;
t=[0:ts:200*tau-ts];
%t=[-t(end:-1:1) t];
y=exp(-abs(t/tau));
figure
plot(t,y)
xlabel('t')
fs=1/ts;
N=length(t);
f=(0:N-1)/N*fs;
ygorro=ts*fft(y);
%ygorro=ts*(fft(fftshift(y)));
figure
loglog(f,abs(ygorro))
title('abs')
%FTteor=tau./(1+j*2*pi*f*tau)+tau./(1-j*2*pi*f*tau);
FTteor=tau./(1+j*2*pi*f*tau);
hold on
loglog(f,abs(FTteor))
figure
loglog(f,real(ygorro),f,real(FTteor))
title('real')

