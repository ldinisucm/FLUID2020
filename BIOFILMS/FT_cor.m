tau=1;
b=4
M=2^24;
N=2*M;
deltaM=b/M;
tM=[0:deltaM:b];
h=exp(-(tM/tau));
figure
plot(tM,h)
xlabel('t')
h=[h zeros(1,N-M-1)];
fn=1/deltaM*[0:N/2-1]/N;
%t=[-t(end:-1:1) t];
theta=2*pi*[0:N/2-1]/N;
alpha0=-(1-cos(theta))./theta.^2+j*(theta-sin(theta))./theta.^2;
W=2*(1-cos(theta))./theta.^2;
dft_aux=fft(h);
dft=dft_aux(1:N/2);
I=W.*dft+alpha0*h(1)+exp(j*2*pi*fn*b).*conj(alpha0);
I=deltaM*I;
%Obtener las frecuencias wM
factor=N/M;
I=I(1:factor:end);
fn=fn(1:factor:end);
dft=dft(1:factor:end);
%teor
FTteor=tau./(1+j*2*pi*fn*tau);
figure
loglog(fn,abs(I),fn,abs(FTteor))
title('abs')
figure
loglog(fn,real(I),fn,real(FTteor))
title('real')
legend('I','Teor')
figure
plot(fn,deltaM*real(dft),fn,real(I),fn,real(FTteor))
title('fft vs I vs Teor. Real')
figure
plot(fn,deltaM*imag(dft),fn,imag(I),fn,imag(FTteor))
title('fft vs I vs Teor. Im')
figure
loglog(fn,0.5*(deltaM*real(dft)+real(I)),fn,real(FTteor))
title('truqui')
