%cargar los datos experimentales
tiempo=5;%12 24 48 72 96
datos_exp_fun('Static',tiempo)



x=1.3
tau0=1
N=100000;
omega=2*pi*10.^[-2:.01:2];


eta0=.0018
G0=6e-6
F=.45;
fc=10;
wc=2*pi*fc;
for i=1:length(omega)
  w=omega(i);
  tau=linspace(tau0,100000*tau0,N);
  Dtau=(max(tau)-min(tau))/N;
  X=w/sinh(pi*w/wc);
  gp(i)=G0*Dtau*sum(tau.^(-x)*besseli(0,F*X).*real(j*w*tau./(1+j*w*tau)));
  gpp(i)=G0*Dtau*sum(tau.^(-x)*besseli(0,F*X).*imag(j*w*tau./(1+j*w*tau)))+eta0*(w);
end
f=omega/2/pi;
loglog(f,gp)
loglog(f,gpp)
