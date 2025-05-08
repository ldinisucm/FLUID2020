%cargar los datos experimentales
tiempo=4;%12 24 48 72 96
[fexp,gpexp,gppexp]=datos_exp_fun('Static',tiempo);



x=1.1
tau0=.1
N=1e5;
omega=2*pi*10.^[-2:.01:2];


eta0=.020
G0=6e-1
fc=6;
wc=2*pi*fc;
for i=1:length(omega)
  w=omega(i);
  tau=linspace(tau0,100000*tau0,N);
  Dtau=(max(tau)-min(tau))/N;
  gp(i)=G0*Dtau*sum(tau.^(-x)/(1+(w/wc)^4).*real(j*w*tau./(1+j*w*tau)));
  gpp(i)=G0*Dtau*sum(tau.^(-x)/(1+(w/wc)^4).*imag(j*w*tau./(1+j*w*tau)))+eta0*(w);
%  gp(i)=G0*Dtau*sum(tau.^(-x)/(1+tiempo*heaviside(w-wc)*(w-wc)^2).*real(j*w*tau./(1+j*w*tau)));
%  gpp(i)=G0*Dtau*sum(tau.^(-x)/(1+tiempo*heaviside(w-wc)*(w-wc)^2).*imag(j*w*tau./(1+j*w*tau)))+eta0*(w);
%  gp(i)=G0*Dtau*sum(tau.^(-x).*real(j*w*tau./(1+j*w*tau)));
%  gpp(i)=G0*Dtau*sum(tau.^(-x).*imag(j*w*tau./(1+j*w*tau)));
end
f=omega/2/pi;
loglog(f,gp)
loglog(f,gpp)
legend('G''','G''''',"G' theor","G'' theor")
%figure
%loglog(fexp,gpexp./gppexp,'o',f,gp./gpp)
