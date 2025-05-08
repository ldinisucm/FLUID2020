%cargar los datos experimentales
tiempo=5;%12 24 48 72 96
[fexp,gpexp,gppexp]=datos_exp_fun('Static',tiempo);



x=1.3
tau0=1
N=100000;
omega=2*pi*10.^[-2:.01:2];


eta0=.0020
G0=8e-3
fc=1;
wc=2*pi*fc;
tiempo=10;
for i=1:length(omega)
  w=omega(i);
  tau=linspace(tau0,100000*tau0,N);
  Dtau=(max(tau)-min(tau))/N;
  gp(i)=G0*Dtau*sum(tau.^(-x)/(1+(w/wc)^2).*real(j*w*tau./(1+j*w*tau)));
  gpp(i)=G0*Dtau*sum(tau.^(-x)/(1+(w/wc)^2).*imag(j*w*tau./(1+j*w*tau)))+eta0*(w);
%  gp(i)=G0*Dtau*sum(tau.^(-x)/(1+tiempo*heaviside(w-wc)*(w-wc)^2).*real(j*w*tau./(1+j*w*tau)));
%  gpp(i)=G0*Dtau*sum(tau.^(-x)/(1+tiempo*heaviside(w-wc)*(w-wc)^2).*imag(j*w*tau./(1+j*w*tau)))+eta0*(w);
%  gp(i)=G0*Dtau*sum(tau.^(-x)*exp(-w^2/wc^2).*real(j*w*tau./(1+j*w*tau)));
%  gpp(i)=G0*Dtau*sum(tau.^(-x)*exp(-w^2/wc^2).*imag(j*w*tau./(1+j*w*tau)))+eta0*(w);
end
f=omega/2/pi;
loglog(f,gp,'r')
loglog(f,gpp,'r--')
legend('G''','G''''',"G' theor","G'' theor")
%figure
%loglog(fexp,gpexp./gppexp,'o',f,gp./gpp)
%figura formateada para artículo
figure
loglog(f,gp,'r')
hold on
loglog(f,gpp,'r--')
loglog(fexp,gpexp,'ro','MarkerFaceColor','r')
loglog(fexp,gppexp,'ro','MarkerFaceColor','w')
xlabel('frequency_{}(s^{-1})')
%ylabel('G','Rotation',0)
ylabel("G',G''_{}(Nm^{-2})")
legend("G' theor","G'' theor","G' exp","G'' exp",'Location','northwest')
text(0.02,0.5e-5,"96h",'FontSize',14,'FontWeight','bold')

