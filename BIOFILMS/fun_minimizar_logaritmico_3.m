function dist=fun_minimizar_logaritmico(parametros,omega_exp_Gp,omega_exp_Gpp,Gp_exp,Gpp_exp)
%Calcula la distancia cuadrática entre las funciones y es lo que tratamos de minimizar
x=parametros(1);
fc=parametros(2);
wc=2*pi*fc;
Gp0=parametros(3);
eta=parametros(4);
Gpp0=parametros(5);
F=parametros(6);

X=omega_exp_Gp./sinh(pi*omega_exp_Gp./wc);
I0=besseli(0,F*X);
Gp_teor=Gp0*omega_exp_Gp.^(x-1).*I0.^(-x);
Gp_teor=log(Gp_teor)
%tam_Gp_teor=size(Gp_teor)

X=omega_exp_Gpp./sinh(pi*omega_exp_Gpp./wc);
I0=besseli(0,F*X);
Gpp_teor=Gpp0*(eta*omega_exp_Gpp+Gp0*omega_exp_Gpp.^(x-1).*I0.^(-x));
Gpp_teor=log(Gpp_teor);
%tam_Gpp_teor=size(Gpp_teor)
pesoGp=1;
pesoGpp=1;
sumpesos=pesoGp+pesoGpp;
dist=pesoGp/sumpesos*sum((Gp_exp-Gp_teor).^2)+pesoGpp/sumpesos*sum((Gpp_exp-Gpp_teor).^2)
