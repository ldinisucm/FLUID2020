%cargar los datos experimentales
tiempo=5;%12 24 48 72 96
[omega_exp,Gp_exp,Gpp_exp]=datos_exp_fun_minimizar('Static',tiempo)
%Comprobamos si hay algún dato negativo y lo quitamos
Gp_prune=0;
wp=0;
Gpp_prune=0;
wpp=0;
jp=1;
jpp=1;
display('Quitando valores negativos')
for i=1:length(Gp_exp)
	if Gp_exp(i)>=0
		Gp_prune(jp)=Gp_exp(i);
		wp(jp)=omega_exp(i);
		jp=jp+1;
	end
end
for i=1:length(Gpp_exp)
        if Gpp_exp(i)>=0
		Gpp_prune(jpp)=Gpp_exp(i);
		wpp(jpp)=omega_exp(i);
		jpp=jpp+1;
	end
end
Gp_exp=Gp_prune'
omega_exp_Gp=wp'
Gpp_exp=Gpp_prune'
omega_exp_Gpp=wpp'



%Guess inicial de los parámetros
%Buenos para 72h
eta=0.6;
x=1.1;
fc=6
Gp0=0.00002;
Gpp0=.24;
F=1.2;
%Buenos para 12h
eta=1;
x=1.1;
fc=60
Gp0=0.0002;
Gpp0=.4;
F=0.04;
%Buenos para 24h
eta=1;
x=1.05;
fc=15
Gp0=0.00002;
Gpp0=.4;
F=0.1;
%Buenos para 48h
eta=1;
x=1.08;
fc=35
Gp0=0.00002;
Gpp0=.4;
F=0.5;
%Buenos para 96h
eta=0.09;
x=1.08;
fc=1.2
Gp0=0.00002;
Gpp0=.3;
F=5.5;
parametros=[x fc Gp0  eta Gpp0 F];

tolerancia=1e-8
iter=500000
options=optimset('TolFun',tolerancia,'TolX',tolerancia,'MaxFunEvals',iter,'MaxIter',iter);
%options=optimset;
parmFit=fminsearch('fun_minimizar_3',parametros,options,omega_exp_Gp,omega_exp_Gpp,Gp_exp,Gpp_exp)

x=parmFit(1);
fc=parmFit(2);
wc=2*pi*fc;
Gp0=parmFit(3);
eta=parmFit(4);
Gpp0=parmFit(5);
F=parmFit(6);

omega=[0.01 0.1 1:10 10:5:100 100:10:200];
X=omega./sinh(pi*omega./wc);
I0=besseli(0,F*X);
Gp_teor=Gp0*omega.^(x-1).*I0;

X=omega./sinh(pi*omega./wc);
I0=besseli(0,F*X);
Gpp_teor=Gpp0*(eta*omega+Gp0*omega.^(x-1).*I0);

loglog(omega,Gp_teor,omega,Gpp_teor)
