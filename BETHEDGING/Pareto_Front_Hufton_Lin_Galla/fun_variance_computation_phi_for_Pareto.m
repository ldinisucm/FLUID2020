function [objective]=fun_variance_computation_phi_for_Pareto(pis, extra_param)

%% Paramètres
% kappa1, kappa2, kA1.. sont les paramètres de la simulation

%% extract parameters from input
kappa1=extra_param(1);
kappa2=extra_param(2);
pi1=pis(1);
pi2=pis(2);
kA1=extra_param(3);
kB1=extra_param(4);
kA2=extra_param(5);
kB2=extra_param(6);
alpha=extra_param(7);
numpuntos=extra_param(8);

%Compute theoretical distribution according to Jérémie's notes
%compute phi plus minus
%Delta1=kA1-kA2; %wrong, but anyway kA2=0,kB1=0 then nothing changed.
%Delta2=kB1-kB2;
Delta1=kA1-kB1; %Correct
Delta2=kA2-kB2;
phiv1plus=(Delta1-pi1-pi2+sqrt((Delta1-pi1-pi2)^2+4*pi2*Delta1))/(2*Delta1);
phiv1minus=(Delta1-pi1-pi2-sqrt((Delta1-pi1-pi2)^2+4*pi2*Delta1))/(2*Delta1);
phiv2plus=(abs(Delta2)+pi1+pi2+sqrt((abs(Delta2)+pi1+pi2)^2-4*pi2*abs(Delta2)))/(2*abs(Delta2));
phiv2minus=(abs(Delta2)+pi1+pi2-sqrt((abs(Delta2)+pi1+pi2)^2-4*pi2*abs(Delta2)))/(2*abs(Delta2));

phi1plus=max(phiv1plus,phiv1minus);
phi2plus=min(phiv2plus,phiv2minus);

phi1minus=min(phiv1plus,phiv1minus);
phi2minus=max(phiv2plus,phiv2minus);

g=kappa1/Delta1/(phi1plus-phi1minus);
h=-kappa2/Delta2/(phi2minus-phi2plus);
Q1=kappa2/(kappa1+kappa2);
Q2=kappa1/(kappa1+kappa2);
%Normalization of P1, definition of p11plus, the regular part
%numpuntos=200; %number of points for integration
phi=linspace(phi2plus,phi1plus,numpuntos);
Dphi=phi(2)-phi(1);
integrand1=(phi1plus-phi).^(g-1).*(phi-phi2plus).^h.*(phi-phi1minus).^(-g-1).*(phi2minus-phi).^(-h);
int1=Dphi*sum(integrand1);
%int1=trapz(Dphi*integrand1);
N1=Delta1*Q1/int1;
P1=N1*integrand1/Delta1;
Norm_P1=Dphi*sum(P1);

%Normalization of P2, definition of p22plus the regular part
integrand2=(phi-phi2plus).^(h-1).*(phi1plus-phi).^g.*(phi-phi1minus).^(-g).*(phi2minus-phi).^(-h-1);
int2=Dphi*sum(integrand2);
N2=abs(Delta2)*Q2/int2;
P2=N2*integrand2/abs(Delta2);
Norm_P2=Dphi*sum(P2);
Norm_P=Norm_P1+Norm_P2;

%Average growth computation
%first part

lambda1=sum(P1.*(kA1*phi+kB1*(1-phi)))*Dphi;

%second part
lambda2=sum(P2.*(kA2*phi+kB2*(1-phi)))*Dphi;

lambda=lambda1+lambda2;

%Variance computation


v1=Delta1*phi.*(1-phi)-pi1*phi+pi2*(1-phi);
q=(v1.*P1).^2./(kappa1*P1+kappa2*P2);

f1=(kA1*(phi)+kB1*(1-phi)-lambda);
norm_f1=sum(P1.*f1)*Dphi;
  f2=(kA2*(phi)+kB2*(1-phi)-lambda);
norm_f2=sum(P2.*f2)*Dphi;
norm_f1f2=norm_f1+norm_f2;
I=cumsum(P1.*f1+P2.*f2)*Dphi;
integrand=I.^2./q;

varlambda_theor=Dphi*sum(integrand(2:end-1));

integrand2=f1.^2.*P1+f2.^2.*P2;
lambda_theor_J=lambda;

%output values
Var_theor=varlambda_theor;
lambda_theor=lambda_theor_J;

objective=-alpha*lambda_theor+(1-alpha)*Var_theor;

