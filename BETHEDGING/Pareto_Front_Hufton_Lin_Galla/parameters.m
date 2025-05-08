%Parameter definitions
set='a'


kappa_pi_set='DLc'

switch kappa_pi_set
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
	A=0.5
	beta0=10
	beta_final=5000
pis_guess=[1 1]; %initial guess
n_iterations=50;
        %pi1=0.4;pi2=0.4;
case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
  A=0.5
  beta0=10
  beta_final=10000
  %beta_final=1000
pis_guess=[kappa1 kappa2]; %initial guess
n_iterations=60;
  %pi1=0.24;
  %pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
  %pi1=0.064;pi2=0.064;
  A=0.5
  beta0=10
  beta_final=1000
pis_guess=[kappa1 kappa2]; %initial guess
n_iterations=50;
case {'DLa'}
kappa1 = 0.01; kappa2 = 0.033;
%pi1 = kappa1; pi2 = kappa2;
A=0.005
  beta0=5
  beta_final=2000
pis_guess=[kappa1 kappa2]; %initial guess
n_iterations=100;
case {'DLb'}
kappa1 = 0.1; kappa2 = 0.33;
%pi1 = kappa1; pi2 = kappa2;
A=0.05
  beta0=5
  beta_final=2000
pis_guess=[kappa1 kappa2]; %initial guess
n_iterations=100;
case {'DLc'}
kappa1=1; kappa2 = 3.3;  
A=0.5
  beta0=10
  beta_final=1000
pis_guess=[kappa1 kappa2]; %initial guess
n_iterations=50;
end

switch set
case {'a'}
%set a) Hufton
kA1 = 2; kB2 = -0.2;
kA2=-2; kB1=0.2;
case {'b'}
%set b) Hufton
kA1 = 0.5;
kB1=0.0001;
kA2=0.0001;
kB2=0.3250;
otherwise
        kA1=2;
        kB2=0.2;
        kA2=0;
        kB1=0;
end
numpuntos_ave=1e6;
numpuntos_ave=5e5;
numpuntos_var=1000;

