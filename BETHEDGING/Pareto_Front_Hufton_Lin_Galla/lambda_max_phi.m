%Computes "Kelly's" like point

%Parameter definitions
set='a'


kappa_pi_set='fig2c'

switch kappa_pi_set
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
        %pi1=0.4;pi2=0.4;
	    A=0.5
        beta0=10
        beta_final=200000
pis_guess=[1 1] ; %initial guess
n_iterations=200;

case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
    A=0.5
  beta0=10
  beta_final=1000
pis_guess=[kappa1 kappa2];
n_iterations=100;
  %pi1=0.24;
  %pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
  A=0.5
   beta0=10
  beta_final=1000
pis_guess=[kappa1 kappa2]; %initial guess
n_iterations=50;

  %pi1=0.064;pi2=0.064;
case {'DL'}
kappa1 = 0.1; kappa2 = 0.33;
%pi1 = kappa1; pi2 = kappa2;
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

%kappa_set=[kappa1 kappa2];
%ks=[kA1 kB1 kA2 kB2];

Lambda=[];
Var=[];
PI0=[];
options=optimset();
%pis_guess=[0.3841    0.2777]
	extra_param=[kappa1 kappa2 kA1  kB1 kA2 kB2 numpuntos_ave];
	disp('Running fminsearch')
	tic
  	%[optimal_pis,funcion_objetivo] = fminsearch('fun_lambda_for_minimization',pis_guess,options,extra_param)
%  	[optimal_pis,funcion_objetivo] = simulated_annealing('fun_lambda_for_minimization',pis_guess,extra_param)
  	[optimal_pis,funcion_objetivo] = simulated_annealing_var('fun_lambda_for_minimization_phi',A,beta0,beta_final,n_iterations,pis_guess,extra_param)
	toc
	PI0=[PI0;optimal_pis];
	disp('Computing Lambda for optimal pi''s')
	L=fun_lambda_phi([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2],numpuntos_ave);

	Lambda=[Lambda L];

