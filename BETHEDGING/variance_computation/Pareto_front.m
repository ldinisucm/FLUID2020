%Computes Pareto front from theoretical expressions for mean and variance

%Parameter definitions
set='a'


kappa_pi_set='fig2b'

switch kappa_pi_set
case {'fig2c'}
        kappa1 = 10; kappa2 = 10;
        %pi1=0.4;pi2=0.4;
case {'fig2b'}
  kappa1=1.0;kappa2=1.0;
  %pi1=0.24;
  %pi2=pi1;
case {'fig2a'}
  kappa1=0.1;kappa2=0.1;
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

%kappa_set=[kappa1 kappa2];
%ks=[kA1 kB1 kA2 kB2];

Lambda=[];
Var=[];
PI0=[];
options=optimset();
%for alpha=1:-.02:.94
pis_guess=[kappa1 kappa2]; %initial guess
for alpha=1:-.1:.3
	alpha
	extra_param=[kappa1 kappa2 kA1  kB1 kA2 kB2 alpha];
	disp('Running fminsearch')
  	optimal_pis = fminsearch('fun_mean_variance_for_Pareto',pis_guess,options,extra_param)
	PI0=[PI0;optimal_pis];
	disp('Computing Lambda and VarLambda for optimal pi''s')
	[V,L] = fun_variance_computation([kappa1 kappa2 optimal_pis(1) optimal_pis(2)],[ kA1  kB1 kA2 kB2]);
	Lambda=[Lambda L];
	Var=[Var V];
	pis_guess=optimal_pis;

end
plot(Lambda,Var,'o-')
xlabel('Mean g.r.')
ylabel('Variance')

