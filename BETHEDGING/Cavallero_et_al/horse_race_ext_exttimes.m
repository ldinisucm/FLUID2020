%Horse race with correlations
%november 19
%Simulates horse races outcomes with probabilities conditioned on last outcome

%%%%%%%%%Definition of probabilities

%Conditional probabilities P(i,j)=P(i|j) prob. of horse i winning provided horse j won in the previous round
%P=[p(1|1) p(1|2) p(1|3)
%   p(2|1) P(2|2) p(2|3)
%   p(3|1) P(3|2) p(3|3)];
   
P=[0.2  .2 .3
   .4   .5 .4
   .4   .3 .3];
P=[0.4  .4 
   .6   .6];

q=length(P); %Number of horses   
%Probabilities for first race. They shouldn't really matter for the long run
px0=[.4 
     .6];
%%%%%%%%%Definition of odds
%odds you get for winning
odds=zeros(q,1);
odds=[4 4 2]';
odds=[2.3 1.9 ]';
odds=[2.3 2 ]';
%%%%%%%%%Definition of bets
%Conditional bets b(i,j)=b(i|j). How much to bet on horse i provided horse j won in the previous round
b=zeros(size(P));
%Some winning but not optimal betting strategy
b=[.25 .2 .2
   .5  .35 .6
   .25 .45 .2 
];

b=P; %Uncomment this line to simulate Kelly's case
%b=1./odds.*ones(q); %Uncomment this one to play the null strategy
%b=1./odds.*ones(q)+0.2*rand(q); %Uncomment this one to play something CLOSE to the null strategy
%norm=sum(b);
%b=b./norm; %normalized
%Uncorrelated
b=b(:,1).*ones(size(b))


%%%%Rest of parameters and variables
num_races=100; %number of races   

umbrales=[-0.05 -.4];
num_umbrales=length(umbrales);
p_ext=zeros(size(umbrales));
pext_th=p_ext;
for k=1:num_umbrales
	%%The realization loop
	realizations=4000;
	extintions=0;
	umbral=umbrales(k);
	t_ext=zeros(1,realizations);
	for l=1:realizations
	%initialize variables
	S=0; %initial (log) capital
	log_capital=0; %Stores the logarithm of the capital for 1 race
	horse=zeros(1,num_races); %Stores which horse won last race. 
	horse(1)=1; %First race was won by 1. Does not matter since 1 race probabilities are the same no matter what column you take;
	%In matlab, first element of a vector is 1, not 0
	
	%%%Start with first race probabilities.
	p=zeros(q,1); %Stores the probabilities of the q horses for each race. The columno corresponding to the horse that won before;
	p=px0; %first race.
	
	%%The races loop
	
	for i=2:num_races+1
	     r=rand; %random number from 0 to 1;
	     %Who wins this race with that r?
	     cumulative_p=cumsum(p); %Cumulative sum. Every element is the sum of the previous p values [p1 p1+p2 p1+p2+p3, etc....]
	     check=(r<=cumulative_p); %check is a vector with elements 0 or 1 depending whether the condition is met or not. The first appearance of 1 is the horse that wins
	     [m,where_is_the_first_one]=max(check); 
	     %max returns in the second variable (where_is_the_first_one) the index where the maximum (1) of the vector is met. Matlab function "max" is such that if the maximum occurs in several places, the first index is returned. What a lucky coincidence!
	     horse(i)=where_is_the_first_one; %horse(i) wins this race
	     %Compute the capital
	     log_capital=log(odds(horse(i))*b(horse(i),horse(i-1)));
	     %I am assuming bets for horse=1 in the first race. Shouldn't change much allowing for different bets in the first race.
	     S(i)=log_capital;
	     p=P(:,horse(i));%Defines probabilities for next race 
	end 
	accrued_capital=cumsum(S); %The cumulative capital
	t_ext(l)=0;
	test=accrued_capital<umbral;
	if any(test)
		extintions=extintions+1;
		[maximo,t_ext(l)]=max(int8(test));
	end
	
	end %realizations
	t_ext=nonzeros(t_ext)-1; %we take into account that initial time is step 1, not 0 
	figure
	histogram(t_ext,'Normalization','pdf')
	p_ext(k)=extintions/realizations
	
	
	%%%Stationary solution
	%compute stationary probability for correlated games
	sufficiently_big_number=1000; %Number of races to compute stationary prob.
	pi_st=1/q*ones(q,1); %initial guess prob. Shouldn't really matter.
	for i=1:sufficiently_big_number
	  pi_st=P*pi_st;
	end
	%Expected returns for the games played
	W_theory=sum(P.*log(odds.*b)*pi_st);
	pcomp=P(:,1)
	bcomp=b(:,1)
	W_th_uncorr=sum(pcomp.*log(odds.*bcomp))
	Var=sum(pcomp.*log(odds.*bcomp).^2)-W_th_uncorr^2
	
	hold on
	Dt=0.001;
	t=[0:Dt:max(t_ext)+1];
	x0=abs(umbral);
	sigma=sqrt(Var);
	mu=W_th_uncorr;
	g=x0/sigma/sqrt(2*pi)./t.^(3/2).*exp(-(x0+mu*t).^2./(2*sigma^2*t));
	%plot(t,g,'r')
	%Norma_g=Dt*sum(g(2:end))
	Norma_g=exp(-4*x0*mu/2/sigma^2) %expresion teórica
	plot(t,g/Norma_g)
	xlabel('extinction run')
	legend('simulations','theory')
	title(['thres.= ' num2str(-x0)])
	tfinal=max(t);
	F=.5*(1+erf((x0+mu*tfinal)/sqrt(2*sigma^2*tfinal))-exp(-2*x0*mu/sigma^2)*(1+erf((mu*tfinal-x0)/sqrt(2*sigma^2*tfinal))));
	pext_th(k)=1-F;
	p_ext(k)
        pext_th(k)	

	%%%Kelly's expected gains
	%Compute Kelly's solution.
	b_Kelly=P; %Proportional betting
	W_Kelly=sum(P.*log(odds.*b_Kelly)*pi_st);
end %%umbrales
figure
plot(umbrales,p_ext,'-o',umbrales,pext_th)
legend('simulation','theory')
xlabel('threshold')
ylabel('Prob. ext')

%Check horse stationary probabilities by computing the histogram of horse wins
show_histogram=0; %change to 0 if you don't want the histogram to be displayed every time
if show_histogram
	figure
	bar(1:q,pi_st,.25)
	hold on
	h=histogram(horse,'Normalization','probability','BinWidth',0.2);
	xlabel('horse')
	ylabel('probability')
	legend('St prob','observed freq')
	disp('----Checking stationary distribution----')
	disp('Stationary distribution')
	pi_st
	disp('Histogram results')
	where_is_the_relevant_data=h.Values>0; %Just a trick to display only the relevant values
	indices=[];
	for i=1:length(where_is_the_relevant_data)
		if where_is_the_relevant_data(i)
		  indices=[indices i];
		end
        end
	h.Values(indices)'
end
%%%%Results%%%%%

disp('---Expected Capital growths---')
W_Kelly
W_theory

%%Plotting stuff
%figure
%plot(horse) %histogram is more interesting in fact

figure
race=1:num_races+1; %race vector for plotting
expected_W=W_theory*race;%The expected average W for the actual b played
Kelly_W=W_Kelly*race; %The expected average W for Kelly's prop. betting b=P
plot(race,accrued_capital,race,expected_W,race,Kelly_W)
legend('Sim. capital','Expected W','Kelly''s')
xlabel('race number')
ylabel('log-Capital')

%figure
%plot(race,exp(accrued_capital),race,exp(expected_W))

