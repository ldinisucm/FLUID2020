
odds=[2.3 2]';
odds=[1/(1-.4) 1/0.4]';

b1_vec=[0.54:.01:0.96];
Pext=zeros(size(b1_vec));
W_theor=zeros(size(b1_vec));
delta_0_vec=zeros(size(b1_vec));
delta_1_vec=zeros(size(b1_vec));
for k=1:length(b1_vec)
	%Initialization
	%b=P;
	b=[b1_vec(k) 1-b1_vec(k)]'
	logC0=0;
	tray=[0;1]
	num_pasos=16
	p=.2;
	P=[1-p;p] %prob 0, prob 1
	delta_0=log(b(1)*odds(1));
	delta_1=log(b(2)*odds(2));
	delta_0_vec(k)=delta_0;
	delta_1_vec(k)=delta_1;
	
	
	%delta_plus=0.6;
	%delta_minus=-0.6;
	threshold=-0.1;
	if logC0+delta_0<threshold
		P(1)=0;
	end
	if logC0+delta_1<threshold
		P(2)=0;
	end
	P
	
	for i=2:num_pasos
		aux=[tray;tray];
		col_nueva=[zeros(length(tray),1);ones(length(tray),1)];
		tray=[col_nueva aux];
		[f,c]=size(tray);
		Pold=P;
		P=zeros(f,1);
		for j=1:f
			%tray(j,:)	
			logC=logC0+sum(tray(j,:))*delta_1+(c-sum(tray(j,:)))*delta_0;
	 		if logC>threshold
				index=dec(tray(j,:));
				previous_index=dec(tray(j,1:i-1));
				P(index)=(p*tray(j,i)+(1-p)*(1-tray(j,i)))*Pold(previous_index);
				P(index);
			end
		end
		%P
	end
	Pext(k)=1-sum(P)
	W_theor(k)=[1-p p]*log(b.*odds)
end
figure
plot(W_theor,Pext,'-o')
%plot(W_theor,Pext)
xlabel('\langleW\rangle')
ylabel('P_{ext}')
figure
b=1-b1_vec;
%plot(b,W_theor,b,delta_0_vec,b,delta_1_vec)
plot(b,W_theor)
xlabel('b')
ylabel('W_theor')
figure
plot(W_theor,delta_0_vec,W_theor,delta_1_vec)
%figure with trade-off branch
Wtoff=W_theor(b1_vec<1-p & W_theor>=0 )
Ptoff=Pext(b1_vec<1-p & W_theor>=0 )
figure
plot(W_theor,Pext,'r')
hold on
plot(Wtoff,Ptoff,'b','LineWidth',2)
xlabel('\langleW\rangle')
ylabel('P_{ext}')

