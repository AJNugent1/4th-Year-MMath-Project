%This code is created for use in the MMath Project: Investigating the
%potetntial of early warning signals in disease elimination. 

%Here the error in approximations against analytic results is calculated
%for an increasing number of realisations. 

clear 

%Simulation Parametes 
Repetitions = 200; %number of simulations
N=10000; %population size
beta = 0.9; %infectivity constant
gamma = 1; %recovery rate
max_time = 10; %length of each simulation
delta=0.1; %incidence aggregation timestep 

%Equation Free Method Parameters 
M = N/2;
sections=50;
sectionwidth=M/sections; 
topf=zeros(1,sections);
topd=zeros(1,sections);
bot=zeros(1,sections);
L=linspace(delta,max_time,max_time/delta);

%Analytic Results 
LINSPACE = linspace(0,M,sections);
for y=1:sections
    x = LINSPACE(y)/N;
    F_x(y) = -N*(gamma*x - beta*x*(1-x)) ;
end
potential_ana = - cumtrapz(LINSPACE,F_x);

Error_potential_L2 = zeros(1,Repetitions);
Error_drift_L2 = zeros(1,Repetitions);
Error_potential_sup = zeros(1,Repetitions);
Error_drift_sup = zeros(1,Repetitions);

for R=1:Repetitions
    clear F_est potential_est 
    %Note that topf and bot are NOT cleared, as they will be cumulative
    %over the repetitions to give more and more accurate results. 
    
    if rem(R,Repetitions/10)==0
        R
    end

%Running each simulation 
clear S I T incidence
startplace=N*0.5*rand; %random starting point with more than half of the population susceptible
[S,I,T,incidence]=gillespie_SIS(startplace,N,beta,gamma,max_time,delta,0); 

Lin = linspace(0,max_time,max_time*1000);
I_Lin = interp1(T,I,Lin);

%The equation free method 
for t=1:length(I_Lin) - 1
    for s=1:sections
        if I_Lin(t)>=(s-1)*sectionwidth && I_Lin(t)<=s*sectionwidth %determining which section Z(t) lies in 
            topf(s)=topf(s) + ( I_Lin(t+1) - I_Lin(t) )/( 0.001 ); %dI/dT 
            bot(s)=bot(s)+1;
        end
    end
end

for s=1:sections
    if bot(s)>0 
    F_est(s)=topf(s)/bot(s);
    else if s>1
        F_est(s)=F_est(s-1);
        else 
            F_est(s)=0;
        end
    end   
end

LM = linspace(0,M,sections);
potential_est = - cumtrapz(LM,F_est);

Error_potential_L2(R) = norm(abs(potential_est - potential_ana));
Error_drift_L2(R) = norm(abs(F_est - F_x));

Error_potential_sup(R) = max(abs(potential_est - potential_ana));
Error_drift_sup(R) = max(abs(F_est - F_x));
end

figure (1)
subplot(2,1,1)
plot([1:1:Repetitions],Error_potential_L2,[1:1:Repetitions],Error_potential_sup)
legend('Supremum Norm','L^2 Norm')
title('Erros in Potential Function')
subplot(2,1,2)
plot([1:1:Repetitions],Error_drift_L2,[1:1:Repetitions],Error_drift_sup)
legend('Supremum Norm','L^2 Norm')
title('Erros in Drift Function')