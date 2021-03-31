%This code is created for use in the MMath Project: Investigating the
%potetntial of early warning signals in disease elimination. 

%Gillespie simulation of the SIS model 

%{
Inputs: 
I0 - initial number of infected individuals 
N - total population size 
beta - infectivity constant 
gamma - recovery rate 
max_time - time for which the simulation runs 
dt - timestep for incidence aggregation
plots - if 0 then does not plot any graphs, if 1 plots output graphs 

Outputs:
Timeseries for Susceptibles, Infected and Incidence aggregated over a
timestep of length dt 
%}

function [S,I,T,incidence] = gillespie_SIS(I0,N,beta,gamma,max_time,dt,plots)

clear S I T incidence

%Starting values for S - susceptibles, I - infected 
S(1)=N - I0;
I(1)=I0; 

%Initialising time (t) and time vector (T). Index is the position of time t
%in the vector T
T(1)=0;
t=0;
ind=1;

%Time values at which incidence is calculated 
L=linspace(dt,max_time,max_time/dt);
incidence = zeros(1,max_time/dt);

%Main iteration 
while t<max_time + dt
    
    %Calculating the total rate of reactions 
    rate(1) = beta*S(ind)*I(ind)/(N); %Rate of infections 
    rate(2) = gamma*I(ind); %Rate of recoveries 
    lambda = sum(rate);
    
    if lambda>0 % i.e. there is still happening on in the simulation
        
        %Calculating the timestep 
        r1=rand;
        DT = -log(r1)/lambda;
        t = t+DT;
        T(ind+1)=t;
        
        %Calculating the next event 
        r2 = rand*lambda;
        event = find(cumsum(rate)>=r2,1);
        if event == 1 %the next event is an infection
            S(ind+1)=S(ind)-1;
            I(ind+1)=I(ind)+1;  
            inc_position = max(1,min(ceil(t/dt),length(incidence)));
            incidence(inc_position) = incidence(inc_position)+1; %recording the infection in the incidence log 
        else 
            if event == 2 % the next event is a recovery
            S(ind+1)=S(ind)+1;
            I(ind+1)=I(ind)-1;
            end
        end
        
    else % i.e. there is nothing going on in the simulation, the disease has died out
        t=max_time + dt;
        T(ind+1)=t;
        S(ind+1)=S(ind);
        I(ind+1)=I(ind);
    end
    
    ind=ind+1;
end

incidence(length(incidence))=incidence(length(incidence))/2; 
%The final value for incidence collects double the data 

if plots==1
figure('Name','Susceptible and Infected Timeseries')
plot(T,S,'b','DisplayName','Susceptible')
hold on 
plot(T,I,'r','DisplayName','Infected')
hold off
legend
figure ('Name','Incidence Timeseries')
plot(L,incidence,'DisplayName','Incidence')
legend
end

end