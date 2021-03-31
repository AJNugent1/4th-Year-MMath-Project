%This code is created for use in the MMath Project: Investigating the
%potetntial of early warning signals in disease elimination. 

%Here the drift function and potential surface for prevalence data is
%estimated from simulated data. 

clear 


%%%%%%%%%%%%%%%%%%%%%%%%%%% Model setup %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Parameters for the gillespie simulation 
Repetitions = 100; %number of simulations
N=10000; %population size
beta0 = 1.3; %infectivity constant
beta1 = 0.9;
gamma = 1; %recovery rate
max_time = 10; %length of each simulation
delta=0.1; %incidence aggregation timestep 

%The process is repeated with slowly varing beta values
nb = 5; %the number of different beta values 

poly_fit_data = zeros(nb,1);
Beta_vec = zeros(nb,1);

for b = 1:nb % b controls the value of beta
    
clear M topf topd F_est D_est potential_est

if nb>1 %If we are testing multiple beta values use this formula to find beta
beta = beta0 - (beta0-beta1)*((b-1)/(nb-1))
else 
    beta = beta0 %If we are testing only one beta value, beta=beta0
end

Beta_vec(b) = beta; %Storing the values of beta for plotting later

%Resetting the equation free method for each beta value 
M = N/2;
sections=50;
sectionwidth=M/sections; 
topf=zeros(1,sections);
topd=zeros(1,sections);
bot=zeros(1,sections);
L=linspace(delta,max_time,max_time/delta);

%Displaying how many repetitions have been completed 
for r=1:Repetitions
if rem(r,Repetitions/10)==0
   r;
end


%%%%%%%%%%%%%%%%%%% Running each simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Running each simulation 
clear S I T incidence
startplace=N*0.5*rand; %random starting point with more than half of the population susceptible

%Run the simulation
[S,I,T,incidence]=gillespie_SIS(startplace,N,beta,gamma,max_time,delta,0); 

%Linearly space data points 
Lin = linspace(0,max_time,max_time*1000);
I_Lin = interp1(T,I,Lin);


%%%%%%%%%%%%%%%% Estimating the derivative in each section %%%%%%%%%%%%%%%


%The equation free method 
for t=1:length(I_Lin) - 1
    for s=1:sections
        if I_Lin(t)>=(s-1)*sectionwidth && I_Lin(t)<=s*sectionwidth %determining which section Z(t) lies in 
            topf(s)=topf(s) + ( I_Lin(t+1) - I_Lin(t) )/( 0.001 ); % \Delta t = 0.001
            %topd(s)=topd(s) + ((incidence(t+10)-incidence(t))^2)/( );
            bot(s)=bot(s)+1;
        end
    end
end
end %End of repetitions 


%%%%%%%%%%%%%%%%%%%%%%% Averaging EFM results %%%%%%%%%%%%%%%%%%%%%%%%%%


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

%Calculating approximation to potential surface
LM = linspace(0,M,sections);
potential_est = - cumtrapz(LM,F_est);

%Analytic Results 
LINSPACE = linspace(0,M,sections);
for y=1:sections
    x = LINSPACE(y)/N;
    F_x(y) = -N*(gamma*x - beta*x*(1-x)) ;
end
potential_ana = - cumtrapz(LINSPACE,F_x);


%%%%%%%%%%%%%%%%%%%%%%%% Creating figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure (1)
plot(LM/N,F_est,'b','DisplayName','Gillespie Simulated Drift')
hold on
plot(LINSPACE/N,F_x,'r--','DisplayName','SDE Analytic Drift')
hold off
legend
%axis([0 5000 -3500 300])
%Movi_1(b) = getframe(1);
 
Z(b,:) = potential_est;
Z2(b,:)= potential_ana;

figure (2)
plot(LM/N,potential_est,'b','DisplayName','Simulated Potential Surface')
hold on
plot(LINSPACE/N,potential_ana,'r--','DisplayName','Analytic Potential Surface')
hold off
legend
%axis([0 5000 -400000 2000000])
%Movi_2(b) = getframe(2);

%Creating the 3d plot
ub = max(max(Z));
lb = min(min(Z));
X = linspace(0,M/N,sections);
Xalt = [X M/N 0];
Y = [beta]*ones(1,sections);
Yalt = [Y beta beta];
Zalt = [Z(b,:) lb lb];
figure (3)
hold on
f1 = fill3(Xalt,Yalt,Zalt,'White');
xb=X( find(Z(b,:)==min(Z(b,:)),1) );
yb=beta;
zb=min(Z(b,:));
if beta>1
f1.EdgeColor = 'b';
plot3(xb,yb,zb,'b*')
else 
f1.EdgeColor = [0.8500 0.3250 0.0980];
plot3(xb,yb,zb,'Color',[0.8500 0.3250 0.0980],'Marker','*');
end
plot3(X,Y,Z2(b,:),'r--')
xb2=X( find(Z2(b,:)==min(Z2(b,:)),1) );
zb2=min(Z2(b,:));
plot3(xb2,yb,zb2,'r*')

%Calculating data on the polynomial fitting 
[coef,data] = poly_estimate(potential_est,LM/N,3,0.1,1);
poly_fit_data(b) = 6*coef(1)*xb + 2*coef(2);

end

%Creating the 3d plot - part 2
figure (3)
M=delta*beta0*N/4;
fill = plot3([M/sections;M/sections;M;M;M/sections],[1;1;1;1;1],[lb;ub;ub;lb;lb],'Color',[17 17 17]/100);
fill.LineStyle = '--';
xlabel('Prevalence')
ylabel('R_0')
zlabel('Potential Surface')
xlim([0 1/2])
zlim([lb ub])
set(gca,'Ydir','reverse')
grid on 

%Displaying Cubic Approximation Results 
figure (5)
hold on 
plot(Beta_vec,poly_fit_data,'DisplayName','Prevalence')
xlabel('R_0') 
ylabel('Power of x^2 in polynomial fitting')

%Creating the 3d plot - older version 
%{
X = [1/(2*sections):1/(2*sections):1/2];
Y = [beta1:-(beta1-beta0)/(nb-1):beta0];
Y = flip(Y);
figure (3)
hold on
w1 = waterfall(X',Y',Z);
CData_size1 = size(w1.CData,1);
for bb=1:nb
    xb=X( find(Z(bb,:)==min(Z(bb,:)),1) );
    yb=beta0 - (beta0-beta1)*((bb-1)/(nb-1));
    zb=min(Z(bb,:));
    figure (3)
    plot3(xb,yb,zb,'r*')
    
    if yb>1
        w1.CData(:,bb) = ones(CData_size1,1)*2*10^6;
    else
        w1.CData(:,bb) = ones(CData_size1,1)*8*10^6;
    end
end
xlabel('Prevalence')
ylabel('R_0')
zlabel('Potential Surface')
w1.CData(2,1) = 9*10^6;
ub = max(max(Z));
lb = min(min(Z));
zlim([lb ub])
fill = plot3([1/(2*sections);1/(2*sections);1/2;1/2;1/(2*sections)],[1;1;1;1;1],[lb;ub;ub;lb;lb],'Color',[17 17 17]/100);
fill.LineStyle = '--';
set(gca,'Ydir','reverse')
grid on 
%}

%Writing the movie 
%{
v_1 = VideoWriter('Prevalence Comparison for drift, N=10000, dt=1, beta 1.3-0.8.avi');
open(v_1)
writeVideo(v_1,Movi_1)
close(v_1)

v_2 = VideoWriter('Prevalence Comparison for potential N=10000, dt=1, beta 1.3-0.8.avi');
open(v_2)
writeVideo(v_2,Movi_2)
close(v_2)
%}