%This code is created for use in the MMath Project: Investigating the
%potetntial of early warning signals in disease elimination. 

%Here the drift function and potential surface for incidence data is
%estimated from simulated data. 

clear 


%Parameters for the gillespie simulation 
Repetitions = 100; %number of simulations
N=10000; %population size
beta0 = 1.3; %infectivity constant
beta1 = 0.9;
gamma = 1; %recovery rate
max_time = 10; %length of each simulation
delta=0.1; %incidence aggregation timestep 

%Repeating the model with varing beta values 
nb = 1; %the number of different beta values 
poly_fit_data = zeros(nb,1);
Beta_vec = zeros(nb,1);
for b = 1:nb
clear M topf topd F_est D_est potential_est
if nb>1
beta = beta0 - (beta0-beta1)*((b-1)/(nb-1))
else 
    beta = beta0
end
Beta_vec(b) = beta;

%Resetting the equation free method for each beta value 
M = (delta*beta*N/4)*1; %The potetial surface is constructed for the range [0,M]
sections=50; %Number of sections into which the range of incidence is split for the equation free method 
sectionwidth=M/sections; 
topf=zeros(1,sections);
topd=zeros(1,sections);
bot=zeros(1,sections);
L=linspace(delta,max_time,max_time/delta);

%Main iteration 
for r=1:Repetitions
%Showing how many repetitions have been completed   
if rem(r,Repetitions/10)==0
   r;
end
   
%Running each simulation 
clear S I T incidence
%Aim to have the start point unformly distributed on [0,M]
startplace=N*0.5*rand; %random starting point with more than half of the population susceptible
[S,I,T,incidence]=gillespie_SIS(startplace,N,beta,gamma,max_time,delta,0); 


%The equation free method 
%Use this code for the forward difference
for t=1:length(incidence)-1
    for s=1:sections
        if incidence(t)>=(s-1)*sectionwidth && incidence(t)<=s*sectionwidth %determining which section Z(t) lies in 
            topf(s)=topf(s) + (incidence(t+1)-incidence(t))/(1*delta);
            bot(s)=bot(s)+1;
        end
    end
end

%Use this code for the central difference 
% % %The equation free method 
% % for t=2:length(incidence)-1
% %     for s=1:sections
% %         if incidence(t)>=(s-1)*sectionwidth && incidence(t)<=s*sectionwidth %determining which section Z(t) lies in 
% %             topf(s)=topf(s) + (incidence(t+1)-incidence(t-1))/(2*delta);
% %             bot(s)=bot(s)+1;
% %         end
% %     end
% % end

end %end of repetitons/main iteration 

%Calculating the simulated estimate to the drift function 
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

%Analytic Results 
LINSPACE = linspace(0,delta*beta*N/4,sections);
for y=1:sections
    u = LINSPACE(y);
    p = (1 - (4*u)/(beta*delta*N))^0.5; %useful in the following computations 
    F_u(y) = beta*N*delta*p*( u/(N*delta) - gamma/2 + (gamma/2)*p)  - beta*delta*(gamma/2 - (gamma/2)*p + u/(N*delta));
   % F_u_alt(y) = - gamma*u + beta*u*(1 - u/(N*gamma*dt));
    %V(y) = b*(b*N/6)*l*(p^3) + (b*(b*N)*(b*N)/60)*(p^5) + (b*g - b*g*N/2)*(b*N/6)*(p^3) + (2*g + 2*b/N)*0.5*l*l - (b*g*N/2 - b*g)*l;
end

%Calculating the Potential surface
LM = linspace(0,M,sections);
potential_est = - cumtrapz(LM,F_est);
potential_ana = - cumtrapz(LINSPACE,F_u);

Z(b,:) = potential_est;
Z2(b,:)= potential_ana;

%Creating Figures 
figure (1)
plot(LM,F_est,'b','DisplayName','Simulated Drift')
hold on
plot(LINSPACE,F_u,'r--','DisplayName','Analytic Drift')
hold off
%hold off
legend
%axis([0 35000 -2500 2000])
%Movi_1(b) = getframe(1);

figure (2)
plot(LM,potential_est,'b','DisplayName','Simulated Potential Surface')
hold on
plot(LINSPACE,potential_ana,'r--','DisplayName','Analytic Potential Surface')
hold off
legend
%axis([0 35000 -3000000 1000000])
%Movi_2(b) = getframe(2);


%Creating the 3d plot - version 2 
ub = max(max(Z));
lb = min(min(Z));
X = linspace(0,M,sections);
Xalt = [X M 0];
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
[coef,data] = poly_estimate(potential_est,LINSPACE,3,0.1,1);
poly_fit_data(b,c) = 6*coef(1)*xb + 2*coef(2);
figure (4)
plot(xb,zb,'k*')

end

figure (3)
xlabel('Incidence')
ylabel('R_0')
zlabel('Potential Surface')
M = delta*beta0*N/4;
xlim([0 M])
zlim([lb ub])
set(gca,'Ydir','reverse')
grid on 

%Displaying Cubic Approximation Results 
figure (5)
hold on 
plot(Beta_vec,poly_fit_data,'DisplayName','Incidence')
xlabel('R_0')
ylabel('Power of x^2 in polynomial fitting')

%Writing the Movie: used when testing varying values of beta 
%{
v_1 = VideoWriter('Blank','MPEG-4');
open(v_1)
writeVideo(v_1,Movi_1)
close(v_1)

v_2 = VideoWriter('Blank','MPEG-4');
open(v_2)
writeVideo(v_2,Movi_2)
close(v_2)
%}