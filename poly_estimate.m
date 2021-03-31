%Taking a polynomial estimate of the bottom 10% of the potential surface

function [coef,data] = poly_estimate(v,x,degree,proportion,plots)
% v is the potential surface 
% x is the x values at which the potential surface is calculated 
% d is the degree of the polynomial fitting, usually 2 
% p is the proportion of the surface that will be fitted 

%First determine where the fitting is to be taken 
v_min = min(v);
v_max = max(v);
position_min = find(v==v_min,1);
x_min=x(position_min); %The location of the minimum of v 
range = v_max - v_min;
fitting_range = proportion*range; 
fitting_height = v_min + fitting_range; %The fitting is taken over the vertical range (v_min,fitting_height)

%Now determine what values of v are closest to fitting_height 
v_before_min = v(1:position_min);
v_after_min = v(position_min:end);
left = find( abs(v_before_min-fitting_height)==min(abs(v_before_min-fitting_height)),1 );
right = find( abs(v_after_min-fitting_height)==min(abs(v_after_min-fitting_height)),1 );

%The polynomial fitting will be done on x(LEFT:RIGHT),v(LEFT:RIGHT)
LEFT = left;
RIGHT = position_min + right - 1;

[coef] = polyfit(x(LEFT:RIGHT),v(LEFT:RIGHT),degree);
data = coef(1);

poly = polyval(coef,x(LEFT:RIGHT));

if plots==1
figure (4)
hold on
plot(x,v)
plot(x(LEFT:RIGHT),poly,'r')
hold off
end
end