%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Analytical estimate divergence of 3D wing

%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc

s = 10; % Span of wing, m

rho = 1.225; % Density of air Kg/m3

for n = 1:10

 eta(n) = n/10;
        
GJ_eta = 20000 * sqrt(1-(eta(n)/1.5)^4); % Nm2/deg

c = 2*((1-eta(n))*0.6 + 0.4);

e = 0.1105 * c; % eccentricty, distance of elastic axis from AC positive towards aft

ellipf = sqrt(1-eta(n)^2); % Elliptical function for CLalpha and Cmac 3D effects

alpha_0 = 4 - eta(n); % Geometric twist

U_div(n) = 1/(2*s)*(sqrt(180*GJ_eta/(rho*e*c*ellipf)));
        
       
end

U_div

plot(eta.*10,U_div)
title('Divergence Speed Along Span')
xlabel('Distance from root along span, y(m)')
ylabel('Divergence Speed, (m/s)')
line([1;9],[246.96;246.96],'linestyle','--','color','g')
% line([6;6],[0;450],'linestyle','--','color','g')




