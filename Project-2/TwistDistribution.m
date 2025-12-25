
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Twist Distribution along the span - 3D Wing
%   Assumed Mode method
%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear
syms eta;

GJ_eta = 20000 * sqrt(1-(eta/1.5)^4)*(180/pi); % Nm2/rad

s = 10; % Span of wing, m

rho = 1.225; % Density of air Kg/m3

U = 200 ; % Air speed, m/s

q = .5 * rho * U^2; % Dynamic pressure

c = 2*((1-eta)*0.6 + 0.4);

ellipf = sqrt(1-eta^2); % Elliptical function for CLalpha and Cmac 3D effects

e = 0.1105 * c; % eccentricty, distance of elastic axis from AC positive towards aft

alpha_0 = (4 - eta)*(pi/180); % Geometric twist, in rad

Cmac = 0.05;


for nummodes = 1:3 % Number of modes to represent pitching 
        
        %f(nummodes) = eta^(nummodes); % Polynomial Modes
        f(nummodes) = sin((2*nummodes -1)*pi*eta/2); % Normal Modes, Cantilever beam in torsion
        f_bar(nummodes) = diff(f(nummodes),eta)
                  
for i=1:nummodes
    
    for j=1:nummodes
        
        E_expression = 1/s*(GJ_eta*f_bar(i)*f_bar(j));
        Kbar_expression = s*q*(2*pi*ellipf*e*c*f(i)*f(j));
        
        E(i,j) = double(int(E_expression,0,1));    
         Kbar(i,j) = double(int(Kbar_expression,0,1));
    end

 F_expression = s*q*ellipf*(e*c*2*pi*alpha_0*f(i)+ c*c*Cmac*f(i));
 
 F(i) = double(int(F_expression,0,1)); 
     
end
x = 0:0.1:1;

theta_bar = inv(E-Kbar)*F'

theta_eta = f*theta_bar

theta_new = subs(theta_eta,x)

hold on
plot(x.*10,theta_new*(180/pi))
xlabel('Distance from root along span, y(m)')
ylabel('Pitch angle,\theta(\circ)')
title('Twist Distribution along span for Torsional Modes')
end




