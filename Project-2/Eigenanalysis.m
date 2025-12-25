
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Eigen analysis of divergence
%   Assumed Mode Method

%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
syms eta;

GJ_eta = 20000 * sqrt(1-(eta/1.5)^4)*(180/pi); % Nm2/rad

s = 10; % Span of wing, m

rho = 1.225; % Density of air Kg/m3

c = 2*((1-eta)*0.6 + 0.4); % Trapezoidal Wing, variation of chord, 0.8m at tip and 2m at root

ellipf = sqrt(1-eta^2); % Elliptical function for CLalpha and Cmac 3D effects

e = 0.1105*c; % eccentricty, distance of elastic axis from AC positive towards aft

alpha_0 = (4 - eta)*(pi/180); % Geometric twist


for nummodes = 1:3 % Number of modes to represent pitching 

        %f(nummodes) = eta^(nummodes) % Polynomial Modes
        f(nummodes) = sin((2*nummodes -1)*pi*eta/2); % Normal Modes, Cantilever beam in torsion
        f_bar(nummodes) = diff(f(nummodes),eta);
end     
            
for i=1:nummodes
    
    for j=1:nummodes
        
        E_expression = (1/s)*(GJ_eta*f_bar(i)*f_bar(j));
        K_bar_expression = (s)*(2*pi*ellipf*e*c*f(i)*f(j));
        

         E(i,j) = double(int(E_expression,0,1)) ;   
         K_bar(i,j) = double(int(K_bar_expression,0,1));

    end

end

% Multi Modes

% GENERALIZED EIGENPROBLEM AX=λBX
  for i=1:nummodes 
      
      A = E(1:i,1:i);
      B = K_bar(1:i,1:i);
      lam = eig(A,B);
      %lam = eig(inv(K_bar(1:i,1:i))*E(1:i,1:i)) 
      % E = structural stiffness,  K_bar = -K/q matrix 
      % λ is the dynamic pressure (q)

% AIRSPEED CALCULATIONS
Udiv = sqrt(2*lam/rho);

U_div(i) = min(Udiv)


end
