
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Coupled Flutter analysis
%   Eigen Value Method

%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Inputs

M_h = 2.5; %Kg
M_theta = 0.77; % Kg
s = 0.61; % m
c = 0.156; % m
b = c/2; % m
rho = 1.225; % Kg/m^3
I_EA = 0.0011; %Kgm^2
K_h = 1000; %N/m
K_theta = 0.3; %Nm/rad
x_theta = 0.093; % m
EA = 0.27*c; %m
zeta = 0.05; %Damping Ratio
a_h = -0.46; % 

W_h = sqrt(K_h/M_h) % Heave Frequency, rad/s
W_theta = sqrt(K_theta/I_EA) % Pitch Frequency, rad/s

%  W_h =0;
%  W_theta = 16.51;

D_theta = 2*zeta*W_theta*I_EA;
D_h = 2*zeta*W_h*M_h;

syms U; %m/s

% Wagner Terms
A1 =0.165 ;
A2= 0.335;
b1 = 0.0455*U/b;
b2 = 0.3*U/b;

Udiv = sqrt(K_theta/((rho*b^2*s*2*pi)*(a_h+0.5)*(1-A1-A2)))

% M, C, and K matrices

M = [M_h+b^2*s*rho*pi,                   (M_theta*c*x_theta/2)-(b^2*s*rho*pi*b*a_h),        0;
    (M_theta*c*x_theta/2)-(b^2*s*rho*pi*b*a_h),  I_EA+((b^2*s*rho*pi)*(b^2*(1/8+a_h^2))),   0;
    0,                                        0,                                    1];

C = [D_h+rho*U*b*s*2*pi*(1-A1-A2),            (b^2*s*rho*pi*U)+(rho*U*b*s*2*pi)*(1-A1-A2)*b*(0.5-a_h),                                       (rho*U*b*s*2*pi)*(A1*b1+A2*b2);
    -((rho*U*b^2*s*2*pi)*(a_h+0.5)*(1-A1-A2)), (D_theta+(b^2*s*rho*pi*U)*b*(0.5-a_h)-((rho*U*b^2*s*2*pi)*(a_h+0.5)*(1-A1-A2)*b*(0.5-a_h))),   -((rho*U*b^2*s*2*pi)*(a_h+0.5)*(A1*b1+A2*b2));
    -1,                                              -b*(0.5-a_h),                                                                                (b1+b2)];                              


K = [K_h,  (rho*U^2*b*s*2*pi)*(1-A1-A2),                              (rho*U*b*s*2*pi)*b1*b2*(A1+A2);
    0,      (K_theta-((rho*U^2*b^2*s*2*pi)*(a_h+0.5)*(1-A1-A2))),      -((rho*U*b^2*s*2*pi)*(a_h+0.5)*b1*b2*(A1+A2));
    0,       -U,                                                                 b1*b2           ];

syms S;
A_s = S^2*[M]+S*[C]+[K];
poly = det(A_s)
roots = solve(poly,S)


figure(1)
title('Roots Vs Airspeed')
hold on

u = 0:.1:15;
fnum = 4;
a1 = 0;
a2 = 2;
   
for i = 1:length(u)

r = double(subs(roots,U,u(i)));
Z_1(i) = real(r(1));
Z_2(i) = real(r(2));

% Specific to this case since frequency interchanges at this point
if(i>1)
    if(abs(W_d1(i-1)-abs(imag(r(fnum+a1))))>= 2)
        a1 = 2;
        a2 = 0;
    end
end

W_d1(i) = abs(imag(r(fnum+a1)));
W_d2(i) = abs(imag(r(fnum+a2)));
DC1(i) = -real(r(fnum+a1));
DC2(i) = -real(r(fnum+a2));

u(i)
r 
end

plot(u,Z_1)
plot(u,Z_2)
plot(u,W_d1)
plot(u,W_d2)
plot(u,DC1)
plot(u,DC2)
hold off

ind = find(DC1 < 0,1);
fi = DC1(ind-1:ind);
ind2 = find(DC2 < 0,1);
fi2 = DC2(ind2-1:ind2);

U_f1 = Udiv; % SImply assigning a higher value than flutter speed at first, Divergence speed
U_f2 = Udiv;

if(fi ~= 0)
uf1 = u(ind-1:ind);
U_f1 = interp1(fi,uf1,0);
end

if(fi2 ~= 0)
uf2 = u(ind2-1:ind2);
U_f2 = interp1(fi2,uf2,0);
end

U_f = min(U_f1,U_f2)

grid on
figure(2)
plot(u,W_d1,u,W_d2)
title('Frequency Vs Airspeed')
xlabel('U, (m/s)')
ylabel('Frequency, \omega_d (rad/s)')
grid on 
hold on
line([U_f;U_f],[0;25],'linestyle','--','color','g')
line([0;15],[W_h;W_h],'linestyle','--')
line([0;15],[W_theta;W_theta],'linestyle','--','color','r')

figure(3)
plot(u,DC1,u,DC2)
title('Decay Rate Vs Airspeed')
xlabel('U, (m/s)')
ylabel('DEcay Rate, \beta')
grid on
line([0;15],[0;0],'linestyle','--','color','k')
line([U_f;U_f],[-5;20],'linestyle','--','color','g')
