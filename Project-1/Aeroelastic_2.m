%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Aeroelastic Analysis -2
% Plotting Decay Rate, Damping Coefficient, Stiffness COefficient

%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OmegaDamp = [23.80,22.28,21.23,20.67,18.87,16.71,12.62];

ZiOmega = [0.607,1.20,1.72,2.21,2.90,4.32,5.92];

DampCoeff = [0.0013,0.0027,0.0039,0.0049,0.0066,0.0095,0.0130];

StiffCoeff = [0.624,0.548,0.499,0.475,0.402,0.328,0.214];

StiffCoeffAnalytical = [0.624,0.5286,0.4573,0.3695,0.2653,0.1446,0.0075];

OmegaN = [23.81,22.31,21.30,20.79,19.11,17.26,13.94];

OmegaNAnalytical = [23.81,21.921,20.389,18.32,15.53,11.46,2.611];

U = [0, 4.6, 5.8, 7.0, 8.2, 9.4, 10.6];

for i=1:length(U)

Dae(i) = 0.0013 + 2.005*10^(-4)*U(i);

end

figure(11)
plot(U.^2, OmegaDamp.^2,'*')
title('Variation of Damped Frequency with Airspeed')
xlabel('Air Speed, U^2(m^2/s^2)')
ylabel('Damped Frequency,\omega_d^2 (rad^2/s^2)')
ylim([0 inf])
xlim([0 200])

figure(12)
plot(U, ZiOmega,'*')
title('Variation of Damping Rate with Airspeed')
xlabel('Air Speed, U(m/s)')
ylabel('Damping Rate,\xi\omega_n(rad/s)')
ylim([0 inf])

figure(13)
plot(U, DampCoeff,'*')
title('Variation of Damping Coefficient with Airspeed')
xlabel('Air Speed, U(m/s)')
ylabel('Damping Coefficient, D_a_e(Nm rad s)')
ylim([0 inf])
xlim([0 15])

figure(14)
plot(U, StiffCoeff,'*')
title('Variation of stiffness coefficient with Airspeed')
xlabel('Air Speed, U(m/s)')
ylabel('Stiffness Coefficient, K_a_e(Nm rad^2)')
ylim([0 inf])
xlim([0 15])

figure(15)
plot(U, StiffCoeff,'o',U, StiffCoeffAnalytical,'*')
title('Variation of stiffness coefficient with Airspeed')
xlabel('Air Speed, U(m/s)')
ylabel('Stiffness Coefficient, K_a_e(Nm rad^2)')
ylim([0 inf])
xlim([0 15])

figure(16)
plot(U, OmegaN,'o',U, OmegaNAnalytical,'*')
title('Variation of frequency with Airspeed')
xlabel('Air Speed, U(m/s)')
ylabel('Natural Fequency \omega_n(rad/s)')
ylim([0 inf])
xlim([0 15])

figure(17)
plot(U.^2, StiffCoeff,'o',U.^2,StiffCoeffAnalytical,'*')
title('Variation of Aeroelastic stiffness with Airspeed')
xlabel('Air Speed, U^2 (m^2/s^2)')
ylabel('Aeroelastic Stiffness,K_a_e(Nm rad^2)')
ylim([0 inf])
xlim([0 200])
hold off

figure(18)
plot(U, DampCoeff,'o',U,Dae,'*')
title('Variation of Damping Coefficient with Airspeed')
xlabel('Air Speed, U (m/s)')
ylabel('Damping Coefficient,D_a_e(Nm rad s)')
ylim([0 inf])
xlim([0 inf])

figure(19)
plot(U, OmegaDamp,'o')
title('Variation of Damped Frequency with Airspeed')
xlabel('Air Speed, U (m/s)')
ylabel('Damped Frequency, \omega_d(rad/s)')
ylim([0 inf])
xlim([0 20])


