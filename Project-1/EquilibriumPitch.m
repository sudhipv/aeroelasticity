%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Equilibrium Pitch
%  
%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

Fs = 1000;               % Sampling Frequency
Fpass = 30;              % Passband Frequency
Fstop = 35;              % Stopband Frequency
Dpass = 0.057501127785;  % Passband Ripple
Dstop = 0.0001;          % Stopband Attenuation
flag  = 'noscale';       % Sampling Flag
% Calculate the order from the parameters using KAISERORD.
[N,Wn,BETA,TYPE] = kaiserord([Fpass Fstop]/(Fs/2),[1 0],[Dstop Dpass]);
% Calculate the coefficients using the FIR1 function.
num = fir1(N, Wn, TYPE, kaiser(N+1, BETA), flag);  % Filter parameters 
GD = N/2;                % Offset delay created by filtering
clear BETA Dpass Dstop Fpass Fs Fstop N TYPE Wn flag


% This program filters out the high frequency content of the data
%
dt = 1/1000;
dx = dt;                   % acquisition time step, ie dt = 1/Fs
y = load('175RPM.txt');                 % input data
n = length(y);

yf = filter(num,1,y);      % filtered output data

% Resets the filtered data to the same time as the input data
for i = 1:n
    x(i) = (i-1)*dx;
    xf(i) = (i-1-GD)*dx;
    
 end


%Cutting Yf to have only stabilised noise and remove artificial content
%from FFT

j = 0;
nf = length(yf);
for i = 1:nf
        
    if(xf(i) > 0.2 || xf(i) == 0.2)
        j= j+1;
        yfnew(j) = yf(i);
        xfnew(j) = xf(i);
    end    
end

mean(yfnew)
var(yfnew)

% Removing the noise from filtered data

yfconv = 33.296 - (yfnew.*16.808);

meanyfconv = mean(yfconv)
varyfconv = var(yfconv)

% Plotting the degree vs time

figure(4)
plot(xfnew,yfconv)
title('Filtered data')
xlabel('time')
ylabel('Pitch Amplitude,(in \circ)')
grid


% Cutting signal to take values after oscillations , Change time for different oscillations

j= 0;
ncut = length(yfconv);
for i = 1:ncut
        
    if(xfnew(i) > 13 || xfnew(i) == 13)
        j= j+1;
        yfcut(j) = yfconv(i);
        xfcut(j) = xfnew(i);
    end    
end

meanyfcut = mean(yfcut)
varyfcut = var(yfcut)

figure(5)
plot(xfcut,yfcut)
title('Filtered data')
xlabel('time')
ylabel('Pitch Amplitude,(in \circ)')
grid



 % PDF of data
 xnorm = -0.8:0.001:0.8;
 norm = normpdf(xnorm,meanyfcut,sqrt(varyfcut));

figure(6)
hist(yfcut,100)
grid
hold on
plot(xnorm,norm.*35,'r')
title('PDF of filtered data using Gaussian distribution')
xlabel('Pitch Amplitude')
ylabel('Probability,P(x)')


zeropitch = -0.0071;
% 0.1399,-0.3794,0.0058,0.0570,0.0578,-0.4737
eqpitch = [-0.0071,0.1436,-0.3705,0.0062,0.0554,0.0600,-0.4740];

eqpitch = eqpitch - zeropitch;

U = [0, 4.6, 5.8, 7.0, 8.2, 9.4, 10.6];

figure(10)
plot(U,eqpitch,'o')
title('Equilibrium Pitch Vs Air Speed')
xlabel('U(m/s)')
ylabel('Equilibrium Pitch, in \circ')
ylim([-1 0.4])
