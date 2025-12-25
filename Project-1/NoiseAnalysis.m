%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Noise Properties analysis
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

%
% This program filters out the high frequency content of the data
%
dt = 1/1000;
dx = dt;                   % acquisition time step, ie dt = 1/Fs
y = load('0RPM.txt');                 % input data
n = length(y);

j = 0;

% Cutting the time signal to take data where no oscillations are present -
% for with flow case
for i=1:n
     x(i) = (i-1)*dx;
    if(x(i) >10 || x(i) == 10)
        j=j+1;
        ycut(j) = y(i);
        xcut(j) = x(i);
end    
end  

yf = filter(num,1,y);      % filtered output data

ncut = length(ycut);
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
        
    if(xf(i) > 0.1 || xf(i) == 0.1)
        j= j+1;
        yfnew(j) = yf(i);
        xfnew(j) = xf(i);
    end    
end
 
% PSD of data
 [Pxx,F] = periodogram (yfnew, [], length(yfnew),1000);
 
 % PDF of data
 xnorm = 1.95:0.001:2.01;
 norm = normpdf(xnorm,mean(yfnew),sqrt(var(yfnew)));

 
sprintf('%s %.5f', 'Mean of filtered data', mean(yfnew))
sprintf('%s %.5f', 'variance of filtered data', var(yfnew))

%Plotting

figure(1)
hist(y,100);
title('Histogram of noisy data')

  
 figure(3)
 plot(x,y)
 hold on
 plot(xf,yf,'r')
 title('Filtered data over noisy data')
 xlabel('Time(s)')
 ylabel('Voltage(volts)')
 
 figure(5)
 plot(F,10*log10(Pxx))
 title('PSD of noise filtered data')
 xlabel('frequency(Hz)')
 ylabel('Voltage^2/frequency')
 xlim([0 20]);

figure(6)
hist(yfnew,100)
grid
hold on
plot(xnorm,norm.*12,'r')
title('PDF of filtered data using Normal distribution')
xlabel('Voltage')
ylabel('Probability,P(x)')

