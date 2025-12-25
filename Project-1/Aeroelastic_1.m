%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Aeroelastic Analysis - 1
%   
%    SUDHI SHARMA P V       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all

Fs = 1000;               % Sampling Frequency
Fpass = 3;              % Passband Frequency
Fstop = 5;              % Stopband Frequency
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
StrFlow = '200RPM.txt';
Flow = '200RPM';
y = load(StrFlow);                 % input data


% Finding Linear Damping Ratio : Change Accordingly
numL1 = 8;  % Linear Damping Ratio number selected from graph


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

 %yfconv = (yfnew - 1.981)/ (-.0595) ;
 yfconv = 33.296 - (yfnew.*16.808);
 
 mean0rpm  = mean(yfconv)
 var0rpm = var(yfconv)

% Plotting the degree vs time

figure(4)
plot(xfnew,yfconv)
title1 = sprintf('Filtered data (%s)',Flow);
title(title1)
xlabel('time(s)')
ylabel('Pitch Amplitude,(\circ)')
grid
hold on

[pospeaks,s] = findpeaks(yfconv,xfnew);

% Aeroelastic Analysis
% Change peak number accordingly below

%175 - 13:40
% 225 - 55:62

switch Flow
    case 'Free Decay'
        starti = 8;
        endi   = 25;
        
    case '100RPM'
        starti = 13;
        endi   = 22;
       
    case '125RPM' 
        starti = 6;
        endi   = 19;
        
    case '150RPM'
        starti = 8;
         endi  = 30;
         
    case '175RPM'
        starti = 7;
        endi   = 31;
        
   case '200RPM'
        starti = 4;
        endi   = 38;

   case '225RPM'
        starti = 26;
        endi   = 31;

end


time = s(starti:endi);
peaks = pospeaks(starti:endi);

plot(time,peaks,'r')
st1 = sprintf('%s',Flow);
% saveas(gcf,st1,'epsc');
hold off

t = length(time);
p = length(peaks);

for i=1:length(peaks)-1
    
delta(i) = real(log(peaks(i)/peaks(i+1)));

zi(i) = delta(i)/ sqrt(4*pi^2 + delta(i)^2);

T(i) = time(i+1) - time(i);

time_instant(i) = time(i);

omegad_cycle(i) = 2*pi/T(i);

end

figure(10)
plot(peaks(1:p-1),zi, 'o')
ylim([0 inf]);
title2 = sprintf('Damping ratio Vs Pitch Amplitude (%s)',Flow);
title(title2)
xlabel('Pitch Amplitude,(\circ)')
ylabel('Damping Ratio,\xi')
st2 = sprintf('ziVsPitch_%s',Flow);
% saveas(gcf,st2,'epsc');

figure(11)
plot(time(1:t-1),zi, 'o')
ylim([0 inf]);
title3 = sprintf('Damping ratio Vs time (%s)',Flow);
title(title3)
xlabel('time,s')
ylabel('Damping Ratio,\xi')
st3 = sprintf('ziVstime_%s',Flow);
% saveas(gcf,st3,'epsc');

figure(12)
plot(time_instant,omegad_cycle, 'o')
ylim([0 40]);
title4 = sprintf('Damped Frequency Vs time (%s)',Flow);
title(title4)
xlabel('time,s')
ylabel('Damped Frequency,\omega_d (rad/s)')
st4 = sprintf('OmegadVstime_%s',Flow);
% saveas(gcf,st4,'epsc');


figure(13)
plot(peaks(1:p-1),omegad_cycle, 'o')
ylim([0 40]);
title5 = sprintf('Damped Frequency Vs Pitch Amplitude (%s)',Flow);
title(title5)
xlabel('Pitch Amplitude,(\circ)')
ylabel('Damped Frequency,\omega_d (rad/s)')
st5 = sprintf('OmegadVspitch_%s',Flow);
% saveas(gcf,st5,'epsc');


ziL = zi(numL1)

TimeL = time(numL1+1)-time(numL1)

Iea = 0.00110
 
omegaD = round((2 * pi/ TimeL),2)

omegaN = round((omegaD/sqrt(1-(ziL)^2)),2)

Ds = round((2 * ziL * omegaN * Iea),5)

Ks = round((omegaN ^2 * Iea),3)

ZiOmegaN = ziL * omegaN

OmegaHz = omegaD/(2*pi)

