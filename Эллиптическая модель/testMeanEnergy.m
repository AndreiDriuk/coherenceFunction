clear
close all
%load('CohEarth_sigma05_fd_002.mat');
load('CohEarth10_8.mat');
lt = 140660; % coefficient in cosh(z/lt) 
lb = 198722; % coefficient in 1-z^2/lb^2
lend = 532000;
z = [-1000000:1000:lend,lend];
kplm =  0.2487;

charge = 1.602E-19;
vel = 299792458;
diel = 8.85E-12;
masse = 9.1E-31;
kplTop2 = kplm^2./cosh(z(z<=0)/lt);
kplBot2 = kplm^2*sqrt(1-z((z>0)&(z<=lb)).^2/lb^2);
kplUnder2 =z(z>lb)*0;
Fs = (0.8:0.005:1.2)*10^9;
Fd = ((-10:0.002:10)*10^8)'; %% t in DFT
ix = round(length(x)/2);

fpl2 = [kplTop2, kplBot2, kplUnder2]*vel^2/(4*pi^2);
f = 10^9;
phi1 = zeros(length(Fs), length(Fd));
phi2 = zeros(length(Fs), length(Fd));
for i = 1:length(Fs)
    for j = 1:length(Fd)
        phi1(i,j) = 1i*2*pi*(Fs(i)+Fd(j)/2)/vel*trapz(z, sqrt(1-fpl2/(Fs(i)+Fd(j)/2).^2)-1);
        phi2(i,j) = 1i*2*pi*(Fs(i)-Fd(j)/2)/vel*trapz(z, sqrt(1-fpl2/(Fs(i)-Fd(j)/2).^2)-1);
    end
end

T0 = 2*10^-8;
fc = 1*10^9;

p02 = zeros(length(Fs), length(Fd));
for i = 1:length(Fs)
    for j = 1:length(Fd)
       p02(i,j) = T0^2/(2*pi)*exp(-T0^2/2*(fc-Fs(i))^2)*exp(-T0^2*Fd(j)^2/4);
    end
end

FsL = length(Fs);
FdL = length(Fd);
FFt = zeros(FdL, FdL);
i1 = 0:FdL-1;
j1= 0:FdL-1;

MeanField = p02.*exp(phi1-phi2).*squeeze(CoherenceEarth(:, ix,:));%p02.*
MeanField1 = exp(phi1-phi2);%p02.*exp(phi1-phi2);%.*squeeze(CoherenceEarth(:, ix,:));
MeanField2= squeeze(CoherenceEarth(:, ix,:));%p02.*exp(1i*tau_delay);%.*squeeze(CoherenceEarth(:, ix,:));; 
MeanField3 = p02.*exp(phi1-phi2);
MeanField4 = p02.*squeeze(CoherenceEarth(:, ix,:));

MeanFieldFFT =  zeros(length(Fs),length(Fd));
Signal = zeros(length(Fs),length(Fd));

MeanFieldSum = zeros(1, length(Fd));
MeanFieldSum1 = zeros(1, length(Fd));
MeanFieldSum2 = zeros(1, length(Fd));
MeanFieldSum3 = zeros(1, length(Fd));
MeanFieldSum4 = zeros(1, length(Fd));
SignalSum = zeros(1, length(Fd));
% 
MeanFieldSumFd = zeros(1, length(Fs));
MeanFieldSumFd1 = zeros(1, length(Fs));
MeanFieldSumFd2 = zeros(1, length(Fs));
MeanFieldSumFd3 = zeros(1, length(Fs));
MeanFieldSumFd4 = zeros(1, length(Fs));
SignalSumFd = zeros(1, length(Fs));

for j =1:length(Fd)
   MeanFieldSum(j) = trapz(Fs, MeanField(:, j));
   MeanFieldSum1(j) = trapz(Fs, MeanField1(:, j));
   MeanFieldSum2(j) = trapz(Fs, MeanField2(:, j));
   MeanFieldSum3(j) = trapz(Fs, MeanField3(:, j));
   MeanFieldSum4(j) = trapz(Fs, MeanField4(:, j));
   SignalSum(j) = trapz(Fs, p02(:, j));
end

for j = 1:length(Fs)
   MeanFieldSumFd(j) = trapz(Fd, abs(MeanField(j, :)));
   MeanFieldSumFd1(j) = trapz(Fd, abs(MeanField1(j, :)));
   MeanFieldSumFd2(j) = trapz(Fd, abs(MeanField2(j, :)));
   MeanFieldSumFd3(j) = trapz(Fd, abs(MeanField3(j, :)));
   MeanFieldSumFd4(j) = trapz(Fd, abs(MeanField4(j, :)));
   SignalSumFd(j) = trapz(Fd, p02(j, :));
end

%% спектр сигнала
 figure
 plot(Fd, SignalSum, 'LineWidth', 2)
 hold on
 plot(Fd, abs(MeanFieldSum), 'LineWidth', 2)
 %plot(Fd, abs(MeanFieldSum1), 'LineWidth', 2)
 %plot(Fd, abs(MeanFieldSum2), 'LineWidth', 2)
 %plot(Fd, abs(MeanFieldSum3), 'LineWidth', 2)
 %plot(Fd, abs(MeanFieldSum4), 'LineWidth', 2)
 
%  figure
%  plot(Fd, SignalSum, 'LineWidth', 2)
%  hold on
%  plot(Fd, unwrap(angle(MeanFieldSum)), 'LineWidth', 2)
%  %plot(Fd, unwrap(angle(MeanFieldSum1)), 'LineWidth', 2)
%  %plot(Fd, unwrap(angle(MeanFieldSum2)), 'LineWidth', 2)
%  %plot(Fd, unwrap(angle(MeanFieldSum3)), 'LineWidth', 2)
%  %plot(Fd, unwrap(angle(MeanFieldSum4)), 'LineWidth', 2)
%  
 figure
 plot(Fs, SignalSumFd, 'LineWidth', 2)
 hold on
 plot(Fs, abs(MeanFieldSumFd), 'LineWidth', 2)
 %plot(Fs, abs(MeanFieldSum1), 'LineWidth', 2)
 %plot(Fs, abs(MeanFieldSum2), 'LineWidth', 2)
 %plot(Fs, abs(MeanFieldSumFd3), 'LineWidth', 2)
 %plot(Fs, abs(MeanFieldSumFd4), 'LineWidth', 2)
 
%  figure
%  plot(Fs, SignalSumFd, 'LineWidth', 2)
%  hold on
%  plot(Fs, unwrap(angle(MeanFieldSumFd)), 'LineWidth', 2)
%  %plot(Fs, unwrap(angle(MeanFieldSum1)), 'LineWidth', 2)
%  %plot(Fs, unwrap(angle(MeanFieldSum2)), 'LineWidth', 2)
%  %plot(Fs, unwrap(angle(MeanFieldSumFd3)), 'LineWidth', 2)
% % plot(Fs, unwrap(angle(MeanFieldSumFd4)), 'LineWidth', 2)

figure()
plot(Fd, p02(round(length(Fs)/2), :))
hold on
plot(Fd, abs(MeanField(round(length(Fs)/2), :)))

figure()
plot(Fs, p02(:, 1))
hold on
plot(Fs, abs(MeanField(:, 1)))