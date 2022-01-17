%clear
%close all
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
tau_delay = zeros(length(Fs), 1);
for i = 1:length(Fs)
    tau_delay(i) =trapz(z, (1./(1-fpl2/Fs(i)^2).^(1/2) -1)/vel);
    for j = 1:length(Fd)
        phi1(i,j) = 1i*2*pi*(Fs(i)+Fd(j)/2)/vel*trapz(z, sqrt(1-fpl2/(Fs(i)+Fd(j)/2).^2)-1);
        phi2(i,j) = 1i*2*pi*(Fs(i)-Fd(j)/2)/vel*trapz(z, sqrt(1-fpl2/(Fs(i)-Fd(j)/2).^2)-1);
    end
end

T0 = 0.9*10^-8;
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
MeanFieldWoutCoh = p02.*exp(phi1-phi2);%.*squeeze(CoherenceEarth(:, ix,:));

MeanFieldFFT =  zeros(length(Fs),length(Fd));
Signal = zeros(length(Fs),length(Fd));

MeanFieldSum = zeros(1, length(Fd));
MeanFieldSum1 = zeros(1, length(Fd));
%MeanFieldWoutCohSum = zeros(1, length(Fd));
SignalSum = zeros(1, length(Fd));

for j = 1:length(Fd)
   MeanFieldSum(j) = trapz(Fs, MeanField(:, j));
   MeanFieldSum1(j) = trapz(Fs, abs(MeanField(:, j)));
   %MeanFieldWoutCohSum(j) = trapz(Fs, MeanFieldWoutCoh(:, j));
   SignalSum(j) = trapz(Fs, p02(:, j));
end

%% спектр сигнала
figure
plot(Fd/10^9, SignalSum, 'LineWidth', 2)
hold on
%plot(Fd, abs(MeanFieldWoutCohSum), 'LineWidth', 2)
plot(Fd/10^9, abs(MeanFieldSum1), 'LineWidth', 2)
legend({'input', 'output'})

set(gca, 'FontSize', 16)
grid minor
xlabel('f_d/f_s')
xlim([-0.8, 0.8])
%ylabel()

%% Средняя энергия
FFTSignalSum = fftshift(fft(SignalSum));
%FFTMeanFieldWoutCohSum = fftshift(fft(MeanFieldWoutCohSum));
FFTMeanFieldSum=fftshift(fft(MeanFieldSum));

fs_d = 1/abs((Fd(end)-Fd(end-1)));
T = abs((Fd(end)-Fd(1)));
t = (-fs_d/2:1/T:fs_d/2+1/T);

figure
plot(t*10^6, abs(FFTSignalSum), 'LineWidth', 2)
hold on
%plot(t, abs(FFTMeanFieldWoutCohSum), 'LineWidth', 2)
plot(t*10^6, abs(FFTMeanFieldSum), 'LineWidth', 2)
set(gca,'Fontsize', 16)


grid minor
xlabel('t, [ms]','Fontsize', 16)
xlim([-0.05, 5*10^(-1)])
legend({'input pulse', 'output pulse'})
set(gca,'Fontsize', 16)
