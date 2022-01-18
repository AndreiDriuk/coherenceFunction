%% Скрипт для расчета средней энергии импульса без координат
clear
close all
addpath('../Данные/')
%load('CohEarth_sigma05_fd_002.mat');
load('CohEarth10_8.mat');
%% параметры трассы
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
%% определение суммарной и разнсотной частот
Fs = (0.8:0.005:1.2)*10^9;
Fd = ([-10:0.002:10]*10^8)'; %% t in DFT
ix = round(length(x)/2); % uncomment for Coherence

%% экспоненциальный множитель с вычетом 1
fpl2 = [kplTop2, kplBot2, kplUnder2]*vel^2/(4*pi^2);
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

%% определение импульса 
T0 = 10^-8; % длительность импульса
fc = 1*10^9; % Центральная частота импульса

p02 = zeros(length(Fs), length(Fd));
for i = 1:length(Fs)
    for j = 1:length(Fd)
       p02(i,j) = T0^2/(2*pi)*exp(-T0^2/2*(fc-Fs(i))^2)*exp(-T0^2*Fd(j)^2/4);
    end
end

%% подынтегральное выражение 
MeanField = p02.*exp(phi1-phi2).*squeeze(CoherenceEarth(:, ix,:));
%MeanFieldWoutCoh = p02.*exp(phi1-phi2); checking

MeanFieldFFT =  zeros(length(Fs),length(Fd));
Signal = zeros(length(Fs),length(Fd));

MeanFieldSum = zeros(1, length(Fd));
%MeanFieldWoutCohSum = zeros(1, length(Fd)); checking
SignalSum = zeros(1, length(Fd));
%% Первоачальной суммирование по суммарной частоте (первый интеграл)
for j = 1:length(Fd)
   MeanFieldSum(j) = trapz(Fs, MeanField(:, j));
   %MeanFieldWoutCohSum(j) = trapz(Fs, MeanFieldWoutCoh(:, j)); %checking
   SignalSum(j) = trapz(Fs, p02(:, j));
end

%% спектр энегрии
%SignalSum - исходный p0.^2
%MeanFieldSum - конечный p0.^2
figure
plot(Fd/10^6, SignalSum, 'LineWidth', 2)
hold on
%plot(Fd, abs(MeanFieldWoutCohSum), 'LineWidth', 2)
plot(Fd/10^6, abs(MeanFieldSum), 'LineWidth', 2)
legend({'input', 'output'})

set(gca, 'FontSize', 16)
grid minor
xlabel('f_d MHz')
%xlim([-0.8, 0.8])

%% Средняя энергия
% FFTSignalSum - исходный импульс
% FFTMeanFieldSum - конечный импульс
FFTSignalSum = abs(fftshift(fft(SignalSum)));
%FFTMeanFieldWoutCohSum = fftshift(fft(MeanFieldWoutCohSum)); checking
FFTMeanFieldSum=abs(fftshift(fft(MeanFieldSum)));

%% переход во временную область и определение параметров
fs_d = 1/abs((Fd(end)-Fd(end-1)));
T = abs((Fd(end)-Fd(1)));
t = (-fs_d/2:fs_d/length(Fd):fs_d/2-fs_d/length(Fd));

figure
plot(t*10^6, FFTSignalSum, 'LineWidth', 2)
hold on
%plot(t, abs(FFTMeanFieldWoutCohSum), 'LineWidth', 2)
plot(t*10^6, FFTMeanFieldSum, 'LineWidth', 2)
set(gca,'Fontsize', 16)
grid minor
xlabel('t, [µs]','Fontsize', 16)
xlim([-0.05, 5*10^(-1)])
legend({'input pulse', 'output pulse'})
set(gca,'Fontsize', 16)

%% проверка сохранности площади под энергией импульса во времени
S1 = trapz(t, FFTSignalSum);
S2 = trapz(t, FFTMeanFieldSum);
disp(['Отношение площадей  исходного и конечного импульса ', num2str(S1/S2*100)])

%% Переход обратно в частотную область  (но уже от модуля энергии импульса)
% 1. Определение параметров
fs1 = 1/(fs_d/length(Fd));
L = length(t);
w = (-fs1/2:fs1/L:fs1/2-fs1/L);

%% обратное Фурье преобразование
SS = fftshift(fft(FFTSignalSum));
SS1 = fftshift(fft(FFTMeanFieldSum));

%% Проверка Равенства Парсеваля
% не работает так
% S11 = trapz(w, abs(SS.^2))/length(w);
% S12 = trapz(w, abs(SS1.^2))/length(w);
% 
% S21 = trapz(t, abs(FFTSignalSum.^2));
% S22 = trapz(t, abs(FFTMeanFieldSum.^2));
% а так работает
S11 = sum(abs(SS.^2))/length(w);
S12 = sum(abs(SS1.^2))/length(w);

S21 = sum( abs(FFTSignalSum.^2));
S22 = sum( abs(FFTMeanFieldSum.^2));

disp(['Отношение квадрата площадей исходного импулса во временной и частотной области ', num2str(S11/S21)])
disp(['Отношение квадрата площадей выходного импулса во временной и частотной области ', num2str(S12/S22)])
