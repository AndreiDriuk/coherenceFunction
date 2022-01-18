%% Скрипт для расчета средней энергии импульса с учетом координат
clear
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


fpl2 = [kplTop2, kplBot2, kplUnder2]*vel^2/(4*pi^2);
f = 10^9;

phi1 = zeros(length(Fs), length(Fd));
phi2 = zeros(length(Fs), length(Fd));

for i = 1:length(Fs)
    for j = 1:length(Fd)
        phi1(i,j) = 1i*2*pi*(Fs(i)+Fd(j)/2)/vel*trapz(z, (sqrt(1-fpl2/(Fs(i)+Fd(j)/2)^2)));
        phi2(i,j) = 1i*2*pi*(Fs(i)-Fd(j)/2)/vel*trapz(z, (sqrt(1-fpl2/(Fs(i)-Fd(j)/2)^2)));
    end
end

T0 = 10^-8;
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

 for ix1 = 1:FsL
      FFt(ix1,:) = exp(-2*pi*1i*i1(ix1)*j1/FdL);
 end


fs_d = 1/abs((Fd(end)-Fd(end-1)));
T = (Fd(end)-Fd(1));
t = (-fs_d/2:1/T:fs_d/2+1/T);
MeanField = zeros(length(Fs), length(x),length(Fd));

for ix = 1:length(x)
%    for iy =1:length(y)
     MeanField(:, ix,:) = p02.*exp(phi1-phi2).*squeeze(CoherenceEarth(:, ix,:));  
%    end
end

SumMeanField = zeros(length(x), length(fd));
for ix = 1:length(x)
    for iy =1:length(t)
     SumMeanField(ix, iy) = trapz(Fs, MeanField(:,ix, iy)); 
    end
end

MeanFieldFFT =  zeros(length(x),length(t));

for ix = 1:length(x)
     MeanFieldFFT(ix, :) = fftshift(fft(SumMeanField(ix,:))); 
end

MeanFieldFFT1 =  zeros(length(x),length(t));

for it = 1:length(t)
     MeanFieldFFT1(:, it) =fftshift(fft(MeanFieldFFT(:,it))); 
end
x_d = 1/abs((x(end)-x(end-1)));
TX = (x(end)-x(1));
tx = (-x_d/2:1/TX:x_d/2);

contour(x/1000, t*10^6, abs(MeanFieldFFT)', 'LineWidth', 2)
ylim([0.0, 0.2])
xlim([-5, 5 ])
xlabel("km", 'Fontsize', 16)
ylabel('Delay, ms', 'Fontsize', 16)
grid on
xticks(-5:1:5)
set(gca,'Fontsize', 16)
% figure(2)
% for i = 1:length(Fs)
%     for j = 1:length(Fd)
%         phi1(i,j) = (Fs(i)+Fd(j)/2)*trapz(z, sqrt(sqrt(1-fpl2/(Fs(i)+Fd(j)/2)^2)));
%         phi2(i,j) = (Fs(i)-Fd(j)/2)*trapz(z, sqrt(sqrt(1-fpl2/(Fs(i)-Fd(j)/2)^2)));
%     end
% end
% 
% %plot(Fd, phi1(5,:)-phi2(5,:))
% 
% phi3 = zeros(length(Fs), length(Fd));
% for i = 1:length(Fs)
%     for j = 1:length(Fd)
%         phi3(i,j) = Fd(j)*trapz(z, (1/2*(1-fpl2/Fs(i)^2).^(-3/4).*fpl2/Fs(i)^2+(1-fpl2/Fs(i)^2).^(1/4)));
%     end
%     
% end
% hold on
% plot(Fd,phi3(1,: )./( phi1(1,:)-phi2(1,:)))