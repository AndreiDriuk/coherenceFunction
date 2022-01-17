load('CohEarth.mat');
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
Fs = (0.6:0.1:1.5)*10^9;
Fd = ((-10:0.05:10)*10^8)'; %% t in DFT


fpl2 = [kplTop2, kplBot2, kplUnder2]*vel^2;
f = 10^9;

phi1 = zeros(length(Fs), length(Fd));
phi2 = zeros(length(Fs), length(Fd));
for i = 1:length(Fs)
    for j = 1:length(Fd)
        phi1(i,j) = 1i*(Fs(i)+Fd(j)/2)/vel*trapz(z, sqrt(sqrt(1-fpl2/(Fs(i)+Fd(j)/2)^2)));
        phi2(i,j) = 1i*(Fs(i)-Fd(j)/2)/vel*trapz(z, sqrt(sqrt(1-fpl2/(Fs(i)-Fd(j)/2)^2)));
    end
end

T0 = 10^-9;
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

 for ix1 = 1:FdL
     FFt(ix1,:) = exp(-2*pi*1i*i1(ix1)*j1/FdL);
 end


fs_d = 1/abs((Fd(length(Fd))-Fd(length(Fd)-1)));
L = length(Fd)/8;
t = 0:fs_d/L:3*2*fs_d-fs_d/L ;%-fs_d:fs_d/L:fs_d-fs_d/L;
MeanField = zeros(length(Fs), length(x),length(Fd), length(t));

for ix = 1:length(x)
%    for iy =1:length(y)
        for it = 1:length(t)
            MeanField(:, ix,:, it) = p02.*exp(phi1-phi2).*squeeze(CoherenceEarth(:, ix, 101,:)).*exp(-1i*Fd*t(it)).';  
        end
%    end
end
sum = zeros(1, length(Fs));

for ifs = 1:length(Fs)
    sum(ifs) = trapz(Fd, squeeze(MeanField(ifs,101,:,1)));
end


MeanFieldAfter = zeros(length(x), length(t));
for ix = 1:length(x)
    for it = 1:length(t)
        MeanFieldAfter(ix, it) = trapz(Fs,trapz(Fd, squeeze(MeanField(:, ix,:, it)),2));
    end
end


f = figure();
contour(x, t, abs(MeanFieldAfter)')
ylabel('td','FontSize',18);
xlabel('x','FontSize',18);
zlabel('We','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
