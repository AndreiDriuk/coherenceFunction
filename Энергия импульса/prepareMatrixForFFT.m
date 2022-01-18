%% Скрипт создает входной файл в котором Функция когерентности представлена на Земле как функция Fs, Fd, x координат
ztop =[-100000, 0];%:10000:0;
lt = 140660; % coefficient in cosh(z/lt) 
lb = 198722; % coefficient in 1-z^2/lb^2
lend = 532000; % from bottom layer to the point of view
kplm =  0.2487; % from NeQuick model
phi = 87.6274*pi/180; % angle between trace
zbot = [0, lb];%[0:10000:lb,lb];
zunder = [lb, lend];%[lb:50000:lend,lend];
z_all = [ztop,zbot,zunder];


fd = ([-10:0.002:10]*10^8)';
sigma2=0.1;


v = 299792458; %
%fplm = 10^7; %
%kplm=0.215335266700290; %
%ks = fs*2*pi/v;%kz(for )
kd = fd*2*pi/v;
AX = 9751;
BY = 145;
%CZ = 76;

ax = sqrt(AX^2*BY^2/((AX*cos(phi))^2+(BY*sin(phi))^2));
by = BY;
le = sqrt(AX^2*BY^2/((AX*sin(phi))^2+(BY*cos(phi))^2));
fs1 = [9.8:0.005:10.2]*10^8;
CoherenceEarth = zeros(length(fs1), 201, length(fd));

for ifs = 1:length(fs1)
    disp(length(fs1))
    disp(ifs)
    name = ['McDnld_phi' num2str(round(phi/pi*180)) '_AX' num2str(AX) '_BY' num2str(BY) '_sigma' num2str(sigma2) '_fs' num2str(fs1(ifs)/1000000000) 'HHz_y0.mat'];
    load(name )
    CoherenceEarth(ifs, :,:) = squeeze(Coherency_Function(:,:,:, end));
end

save CohEarth.mat CoherenceEarth fd z_all x y ax by le -v7.3