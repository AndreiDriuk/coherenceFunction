clear 
ztop = -1000000:10000:0;
lt = 154154;
lb = 147536;
lend = 422000;
kplm = 0.22417;
phi = pi/2;%86.0433/180*pi;
zbot = [0:10000:lb,lb];
zunder = [lb:10000:lend,lend];
z_all = [ztop,zbot,zunder];

fs=10^9;
fd = 0*10^8;
sigma2 = (0.1)^2;
p = 11/3;

v = 299792458; %скорость распространения
%fplm = 10^7; % максимальная плазменная частота
%kplm=0.215335266700290; % максимальная плазменное волновое число
ks = fs*2*pi/v;
kd = fd*2*pi/v;
AX = 25000;
BY = 3000;
CZ = 3000;
AXP = 517;%3630;
BYP = 61;%433;
CZP = BYP;%3000;
ax = sqrt(AX^2*BY^2/((AX*cos(phi))^2+(BY*sin(phi))^2));
by = BY;
le = sqrt(AX^2*BY^2/((AX*sin(phi))^2+(BY*cos(phi))^2));

axp = sqrt(AX^2*BY^2/((AX*cos(phi))^2+(BY*sin(phi))^2));
byp = BY;
lep = sqrt(AX^2*BY^2/((AX*sin(phi))^2+(BY*cos(phi))^2));

x = -5000:5:5000;
y = -500:1:500;

coeff = 1/4*sigma2*le/(ks^2-1/4*kd^2).^2;
coeffp = 1/4*sigma2*lep/(ks^2-1/4*kd^2).^2;

%Aex = 1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/2^(-(p-2)/2)*(zMcDonald((p-2)/2,2*pi*x/ax));
%Aey = 1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/2^(-(p-2)/2)*(zMcDonald((p-2)/2,2*pi*y/by)); 
Aem = zeros(length(x), length(y));
Aep = zeros(length(x), length(y));
for i = 1:length(x)
    Aem(i,:) = 1/(pi^(1/2))/gamma(p-2)*gamma(p/2-1/2)/(2^(-(p-2)/2))*(zMcDonald((p-2)/2,2*pi*sqrt(x(i)^2/ax^2+y.^2/by^2)));
    Aep(i,:) = (1-x(i)^2/axp^2-y.^2/byp^2);
end

%DexMc = 1-Aex/max(Aex);
%DeyMc = 1-Aey/max(Aey);
DeMc = coeff*(1-Aem);
DeP = coeffp*(1-Aep);

f_top = lt*(tanh(ztop/lt)+1);
f_bot = (lt+zbot-zbot.^3/(3*lb^2));
f_under = ones(1, length(zunder)).*(lt+2/3*lb);

Coherence_FunctionMc = zeros(length(x), length(y), length(z_all));
Coherence_FunctionPar = zeros(length(x), length(y), length(z_all));
for ix = 1:length(x)
    for  iy = 1:length(y)
        %Coherence_FunctionMc(ix, iy, :) = exp(-coeff*(DexMc(ix)*DeyMc(iy))*[f_top, f_bot, f_under]);
        Coherence_FunctionMc(ix, iy, :) = exp(-(DeMc(ix, iy))*[f_top, f_bot, f_under]);
        Coherence_FunctionPar(ix, iy, :) = exp(-DeP(ix, iy)*[f_top, f_bot, f_under]);
    end
end
figure(3)
plot(y, DeMc(1001, :),'LineWidth', 2)
hold on
plot(y, DeP(1001,:),'LineWidth', 2)
grid minor
legend('McDonald model','parabolic model')
set(gca,'FontSize', 14)

figure(4)
plot(x, DeMc(:, 501),'LineWidth', 2)
hold on
plot(x, DeP(:,501),'LineWidth', 2)
grid minor
legend('McDonald model','parabolic model')
set(gca,'FontSize', 14)

figure(1)
plot(x, squeeze(Coherence_FunctionMc(:,601,end)),'LineWidth', 2)
hold on
plot(x, squeeze(Coherence_FunctionPar(:,601,end)),'r','LineWidth', 2)
grid minor
legend('McDonald model','parabolic model')
set(gca,'FontSize', 14)
%xlim([-500,500])
xlabel('x m', 'FontSize', 14)
hold off

figure(2)
plot(y, squeeze(Coherence_FunctionMc(1101,:,end)),'LineWidth', 2)
hold on
plot(y, squeeze(Coherence_FunctionPar(1101,:,end)),'r','LineWidth', 2)
grid minor
legend('McDonald model','parabolic model')
set(gca,'FontSize', 14)
%xlim([-200, 200])
xlabel('y m', 'FontSize', 14)
hold off
