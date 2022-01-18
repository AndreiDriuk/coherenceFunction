%% TO DO Закончить расчеты функция рассеяния
clear
phi = 87.6274*pi/180;
AX = 9751;
BY = 145;
sigma2=0.1;
fs=10^9;

name = ['McDnld_phi' num2str(round(phi/pi*180)) '_AX' num2str(AX) '_BY' num2str(BY) '_sigma' num2str(sigma2) '_fs' num2str(fs/1000000000) 'HHz_y0.mat'];

load(name);

fs_d = 1/abs((fd(length(fd))-fd(length(fd)-1)));
fs_x = 1/abs((x(length(x))-x(length(x)-1)));
L = length(fd);
L_d = length(x);
t = (-fs_d/2:fs_d/L:fs_d/2-fs_d/L);
t_d = (-fs_x/2:fs_x/L_d:fs_x/2-fs_x/L_d);

iy = round(length(y)/2);

%fftCoherency = zeros(size(Coherency_Function));
fftCoherency = fftshift(fft2(squeeze(Coherency_Function(:, iy,: ,end))));

%kx = round(length(x)/2);
%ky = round(length(y)/2);

[X,Y]=meshgrid(x, fd);

figure(1)
module = squeeze((abs(Coherency_Function(:,iy,:,end))));
angle1 = squeeze((angle(Coherency_Function(:,iy,:,end))));
contour(x,fd,module');
% title('Модуль частотной функции когерентности SA') 
xlabel('x','FontSize',18);
ylabel('fd ','FontSize',18);
zlabel('|\Gamma|','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)

figure(2)
contour(X,Y,modunwrap(unwrap(angle1, [], 1), 1, 1)');
% title('Модуль частотной функции когерентности SA') 
ylabel('fd','FontSize',18);
xlabel('x м','FontSize',18);
zlabel('\angle \Gamma','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)


[X,Y]=meshgrid(t_d, t);

figure(3)
module = abs(fftCoherency);
angle1 = angle(fftCoherency);
contour(X,Y,((module')/max(max(module'))));
% title('Модуль частотной функции когерентности SA') 
ylabel('t','FontSize',18);
xlabel('kx','FontSize',18);
zlabel('|FFt[\Gamma]|','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)
xlim([-2e-3, 2e-3 ])
ylim([-1e-8, 1e-8 ])

figure(4)
contour(X,Y,modunwrap(unwrap(angle1, [], 1), 1, 1)');
% title('Модуль частотной функции когерентности SA') 
ylabel('t','FontSize',18);
xlabel('kx','FontSize',18);
zlabel('\angle \Gamma','FontSize',18)
colormap ('jet')
set(gca,'Fontsiz',18)
