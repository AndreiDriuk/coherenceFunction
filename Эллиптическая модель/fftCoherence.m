clear
phi = 87.6274*pi/180;
AX = 9751;
BY = 145;
sigma2=0.1;
fs=0.9*10^9;

name = ['McDnld_phi' num2str(round(phi/pi*180)) '_AX' num2str(AX) '_BY' num2str(BY) '_sigma' num2str(sigma2) '_fs' num2str(fs/1000000000) 'HHz_y0.mat'];

load(name);

fs_d = 1/abs((fd(length(fd))-fd(length(fd)-1)));
L = length(fd);
t = (-fs_d/2:fs_d/L:fs_d/2-fs_d/L);

fftCoherency = zeros(size(Coherency_Function));
%for ix =1:length(x)
    %for iy = 1:length(y)
        %for iz = 1:length(z_all)
                fftCoherency(ix, iy, :, iz) = fftshift(fft2(squeeze(Coherency_Function(:, iy,: ,iz))));
        %end
    %end
%end

kx = round(length(x)/2);
ky = round(length(y)/2);

[X,Y]=meshgrid(t,z_all/1000);

figure(1)
module = squeeze((abs(Coherency_Function(1,ky,:,:))));
angle1 = squeeze((angle(Coherency_Function(1,ky,:,:))));
surf(X,Y,module');
% title('Модуль частотной функции когерентности SA') 
ylabel('z [km]','FontSize',18);
xlabel('f_d/f_s','FontSize',18);
zlabel('|\Gamma|','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)

figure(2)
surf(X,Y,angle1');
% title('Модуль частотной функции когерентности SA') 
ylabel('z [km]','FontSize',18);
xlabel('f_d/f_s','FontSize',18);
zlabel('\angle \Gamma','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)
[X,Y]=meshgrid(t,z_all/1000);

module = squeeze((abs(fftCoherency(kx,ky,:,:))));
angle1 = squeeze((angle(fftCoherency(kx,ky,:,:))));
figure(3)
surf(X,Y,module');
 % title('Модуль частотной функции когерентности SA') 
 ylabel('z [km]','FontSize',18);
 xlabel('t [µs]','FontSize',18);
 zlabel('| F[\Gamma]| ','FontSize',18)
 colormap ('jet')
 set(gca,'Fontsize',18)
 view(135,45)
 xlim([-2*10^(-8), 2*10^-8])
 
 
figure(4)
surf(X,Y,unwrap(angle1,[],1)');
ylabel('z [km]','FontSize',18);
xlabel('t [µs]','FontSize',18);
zlabel('\angle F[\Gamma]','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)
 

 module1 = squeeze((abs(fftCoherency(:,ky,:,end))));
 angle2 = squeeze((angle(fftCoherency(:,ky,:,end))));
 [X,Y]=meshgrid(x ,t);
 
figure(5)
surf(X,Y,module1');
 % title('Модуль частотной функции когерентности SA') 
ylabel('t [µs]','FontSize',18);
xlabel('x [m]','FontSize',18);
 zlabel('| F[\Gamma]|','FontSize',18)
 colormap ('jet')
 set(gca,'Fontsize',18)
 view(135,45)

 
 
figure(6)
surf(X,Y,unwrap(unwrap(angle2,[],2),[],1)');
ylabel('t [µs]','FontSize',18);
xlabel('x [m]','FontSize',18);
zlabel('\angle F[\Gamma]','FontSize',18)
colormap ('jet')
set(gca,'Fontsize',18)
view(135,45)