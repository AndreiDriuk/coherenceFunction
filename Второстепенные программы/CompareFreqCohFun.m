clear

addpath('F:\Данные когерентности для статьи\Экватор\Результаты New')
load('phi86_AX50000_sigma0.1_fs1HHz.mat')

figure
plot(fd,abs(squeeze(FUnderFunSA(1,1,:,end))),'LineWidth', 2)
hold on
plot(fd, abs(squeeze(EXPUnderSA(1,1,:,end))),'r','LineWidth', 2)
grid minor
xlabel('z [km]','FontSize',16);
set(gca,'FontSize', 16)