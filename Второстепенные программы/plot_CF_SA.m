% ���������� ������� �������������
% addpath('F:\������ ������������� ��� ������\�������\���������� ��� ks � kd')
% load('phi86_AX25000_sigma001_fs1HHz.mat')
close all
disp('��� ���������� ��������� ������� �������������  ������� 1' )
disp('��� ���������� ���������������� ������� �������������  ������� 2')
number=input('Number = ');
if number == 1
   disp('������� ��� ������ �������� x - ���������� ������� ������')
   disp('���������� �������� ')
   disp(x')
   x1 = input('x = ');
   kx=1;
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   disp('������� ��� ������ �������� y - ���������� ������� ������')
   disp('���������� �������� ')
   disp(y)
   y1 = input('y = ');
   ky=1;
   while y(ky)~=y1
    ky=ky+1;   
   end
   
  [X,Y]=meshgrid(fd/fs,z_all/1000);
  
  figure
  module = squeeze((abs(Coherency_Function(kx,ky,:,:))));
  angle1 = squeeze((angle(Coherency_Function(kx,ky,:,:))));
  surf(X,Y,module');
 % title('������ ��������� ������� ������������� SA') 
  ylabel('z [km]','FontSize',18);
  xlabel('f_d/f_s','FontSize',18);
  zlabel('|\Gamma|','FontSize',18)
  colormap ('jet')
  set(gca,'Fontsize',18)
  view(135,45)
  
  figure
   surf(X,Y,modunwrap(unwrap(angle1, [], 2)',3,2));
   %surf(X,Y,angle1');
   ylabel('z [km]','FontSize',18);
   xlabel(' f_d/f_s','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('���� ��������� ������� ������������� SA') 
   
elseif number==2
   disp('������� ��� ������ �������� ��������� ���������� ���������� ������� ������')
   disp('���������� �������� [���]')
   disp(fd/1000000)
   kd1=1;
   number1 = input('fd = ');
   while fd(kd1)/1000000~=number1
    kd1=kd1+1;   
   end
   
   disp('������� ��� ������ �������� y - ���������� ������� ������')
   disp('���������� �������� ')
   disp(y')
   y1 = input('y = ');
   ky=1;
   while y(ky)~=y1
    ky=ky+1;   
   end
   
   [X,Y]=meshgrid(x/1000,z_all/1000);
   
   figure
  module = squeeze((abs(Coherency_Function(:,ky,kd1,:))));
  angle1 = squeeze((angle(Coherency_Function(:,ky,kd1,:))));
  surf(X,Y,module');
   %title('������ ���������������� ������� ������������� SA') 
   ylabel('z [km]','FontSize',18);
   xlabel('x [km]','FontSize',18);
   zlabel('|\Gamma|','FontSize',18);
   colormap ('jet')
   
   set(gca,'Fontsize',20)
   view(135,45)
   
   figure
   surf(X,Y,modunwrap(unwrap(angle1,[],2)',3,2));
   ylabel('z [km]','FontSize',18);
   xlabel('x [km]','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('���� ���������������� ������� ������������� SA') 
   %%
   
   disp('������� ��� ������ �������� ��������� ���������� ���������� ������� ������')
   disp('���������� �������� [���]')
   disp(fd/1000000)
   kd1=1;
   number1 = input('fd = ');
   while fd(kd1)/1000000~=number1
    kd1=kd1+1;   
   end
   
   disp('������� ��� ������ �������� x - ���������� ������� ������')
   disp('���������� �������� ')
   disp(x')
   x1 = input('x = ');
   kx=1;
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   [X,Y]=meshgrid(y/1000,z_all/1000);
   
   figure
  module = squeeze((abs(Coherency_Function(kx,:,kd1,:))));
  angle1 = squeeze((angle(Coherency_Function(kx,:,kd1,:))));
  surf(X,Y,module');
   %title('������ ���������������� ������� ������������� SA') 
   ylabel('z [km]','FontSize',18);
   xlabel('y [km]','FontSize',18);
   zlabel('|\Gamma|','FontSize',18);
   colormap ('jet')
   
   set(gca,'Fontsize',20)
   view(135,45)
   
   figure
   surf(X,Y,modunwrap(unwrap(angle1,[],2)',3,2));
   ylabel('z [km]','FontSize',18);
   xlabel('y [km]','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('���� ���������������� ������� ������������� SA') 
   
end