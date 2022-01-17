% ���������� ������� �������������
% addpath('F:\������ ������������� ��� ������\�������\���������� ��� ks � kd')
% load('phi86_AX25000_sigma001_fs1HHz.mat')

disp('��� ���������� ��������� ������� �������������  ������� 1' )
disp('��� ���������� ���������������� ������� �������������  ������� 2')
number=input('Number = ');
if number == 1
   disp('������� ��� ������ �������� r - ���������� ������� ������')
   disp('���������� �������� ')
   disp(x')
   x1 = input('r = ');
   kx=1;
   while x(kx)~=x1
    kx=kx+1;   
   end
   
   
  [X,Y]=meshgrid(fd/1000000,z_all/1000);
  
  figure
  module = squeeze((abs(Coherency_Function(kx,:,:))));
  angle1 = squeeze((angle(Coherency_Function(kx,:,:))));
  surf(X,Y,module');
 % title('������ ��������� ������� ������������� SA') 
  ylabel('z [km]','FontSize',18);
  xlabel('f_d [MHz]','FontSize',18);
  zlabel('|\Gamma|','FontSize',18)
  colormap ('jet')
  set(gca,'Fontsize',18)
  view(135,45)
  
  figure
   surf(X,Y,modunwrap(unwrap(angle1, [], 2)',3,2));
   %surf(X,Y,angle1');
   ylabel('z [km]','FontSize',18);
   xlabel(' f_d  [MHz]','FontSize',18);
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
   
   [X,Y]=meshgrid(x/1000,z_all/1000);
   
   figure
  module = squeeze((abs(Coherency_Function(:,kd1,:))));
  angle1 = squeeze((angle(Coherency_Function(:,kd1,:))));
  surf(X,Y,module');
   %title('������ ���������������� ������� ������������� SA') 
   ylabel('z [km]','FontSize',18);
   xlabel('\rho [km]','FontSize',18);
   zlabel('|\Gamma|','FontSize',18);
   colormap ('jet')
   
   set(gca,'Fontsize',20)
   view(135,45)
   
   figure
   surf(X,Y,modunwrap(unwrap(angle1,[],2)',3,2));
   ylabel('z [km]','FontSize',18);
   xlabel('\rho [km]','FontSize',18);
   zlabel('\angle \Gamma','FontSize',18);
   colormap ('jet')
   set(gca,'Fontsize',18)
   view(135,45)
   %title('���� ���������������� ������� ������������� SA') 
   %%
   
end