%addpath('F:\������ ������������� ��� ������\�������\���������� ��� ks � kd')
%load('phi86_AX25000_sigma001_fs1HHz.mat')
close all
disp('��� ���������� ��������� ������� �������������  ������� 1' )
disp('��� ���������� ���������������� ������� �������������  ������� 2')
number=input('Number = ');
if number == 1
   
    disp('������� ��� ������ �������� y - ���������� ������� ������')
    disp('���������� �������� ')
    disp(y)
    y1 = input('y = ');
    ky=1;
    while y(ky)~=y1
        ky=ky+1;   
    end

    [X,Y]=meshgrid(fd/fs, x/1000);
    
%     figure
%     module = squeeze((abs(Coherency_Function(:,ky,:,end))));
%     angle1 = squeeze((angle(Coherency_Function(:,ky,:,end))));
%     surf(X,Y,module);
%   %  [C, h] = contour(X,Y,module,5, 'LineWidth', 2);%,'ShowText', 'on');%,'LabelSpacing', 200, 'TextStep', 100,'TextStepMode','manual')%'ShowText', 'on','LabelSpacing', 500);
%     grid on
%  %   clabel(C,h, 'manual');
%    % title('������ ��������� ������� ������������� SA') 
%     ylabel('x [km]','FontSize',16);
%     xlabel('f_d/f_s','FontSize',16);
%     zlabel(' $|\Gamma| $','FontSize',16, 'Interpreter', 'latex');
%     colormap ('jet')
%     set(gca,'Fontsize',16)
    
%     figure
%     surf(X,Y,unwrap(angle1,[],2));
%     %[C, h] = contour(X,Y,unwrap(unwrap(angle1,[],2)), 'LineWidth', 2);
%     %clabel(C,h, 'manual')
%     grid on
%     % title('������ ��������� ������� ������������� SA') 
%     ylabel('x [km]','FontSize',16);
%     xlabel('f_d/f_s','FontSize',16);
%     zlabel(' arg \Gamma ','FontSize',16);
%     colormap ('jet')
%     set(gca,'Fontsize',16)
    
    
%     figure
%     contourf(X,Y,module,10);
%     % title('������ ��������� ������� ������������� SA') 
%     ylabel('x [km]','FontSize',16);
%     xlabel('f_d [MHz]','FontSize',16);
%     zlabel('|�|','FontSize',16)
%     colormap ('jet')
%     set(gca,'Fontsize',10)
%     %view(135,45)
%     
%    figure
%    contourf(X,Y,unwrap(angle1,pi,2));
%    ylabel('x [km]','FontSize',20);
%    xlabel('f_d [MHz]','FontSize',20);
%    zlabel('angle \Gamma','FontSize',20);
%    colormap ('jet')
%    set(gca,'Fontsize',20)
%    view(135,45)
  
    disp('������� ��� ������ �������� x - ���������� ������� ������')
    disp('���������� �������� ')
    disp(x)
    x1 = input('x = ');
    kx=1;
    while x(kx)~=x1
        kx=kx+1;   
    end

    [X,Y]=meshgrid(fd/fs,y/1000);
    
    figure
    module = squeeze((abs(Coherency_Function(kx,:,:,end))));
    angle1 = squeeze((angle(Coherency_Function(kx,:,:,end))));
    %[C, h] = contour(X,Y,module,5, 'LineWidth', 2);
    %grid on
%    clabel(C,h, 'manual')
    surf(X,Y,module);
    % title('������ ��������� ������� ������������� SA')
    ylim([-0.3 0.3])
    ylabel('y [km]','FontSize',16);
    xlabel('f_d/f_s','FontSize',16);
    zlabel('  |\Gamma| ','FontSize',16);
    colormap ('jet')
    set(gca,'Fontsize',16)
    
    
    figure
    %[C, h] =contour(X,Y,unwrap(unwrap(angle1,[],2)),5, 'LineWidth', 2);
    %clabel(C,h, 'manual')
    grid on
    surf(X,Y,unwrap(unwrap(unwrap(angle1,[],1),[],2),[],2));
    ylim([-0.3 0.3])
    % title('������ ��������� ������� ������������� SA') 
    ylabel('y [km]','FontSize',16);
    xlabel('f_d/f_s','FontSize',16);
    zlabel(' arg \Gamma ','FontSize',16);
    colormap ('jet')
    set(gca,'Fontsize',16)

elseif number==2
    disp('������� ��� ������ �������� ��������� ���������� ���������� ������� ������')
    disp('���������� �������� [���]')
    disp(fd/1000000)
    kd1=1;
    number1 = input('fd = ');
    while fd(kd1)/1000000~=number1
        kd1=kd1+1;   
    end
    
    [X,Y]=meshgrid(x/1000,y/1000);
    module = squeeze((abs(Coherency_Function(:,:,kd1,end))));
    angle1 = squeeze((angle(Coherency_Function(:,:,kd1,end))));
    
    figure(1)
    %[C, h] = contour(X,Y,module',5, 'LineWidth', 2);
    %clabel(C,h, 'manual')
    grid on
    %surf(X,Y,module');
    contour(X,Y,module','ShowText', 'on');
    ylabel('y [km]','FontSize',16);
    xlabel('x [km]','FontSize',16);
    zlabel('  |\Gamma| ','FontSize',16);
    colormap ('jet')
    set(gca,'Fontsize',16)
%     figure(2)
%     module = squeeze((abs(Coherency_Function(:,:,kd1,end))));
%     angle1 = squeeze((angle(Coherency_Function(:,:,kd1,end))));
%     contourf(X,Y,module',10);
%     ylabel('y [km]','FontSize',16);
%     xlabel('x [km]','FontSize',16);
%     zlabel('arg \Gamma','FontSize',16, 'Interpreter', 'latex');
%     colormap ('jet')
%     set(gca,'Fontsize',20)
%     
   %view(135,45)
%     figure
%     contourf(X,Y,unwrap(angle1,[2*pi],1)');
%     ylabel('y [km]','FontSize',16);
%     xlabel('x [km]','FontSize',16);
%     zlabel('angle \Gamma','FontSize',16);
%     colormap ('jet')
%     set(gca,'Fontsize',20)
%     view(135,45)

end
  
  