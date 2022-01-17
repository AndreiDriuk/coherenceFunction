%addpath('F:\Данные когерентности для статьи\Экватор\Результаты для ks и kd')
%load('phi86_AX25000_sigma001_fs1HHz.mat')
close all
disp('Для построения частотной фукнции когерентности  введите 1' )
disp('Для построения пространственной функции когерентности  введите 2')
number=input('Number = ');
if number == 1
   
    disp('введите для какого значения y - переменной вывести данные')
    disp('Допустимые значения ')
    disp(y)
    y1 = input('y = ');
    ky=1;
    while y(ky)~=y1
        ky=ky+1;   
    end

    [X,Y]=meshgrid(fd/1000000,x/1000);
    
    figure
    module = squeeze((abs(Coherency_Function(:,ky,:,end))));
    angle1 = squeeze((angle(Coherency_Function(:,ky,:,end))));
    contourf(X,Y,module);
    % title('Модуль частотной функции когерентности SA') 
    ylabel('x [km]','FontSize',16);
    xlabel('f_d[MHz]','FontSize',16);
    %zlabel(' $|\Gamma| $','FontSize',16, 'Interpreter', 'latex');
    colormap ('jet')
    set(gca,'Fontsize',16)
   
    figure
    surf(X,Y,modunwrap(modunwrap(unwrap(angle1,[], 1 ), 1, 1), 1, 1));
    %[C, h] = contour(X,Y,unwrap(unwrap(angle1,[],2)), 'LineWidth', 2);
    %clabel(C,h, 'manual')
    grid on
    % title('Модуль частотной функции когерентности SA') 
    ylabel('x [km]','FontSize',16);
    xlabel('f_d[MHz]','FontSize',16);
    zlabel(' arg \Gamma ','FontSize',16);
    colormap ('jet')
    set(gca,'Fontsize',16) 

  
    disp('введите для какого значения x - переменной вывести данные')
    disp('Допустимые значения ')
    disp(x)
    x1 = input('x = ');
    kx=1;
    while x(kx)~=x1
        kx=kx+1;   
    end

    [X,Y]=meshgrid(fd/1000000,y/1000);
    
    figure
    module = squeeze((abs(Coherency_Function(kx,:,:,end))));
    angle1 = squeeze((angle(Coherency_Function(kx,:,:,end))));
    contourf(fd/1000000,y/1000,module);
    % title('Модуль частотной функции когерентности SA') 
    ylabel('y [km]','FontSize',16);
    xlabel('f_d[MHz]','FontSize',16);
    colormap ('jet')
    set(gca,'Fontsize',16)
    
    
    figure
    %[C, h] =contour(X,Y,unwrap(unwrap(angle1,[],2)),5, 'LineWidth', 2);
    %clabel(C,h, 'manual')
    grid on
    surf(X,Y,modunwrap(modunwrap(modunwrap(unwrap(angle1,[], 1 ), 1, 1), 1, 1), 1, 2));
    % title('Модуль частотной функции когерентности SA') 
    ylabel('y [km]','FontSize',16);
    xlabel('f_d[MHz]','FontSize',16);
    zlabel(' arg \Gamma ','FontSize',16);
    colormap ('jet')
    set(gca,'Fontsize',16)


elseif number==2
    disp('введите для какого значения частотной разностной переменной вывести данные')
    disp('Допустимые значения [МГЦ]')
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
    contourf(x/1000,y/1000,module');
    ylabel('y [km]','FontSize',16);
    xlabel('x [km]','FontSize',16);
    %zlabel(' $ |\Gamma| $','FontSize',16,'Interpreter','latex');
    colormap ('jet')
    set(gca,'Fontsize',16)

end
  
  