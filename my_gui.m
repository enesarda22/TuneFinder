% digital art
clearvars *
close all
clc
figure
    t = linspace(0,4*pi,500);
    x = cos(t);
    y = sin(t);
    z = (cos(2.*t)).*sin(20.*t);
    a = -pi : pi/2 : pi;                                % Define Corners
    ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
    A = [cos(a+ph); cos(a+ph)]/cos(ph)./5;
    B = [sin(a+ph); sin(a+ph)]/sin(ph)./5;
    C = [-ones(size(a)); ones(size(a))]./5;

    colormap bone
    cmap = colormap;
    temp = cmap(1:2:end,:);
    temp_flip = flip(temp);
    cmap = [temp; temp_flip];
    k=1;
            t_k = t(k);
            x_k = x(k);
            y_k = y(k);
            z_k = z(k);
            Plot_color=cmap(round(k*end/(length(t)+5))+1,:);

            plot3(x_k, y_k, z_k, 'o', 'LineWidth', 3, 'MarkerSize', 15,'Color','r')
            hold on
            surf(A,B,C, 'FaceColor',Plot_color,'EdgeColor','none')                      % Plot Cube
            patch(A', B', C', Plot_color, 'EdgeColor', 'none') 
            plot3(x,y,z,'-', 'LineWidth', 3,'Color',Plot_color)

            hold off
            view([30 + 360*t_k/(4*pi) 50+35*cos(pi/2+t_k/4)])
            axis off
            title('TundeFinder')
            drawnow

set(figure(1),'WindowButtonDownFcn',@mytestcallback);



function []=mytestcallback(src,~)
    t = linspace(0,4*pi,500);
    x = cos(t);
    y = sin(t);
    z = (cos(2.*t)).*sin(20.*t);
    a = -pi : pi/2 : pi;                                % Define Corners
    ph = pi/4;                                          % Define Angular Orientation (‘Phase’)
    A = [cos(a+ph); cos(a+ph)]/cos(ph)./5;
    B = [sin(a+ph); sin(a+ph)]/sin(ph)./5;
    C = [-ones(size(a)); ones(size(a))]./5;

    colormap bone
    cmap = colormap;
    temp = cmap(1:2:end,:);
    temp_flip = flip(temp);
    cmap = [temp; temp_flip];
    while 1
    
        for k=1:length(t)
            t_k = t(k);
            x_k = x(k);
            y_k = y(k);
            z_k = z(k);
            Plot_color=cmap(round(k*end/(length(t)+5))+1,:);

            plot3(x_k, y_k, z_k, 'o', 'LineWidth', 3, 'MarkerSize', 15,'Color','r')
            hold on
            surf(A,B,C, 'FaceColor',Plot_color,'EdgeColor','none')                      % Plot Cube
            patch(A', B', C', Plot_color, 'EdgeColor', 'none') 
            plot3(x,y,z,'-', 'LineWidth', 3,'Color',Plot_color)

            hold off
            view([30 + 360*t_k/(4*pi) 50+35*cos(pi/2+t_k/4)])
            axis off
            drawnow
        end
    end
end