function [] = HW4

    h = 0.01; % Grid spacing

    % set a = 1
    w = @(k_x, k_y) sqrt( 1 - cos( (k_x + k_y) * 0.5 ) .* cos( (k_x - k_y) * 0.5 ) );

    % Brillouin Zone for 2D cubic lattice is squared bounded by [-pi / 2a, pi /
    % 2a]
    BrillouinZone = [-pi / 2 : h : pi / 2]';

    omega = zeros(numel(-pi / 2 : h : pi / 2), numel(-pi / 2 : h : pi / 2));


    for i = 1 : numel(BrillouinZone)
    % k_x loop

        omega(i, :) = w(BrillouinZone(i), BrillouinZone);

    end

    surf(BrillouinZone, BrillouinZone, omega);
    shading interp
    colormap(jet);
    xlabel('$$k_{x}\quad(\frac{1}{a})$$','Interpreter','latex');
    ylabel('$$k_{y}\quad(\frac{1}{a})$$','Interpreter','latex');
    zlabel('$$\frac{\omega (k_{x},k_{y})}{2\sqrt{\frac{C}{M}}}$$','Interpreter','latex');

    mat = zeros(numel(-pi / 2 : h : pi / 2) + 1, numel(-pi / 2 : h : pi / 2) + 1);
    mat(2 : end, 1) = BrillouinZone;
    mat(1, 2 : end) = BrillouinZone;
    mat(2 : end, 2 : end) = omega;
    ColorMosaicFunc(mat);
    xlabel('$$k_{x}\quad(\frac{1}{a})$$','Interpreter','latex');
    ylabel('$$k_{y}\quad(\frac{1}{a})$$','Interpreter','latex');
    colorbar
    
end

function [x,y,z] = ColorMosaicFunc( input )
% input = input data
% takes a 2d matrix, the stuff the XPlora spits out

    z = input;
    y = input(:, 1);
    x = input(1, :);

    x(1) = [];
    y(1) = [];
    z(:,1) = [];
    z(1, :) = [];

    figure;

    imagesc(x,y,z); 
    set(gca,'YDir','normal');
    set(gca, 'box', 'off')

    % remove bottom ticks
    %set(gca,'XtickLabel',[]);
    %set(gca,'xtick',[])

    colormap(jet);
    
end
