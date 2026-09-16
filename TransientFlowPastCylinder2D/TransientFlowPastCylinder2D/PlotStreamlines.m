function PlotStreamlines(Coords, elemCon, u, resolution)

% Coords  : (nNodes x 2) node coordinates
% elemCon : element connectivity
% u       : (nNodes x 2) velocity field
% resolution : grid resolution (e.g. 100)

if nargin < 4
    resolution = 100;
end

x = Coords(:,1);
y = Coords(:,2);
ux = u(:,1);
uy = u(:,2);

%% -------------------------------------------------
% Create structured grid
%% -------------------------------------------------

xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);

[X,Y] = meshgrid( linspace(xmin,xmax,resolution), ...
                  linspace(ymin,ymax,resolution) );

%% -------------------------------------------------
% Interpolate FEM velocity to grid
%% -------------------------------------------------

Ux = griddata(x,y,ux,X,Y,'cubic');
Uy = griddata(x,y,uy,X,Y,'cubic');

%% -------------------------------------------------
% Mask outside triangulation
%% -------------------------------------------------

TR = triangulation(elemCon,Coords);
ID = pointLocation(TR,[X(:) Y(:)]);

mask = reshape(~isnan(ID), size(X));
Ux(~mask) = NaN;
Uy(~mask) = NaN;

%% -------------------------------------------------
% Plot streamlines
%% -------------------------------------------------

figure
hold on
box on

% Optional: background velocity magnitude
Umag = sqrt(Ux.^2 + Uy.^2);
contourf(X,Y,Umag,20,'LineColor','none');
colormap jet
colorbar

% Streamlines
streamslice(X,Y,Ux,Uy);

axis equal
title('Streamlines')
xlabel('x')
ylabel('y')

end
