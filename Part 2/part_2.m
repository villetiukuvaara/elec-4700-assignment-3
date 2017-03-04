%% Part 2 (part_2.m)
% This is basically the same as the last part of assignment 2. The equation
% used to produce the G matrix is the following:
% 
% $$\nabla(\sigma \nabla V) = 0$$
%
% $$\frac{\partial \sigma}{\partial x}\frac{\partial V}{\partial x}
% + \sigma\frac{\partial^2 V}{\partial x^2}
% + \frac{\partial \sigma}{\partial y}\frac{\partial V}{\partial y}
% + \sigma\frac{\partial^2 V}{\partial y^2} = 0$$
%
% To use finite differences (FD), the derivatives can be converted to
% discrete numerical approximations as follows. Note that there is more
% than one possible solution here, since there is more than one
% approximation of the first and second derivatives.
%
% $$ \left( \frac{\sigma_{x+1,y}-\sigma_{x-1,y}}{2\Delta x} \right)
% \left( \frac{V_{x+1,y}-V_{x-1,y}}{2\Delta x} \right)
% + \sigma_{x,y}  \left( \frac{V_{x+1,y} - 2V_{x,y}+V_{x-1,y}}{(\Delta
% x)^2} \right)$$
%
% $$+ \left( \frac{\sigma_{x,y+1}-\sigma_{x,y-1}}{2\Delta y} \right)
% \left( \frac{V_{x,y+1}-V_{x,y-1}}{2\Delta y} \right)
% + \sigma_{x,y}  \left( \frac{V_{x,y+1} - 2V_{x,y}+V_{x,y-1}}{(\Delta y)^2} \right)=0 $$
%
% Finally, to construct the G matrix, it is better to express the above
% expression with coefficients in front of the $V$ terms. The final
% equation is:
%
% $$\left(\frac{1}{4(\Delta x)^2}(\sigma_{x+1,y}-\sigma_{x-1,y}) +
% \frac{\sigma_{x,y}}{(\Delta x)^2}\right) V_{x+1,y}
% +\left(\frac{-1}{4(\Delta x)^2}(\sigma_{x+1,y}-\sigma_{x-1,y}) +
% \frac{\sigma_{x,y}}{(\Delta x)^2}\right) V_{x-1,y}$$
%
% $$+\left(\frac{1}{4(\Delta y)^2}(\sigma_{x,y+1}-\sigma_{x,y-1}) +
% \frac{\sigma_{x,y}}{(\Delta y)^2}\right) V_{x,y+1}
% +\left(\frac{-1}{4(\Delta y)^2}(\sigma_{x,y+1}-\sigma_{x,y-1}) +
% \frac{\sigma_{x,y}}{(\Delta y)^2}\right) V_{x,y-1}$$
%
% $$-2\sigma_{x,y} \left( \frac{1}{(\Delta x)^2}
% + \frac{1}{(\Delta y)^2}\right) V_{x,y} =0$$
%
% Here are some parameters for the simulation. Note the the simulation is
% performed with a scale factor of 1 : 100 nm. This is better to avoid
% numerical issues.

clear all;
close all;

W = 1;
L = 2;
scale = 100e-9; % Scale factor, so simulation doesn't happen at very small numbers
V0 = 1;

dx = 0.025; % Mesh spacing along x
dy = 0.025; % Mesh spacing along y
nx = ceil(L/dx); % Number of points along x
ny = ceil(W/dy); % Number of points along y
dx = L/nx;
dy = W/ny;

%%
% The finite differences method can be implemented using a matrix
% calculation, where $GV=F$. $V$ is the voltages at the discrete points,
% $F$ is the "forcing" matrix that is used to set the boundary conditions,
% and G defines how the voltages are each point are related to the other
% voltages. The G matrix is defined by the equation derived above. First,
% I will construct a matrix C that contains the conductivities at all points.
% There are two rectangular highy resistive areas, and the rest of the area
% has a nominal conductivity of 1.

Lb = 0.4; % Length of the highly resistive regions
Wb = 0.4; % Width of the highly resistive regions
sigma1 = 1; % Nominal conductivity
sigma2 = 1e-2; % Conductivity in highly resistive regions

% Construct the C matrix:
C = sigma1.*ones(ny,nx);
Csubtract = zeros(ny,nx);

for x=1:nx
    for y=1:ny
        xx = x*dx;
        yy = y*dy;
        
        % The resistivity is made high in the rectangular regions:
        if(xx <= (L+Lb)/2 && xx >= (L-Lb)/2 && (yy >= W-Wb || yy <= Wb))
            Csubtract(y,x) = sigma1-sigma2;
        end
    end
end

% Filter (smooth) the condicivity a bit to avoid numerical issues that can
% occur if the derivatives are very large.
Csubtract = imgaussfilt(Csubtract, 1);
C = C - Csubtract;

%%
% Below, the conducitivity is plotted. I performed some filtering/smoothing
% so that the derivative is not very large (approaching infinity) at the
% junction of the two regions. It also also more realistic. Often,
% physical junctions are gradients.

figure(1);
surf(linspace(0,L.*scale,nx),linspace(0,W.*scale,ny),C);
title('Conductivity');
view(30,45);
xlabel('x (m)');
ylabel('y (m)');
grid on;

%%
% Here, the G matrix is generated. The mapCoordinate function takes a
% discrete (x,y) coordinate and converts it to an index in the V array.
% This mapping is performed such that the V array consists of successive
% rows (points with the same y-value).

G = zeros(nx*ny,nx*ny);
F = zeros(nx*ny,1);

dx2 = 1./(dx.^2);
dy2 = 1./(dy.^2);

for x=2:(nx-1)
    for y=2:(ny-1)
        index = mapCoordinate(x,y,nx);
        
        % Apply the equation derived earlier:
        G(index,index) = -2.*C(y,x).*(dx2 + dy2);
        G(index, mapCoordinate(x+1,y,nx)) = dx2.*(0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        G(index, mapCoordinate(x-1,y,nx)) = dx2.*(-0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        
        G(index, mapCoordinate(x,y+1,nx)) = dy2.*(0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
        G(index, mapCoordinate(x,y-1,nx)) = dy2.*(-0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
    end
end

%%
% Next, the F matrix is generated.

% The part of the boundary at V = V0
for x=1:floor(nx/2)
    index = mapCoordinate(x,1,nx);
    G(index,index) = 1;
    F(index) = V0;
    
    index = mapCoordinate(x,ny,nx);
    G(index,index) = 1;
    F(index) = V0;
end

% The part of the boundary where V = 0
for x=floor(nx/2):nx
    index = mapCoordinate(x,1,nx);
    G(index,index) = 1;
    F(index) = 0;
    
    index = mapCoordinate(x,ny,nx);
    G(index,index) = 1;
    F(index) = 0;
end

% The vertical boundaries
for y=1:ny
    index = mapCoordinate(1,y,nx);
    G(index,index) = 1;
    F(index) = V0;
    
    index = mapCoordinate(nx,y,nx);
    G(index,index) = 1;
    F(index) = 0;
end

%%
% After setting up the matrices, getting the solution is trivial. To
% convert back to a matrix where the voltages can be accessed with (x,y)
% coordinates, the reshape function is used.

V = G\F;
V = reshape(V,[],ny)';

%%
% The length and width were $L = 2$ and $W = 1$, which need to be scaled to
% 200 nm and 100 nm, respectively

figure(2);
surf(linspace(0,L.*scale,nx),linspace(0,W.*scale,ny),V);
view(30,45);
xlabel('x (m)');
ylabel('y (m)');
title('Electric Potential (V)');
grid on;

%%
% The electric field is $E=-\nabla V$. Here it is plotted along with the
% voltage.

figure(3);
[Ex,Ey] = gradient(V,dx.*scale,dy.*scale);
Ex = -1.*Ex;
Ey = -1.*Ey;
%contour(linspace(0,L.*scale,nx),linspace(0,W.*scale,ny),V);
%hold on;
quiver(linspace(0,L.*scale,nx),linspace(0,W.*scale,ny),Ex,Ey,4);
xlabel('x (m)');
ylabel('y (m)');
title('Electric Field (V/m)');
axis([0 L.*scale 0 W.*scale]);
grid on;