% This function implements the finite difference (FD) algorithm for the
% case presented in Problem 2. The algorithm is also in part_2.m, where it
% is explained. It is not thoroughly commented again here.
function [C, V, Ex, Ey, Jx, Jy] = nonuniformConductivityFD(L, W, nx, ny, Wb, Lb, sigma1, sigma2)

V0 = 1;
dx = L/nx; % Number of points along x
dy = W/ny; % Number of points along y

% Construct the C matrix:
C = sigma1.*ones(ny,nx);
if(W > 0 && L > 0)
    Csubtract = zeros(ny,nx);
    
    for x=1:nx
        for y=1:ny
            xx = x*dx;
            yy = y*dy;
            
            % The resistivity is made high in the rectangular regions:
            if(xx < (L+Lb)/2 && xx > (L-Lb)/2 && (yy > W-Wb || yy < Wb))
                Csubtract(y,x) = sigma1-sigma2;
            end
        end
    end
    
    % Filter (smooth) the condicivity a bit to avoid numerical issues that can
    % occur if the derivatives are very large.
    Csubtract = imgaussfilt(Csubtract, 1);
    C = C - Csubtract;
end

G = zeros(nx*ny,nx*ny);
F = zeros(nx*ny,1);

dx2 = 1./(dx.^2);
dy2 = 1./(dy.^2);

for x=2:(nx-1)
    for y=2:(ny-1)
        index = mapCoordinate(x,y,nx);
        
        G(index,index) = -2.*C(y,x).*(dx2 + dy2);
        G(index, mapCoordinate(x+1,y,nx)) = dx2.*(0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        G(index, mapCoordinate(x-1,y,nx)) = dx2.*(-0.25.*(C(y,x+1) - C(y,x-1)) + C(y,x));
        
        G(index, mapCoordinate(x,y+1,nx)) = dy2.*(0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
        G(index, mapCoordinate(x,y-1,nx)) = dy2.*(-0.25.*(C(y+1,x) - C(y-1,x)) + C(y,x));
    end
end

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

for y=1:ny
    index = mapCoordinate(1,y,nx);
    G(index,index) = 1;
    F(index) = V0;
    
    index = mapCoordinate(nx,y,nx);
    G(index,index) = 1;
    F(index) = 0;
end

V = G\F;
V = reshape(V,[],ny)';

[Ex,Ey] = gradient(V,dx,dy);
Ex = -1.*Ex;
Ey = -1.*Ey;

Jx = C.*Ex;
Jy = C.*Ey;

end