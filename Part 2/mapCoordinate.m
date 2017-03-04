% This simple function maps an (x,y) coordinate to an index, where the
% mapped array consists of susequent rows (points with the same y value).

function index = mapCoordinate(x,y,xn)

index = (y-1).*xn + x;

end

