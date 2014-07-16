function [x1 y1] = projectionToLine(imlineVertices,pos)
% Calculate point in predefined line closest to arbitrary point
% Line:     y = m*x + b
% Point:    x0,y0
% Projection in line:   x1,y1

m = diff(imlineVertices(:,2))/diff(imlineVertices(:,1));
b = imlineVertices(1,2)-m*imlineVertices(1,1);
x0 = pos(:,1);
y0 = pos(:,2);

x1 = (m*y0+x0-m*b)/(m^2+1);
y1 = (m^2*y0+m*x0+b)/(m^2+1);

