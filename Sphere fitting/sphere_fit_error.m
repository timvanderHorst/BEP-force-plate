function [ Error ] = sphere_fit_error( center, radius, data )
% This calculates the total error of a sphere fitting using the equation: 
% Error =sum((4*pi*R^2 - 4*pi*((xi - xc)^2 + (yi - yc)^2 + (zi - zc)^2))^2)
% center is the sphere's center, radius = sphere's radius, data = data
% points, error = total error
% http://pdfs.semanticscholar.org/ab00/f9883d9ae2addec8983d9c7044aa9fb1ae3d.pdf

% Input:
% center: 1 x 3 matrix containing sphere center [xc, yc, zc]
% radius: sphere radius
% data: M x 3 matrix containing all data points in cartesian coordinates
% Outputs:
% Error: Total error

Error = 0;
for i = 1 : size(data,1)
    Error = Error + (4*pi*radius^2 - 4*pi*(data(i,:)*data(i,:)'))^2;
end

