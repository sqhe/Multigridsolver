function [ va ] = bound_con( i, j, resol, ind2dm, dm2ind )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
h=1/resol;

ind = ij2ind( i, j, resol );
dm = ind2dm( ind );

x=i*h;
y=j*h;
x=x-0.5;
y=y-0.5;

t=5*pi/3;

r=0.3;

va=-(t*cos(t*(x-r))*sin(x)+sin(t*(x-r))*cos(x))*sin(y)*sin(t*(y-r))+(t*cos(t*(y-r))*sin(y)+sin(t*(y-r))*cos(y))*sin(x)*sin(t*(x-r));

end

