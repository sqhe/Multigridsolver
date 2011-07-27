function [ va ] = bd_con(  i, j, resol, ind2dm, dm2ind )
%UNTITLED35 Summary of this function goes here
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

%va=-2*(x^2+y^2-0.18);
%va=2*(x^2+y^2-x-y+2);
va=(x^2+4*x+1.91)*(y^2+4*y+1.91)*exp(x)*exp(y);
%va=-(x^2+4*x+1.91)*(y^2-0.09)*exp(x)*exp(y)+(y^2+4*y+1.91)*(x^2-0.09)*exp(x)*exp(y);
%va=10*pi/3*(cos(x)*cos(5*pi/3*(x-0.3))*sin(y)*sin(5*pi/3*(y-1))+cos(y)*cos(5*pi/3*(y-1))*sin(x)*sin(5*pi/3*(x-1)))-(50*pi^2/9+2)*sin(x)*sin(5*pi/3*(x-0.3))*sin(y)*sin(5*pi/3*(y-1));


%va=-((2*t*cos(t*(x-r))*cos(x)-(t^2+1)*sin(t*(x-r))*sin(x))*sin(y)*sin(t*(y-r))+(2*t*cos(t*(y-r))*cos(y)-(t^2+1)*sin(t*(y-r))*sin(y))*sin(x)*sin(t*(x-r)));

%va=-(((2*t*cos(t*(x-r))*cos(x)-(t^2+1)*sin(t*(x-r))*sin(x))+2*(t*cos(t*(x-r))*sin(x)+sin(t*(x-r))*cos(x)+sin(t*(x-r))*sin(x)))*exp(x)*exp(y)*sin(y)*sin(t*(y-r))+((2*t*cos(t*(y-r))*cos(y)-(t^2+1)*sin(t*(y-r))*sin(y))+2*(t*cos(t*(y-r))*sin(y)+sin(t*(y-r))*cos(y)+sin(t*(y-r))*sin(y)))*exp(x)*exp(y)*sin(x)*sin(t*(x-r)));

end

