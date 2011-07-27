function [ x ] = pictr5( N )
%UNTITLED34 Summary of this function goes here
%   Detailed explanation goes here

h=1/N;

x=0:h:1;
y=ones(1,N+1);

for i=0:N
    figure(5);
    plot(x,h*i*y);
    hold on;
    figure(5);
    plot(h*i*y,x);
    hold on;
end

x=0.2:0.01:0.8;

y=0.5+sqrt(0.09-(x-0.5).^2);

figure(5);
plot(x,y);
hold on;

y=0.5-sqrt(0.09-(x-0.5).^2);

figure(5);
plot(x,y);
hold on;

end

