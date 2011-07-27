function [ d ] = try_recursive( n )
%UNTITLED33 Summary of this function goes here
%   Detailed explanation goes here

if n==1
    d=1;
else
    d=n+try_recursive(n-1);
end

end

