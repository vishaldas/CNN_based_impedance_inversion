function [ seed ] = getseed(  )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
serieslength =1;
seedmin = 24651;
seedmax = 64656655;
halfrange = (seedmax - seedmin) / 2; %replace max and min by actual values, values must be odd.
seed = (randi(halfrange, 1, serieslength) - 1) * 2 + seedmin;

end

