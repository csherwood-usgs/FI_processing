function s1 = s21m( s, zr, zo ) 
% s21m - Function to convert speed to reference elevation of 1 m
% vk = 0.41; % von Karmans constant...factors out, so omit
% us = (s*vk)./ log( zr ./zo );
% s1 = (us./vk).*log(1. ./zo );

us = s./ log( zr ./zo );
s1 = us.*log( 1. ./zo );