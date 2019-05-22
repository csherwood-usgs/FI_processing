function [hs,t]=fetch_limited(depth,fetch,u10);
% FETCH_LIMITED Calculates fetch-limited waves in shallow water
% using 1984 Shore Protection Manual eqn. 3-39 
%
%  Usage:  [hs,t] = fetch_limited(depth,fetch,u10);
% 
%  Input:   depth = depth (m)
%           fetch = fetch (m)
%             u10 = wind at 10 m height (m/s) [vector]
%
%  Output:     hs = sig. wave height (m) [vector]
%               t = Period (s) [vector]

%  Example:   [hs,t]=fetch_limited(5,10000,1:20);
%              calculates [hs,t] in 5 m water with 10 km fetch for winds
%                                from 1 to 20 m/s

      if (nargin<3),
        help fetch_limited
      end         
      ua=0.71*(u10.^1.23);
      uasq=ua.*ua;
      g=9.81;

      fact1=0.53*(((g*depth)./uasq).^0.75);
      fact2=0.00565*(((g*fetch)./uasq).^0.50);
      hs=(uasq*(0.283).*tanh(fact1).*tanh((fact2./tanh(fact1))))/g;

      fact3=0.833*(((g*depth)./uasq).^0.375);
      fact4=0.0379*(((g*fetch)./uasq).^0.333);
      t=(ua*7.54.*tanh(fact3).*tanh((fact4./tanh(fact3))))/g;

%      hrms=hs/1.41;
