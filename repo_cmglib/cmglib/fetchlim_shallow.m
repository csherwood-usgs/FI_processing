function [Hm,Tm]=fetchlim_shallow(U10,F,depth);
% FETCHLIM_SHALLOW - Calculates fetch-limited waves in shallow water
% using 1984 Shore Protection Manual eqn. 3-39
%
%  Usage:  [Hm,Tm] = fetchlim_shallow(Ua,F,depth);
%
%  Input:   U10 = wind at 10 m height (m/s) [vector]
%           F = fetch (m)
%           depth = depth (m)
%
%  Output:  Hm = sig. wave height (m) [vector]
%           Tm = Period (s) [vector]

%  Example: [Hm,Tm]=fetchlim_shallow(1:20,10000,5);
%              calculates [Hm,t] in 5 m water with 10 km F
%              for winds from 1 to 20 m/s
%  Note: Matches worked Example Problem 5 on p. 3-66 if clacs are
%        changed as show in code to work on Ua, the "wind-stress factor".
%        This routine assumes uncorrected winds at 10 m are input,
%        and calculates Ua internally.

%  csherwood@usgs.gov
%  Last revised 17-Feb-2012

if (nargin<3),
   help fetchlim_shallow
end

ua=0.71*(U10.^1.23);
%ua = Ua; % Using Ua directely matches Example Problem 5
uasq=ua.*ua;
g=9.81;

fact1=0.53*(((g*depth)./uasq).^0.75);
fact2=0.00565*(((g*F)./uasq).^0.50);
Hm=(uasq*(0.283).*tanh(fact1).*tanh((fact2./tanh(fact1))))/g;

fact3=0.833*(((g*depth)./uasq).^0.375);
fact4=0.0379*(((g*F)./uasq).^0.333);
Tm=(ua*7.54.*tanh(fact3).*tanh((fact4./tanh(fact3))))/g;

%      hrms=Hm/1.41;
