function w = prolate(n,isw)
% PROLATE  Prolate spheroidal wavefunction data window
%
% w = prolate(n,isw)
%
% Calculates prolate shperoidal wavefunction data window for
% high resolution fourier analysis.
% Reference: D.J.THOMSON, BELL SYST. TECH. J. 56,1769-1815 (1977)
%
% n is the number of points
% isw is a switch:
%    isw=4 : Return a 4-pi prolate window
%    isw=1 : Return the higher resoluton 1-pi prolate window
%
% Scale factors=integral(BOXCAR)/integral(PROLATE WINDOW) are:
%    4-pi prolate window--1.425658520238489
%    1-pi prolate window--1.057568010371401
% These are the numbers to multiply the spectrum by for comparison
% with other windows.

% Fortran source code from Adam Schulz
% Converted to MATLAB function by Chris Sherwood, 5/9/89

% Chris Sherwood, USGS
% March 17, 1999

     a = ones(n,1);
     w = zeros(size(a));
     c1=2.6197747176990866d-11;
     c2=2.9812025862125737d-10;
     c3=3.0793023552299688d-9;
     c4=2.8727486379692354d-8;
     c5=2.4073904863499725d-7;
     c6=1.8011359410323110d-6;
     c7=1.1948784162527709d-5;
     c8=6.9746276641509466d-5;
     c9=3.5507361197109845d-4;
     c10=1.5607376779150113d-3;
     c11=5.8542015072142441d-3;
     c12=1.8482388295519675d-2;
     c13=4.8315671140720506d-2;
     c14=1.0252816895203814d-1;
     c15=1.7233583271499150d-1;
     c16=2.2242525852102708d-1;
     c17=2.1163435697968192d-1;
     c18=1.4041394473085307d-1;
     c19=5.9923940532892353d-2;
     c20=1.4476509897632850d-2;
     c21=1.5672417352380246d-3;
     c22=4.2904633140034110d-5;
     if(isw == 4),
        d=sqrt(2.d0/.508125548147497d0);
        for i=1:n,
          x=(2*i-1)/n-1;
          u=(1.d0-x)*(1.d0+x);
          d00=d*(((((((((((((((((((((c1*u+c2)*u+c3)*u+c4)*u+c5)*u+c6)*u+...
          c7)*u+c8)*u+c9)*u+c10)*u+c11)*u+c12)*u+c13)*u+c14)*u+c15)*u+c16)...
          *u+c17)*u+c18)*u+c19)*u+c20)*u+c21)*u+c22);
          w(i) = d00*a(i);
        end
        return
      elseif (isw ==1),
        a1=5.3476939016920851d-11;
        a2=2.2654256220146656d-9;
        a3=7.8075102004229667d-8;
        a4=2.1373409644281953d-6;
        a5=4.5094847544714943d-5;
        a6=7.0498957221483167d-4;
        a7=7.7412693304064753d-3;
        a8=5.5280627452077586d-2;
        a9=2.2753754228751827d-1;
        a10=4.3433904277546202d-1;
        a11=2.2902051859068017d-1;
        d=sqrt(2.d0);
        for i=1:n,
          x=(2*i-1)/n-1;
          u=(1.d0-x)*(1.d0+x);
          d00=d*((((((((((a1*u+a2)*u+a3)*u+a4)*u+a5)*u+a6)*u+a7)*u+a8)*u+...
          a9)*u+a10)*u+a11);
          w(i)=d00*a(i);
        end
        return
      else
        error('isw must be 1 or 4')
        return
      end












