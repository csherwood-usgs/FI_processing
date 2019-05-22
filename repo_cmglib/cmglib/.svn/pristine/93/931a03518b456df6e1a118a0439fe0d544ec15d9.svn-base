function [rheight, rlength ] = hw_ripple(Taucwmax,Taucr,d50,Ab)
% hw_ripple - Harris-Wiberg ripple calculations from ROMS translated to MATLAB
% function [rheight, rlength ] = hw_ripple(Taucwmax,Taucr,d50,Ab)
% Units are in MKS

% csherwood@usgs.gov
% 14 Aug 2008
% Coefficients for Wiberg-Harris ripple predictor
coef_a1=0.095;
coef_a2=0.442;
coef_a3=2.280;

tstar=Taucwmax ./(Taucr+eps);
if(tstar<=1.0),
   % defaults when subcritical stresses:
   rheight = 0.015;
   rlength = rheight/0.012;
   return
end
if(tstar>1),
   %
   %  Threshold of motion exceeded - calculate new zoST and zoBF
   %  Calculate saltation roughness according to Wiberg & Rubin (1989)
   %  (Eqn. 11 in Harris & Wiberg, 2001)
   %  (d50 is in m, but this formula needs cm)
   %
   coef_st=0.0204*log(100.0*d50)^2+0.0220*log(100.0*d50)+0.0709;
   zoST=0.056*d50*0.68*tstar/(1.0+coef_st*tstar);
   if (zoST<0),
      fprintf(1,' Warning: zoST<0  tstar, d50, coef_st:\n');
      fprintf(1,'%d %d %f %f\n',i,j,tstar,d50,coef_st);
   end
   %
   %  Calculate ripple height and wavelength.
   %  Use Malarkey & Davies (2003) explict version of Wiberg & Harris.
   %
   coef_b1=1.0/coef_a1;
   coef_b2=0.5*(1.0 + coef_a2)*coef_b1;
   coef_b3=coef_b2^2-coef_a3*coef_b1;
   d0=2.0*Ab;
   if ((d0/d50)>13000.0),             % sheet flow
      rheight=0.0;
      rlength=535.0*d50   ;      % does not matter since rheight=0
   else
      dolam1=d0/(535.0*d50);
      doeta1=exp(coef_b2-sqrt(coef_b3-coef_b1*log(dolam1)));
      lamorb=0.62*d0;
      lamanorb=535.0*d50;
      if (doeta1<20.0),
         dolam=1.0/0.62;
      elseif(doeta1>100.0),
         dolam=dolam1;
      else
         fdo_etaano=-log(lamorb/lamanorb)*log(0.01*doeta1)/log(5.0);
         dolam=dolam1*exp(-fdo_etaano);
      end
      doeta2=exp(coef_b2-sqrt(coef_b3-coef_b1*log(dolam)));
      rheight=d0/doeta2;
      rlength=d0/dolam;
   end
end


