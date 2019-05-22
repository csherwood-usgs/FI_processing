H, T, h, Bf, gam
g = 9.81;
Lo = g*T*T/(2*pi)
Co = g*T/(2.*pi);
  Cgo = 0.5*Co; 
   if(break_depth)
	    % evaluate at break depth
	    Hb = H;
			hb = Hb/gam;
			kh = qkhfs( 2.*pi/T, hb );
			k = kh/hb;
			Ks = 1./ sqrt( tanh(kh)*(1.+2*kh/sinh(2*kh)) );
			Ho = Hb/Ks;
   end
	 if(deep_depth)
	    % evaluate in deep water
	    Ho = H;
			%Hb = 0.39*Math.pow(g,(1./5.))*Math.pow( T*Ho*Ho, (2./5.) ); % Komar eqn 6.6
			Hb = (gam/g).^(1./5.) * (Ho*Ho*Cgo).^(2./5.); % Plant et al (1999) eqn 8
      hb = Hb/gam;
			kh = qkhfs( 2.*pi/T, hb );
			k = kh/hb;
		  Ks = 1./ sqrt( tanh(kh)*(1.+2*kh/sinh(2*kh)) );
    end
	 if (intermed_depth)
	    % evaluate in intermediate depths
	    kh = qkhfs( 2.*pi/T, h );
		  k = kh/h;
			Ks = 1./ sqrt( tanh(kh)*(1.+2*kh/sinh(2*kh)) );
			Ho = H/Ks;
		  %Hb = 0.39*Math.pow(g,(1./5.))*Math.pow( T*Ho*Ho, (2./5.) ); % Komar eqn 6.6
			Hb = (gam/g).^(1./5.) * (Ho*Ho*Cgo).^(2./5.) ; % Plant et al (1999) eqn 8
			hb = Hb/gam;
    end
  I = Bf*(sqrt(Lo/Ho));
  eta = 0.35*Bf*sqrt(Ho*Lo);   % Eqn. 10
  Sinc = 0.75*Bf*sqrt(Ho*Lo);  % Eqn. 11
  SIG =  0.06*sqrt(Ho*Lo);     % Eqn. 12
  R2 = 1.1*(eta+0.5*sqrt(Sinc.^2. + SIG.^2)); % Eqns 6 and 7
   if (I<0.03)
      R2 = 0.043*sqrt(Ho*Lo);
   end
  x = R2/(Math.sin(Math.atan(Bf)));