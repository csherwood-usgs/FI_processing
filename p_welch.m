function [Pxx, freq] = p_welch(x,nfft,Fs,win,noverlap)
% P_WELCH Power Spectral Density estimate via Welch's method.
% [Pxx, freq] = P_WELCH(x,nfft,Fs,win,noverlap)
%
% Version from older Matlab Toolbox, revised for use in Meg
% Palmsten's routine.
%
% The single-sided PSD is returned for real signals.
% No error checking, and limited options.
%
% [Pxx,F] = P_WELCH(X,NFFT,Fs,WIN,NOVERLAP) X is divided into overlap-
% ping sections, then windowed by the WIN parameter, then zero-padded
% to length NFFT.  The magnitude squared of the length NFFT DFTs of the 
% sections are averaged to form Pxx.  Pxx is length NFFT/2+1 for NFFT
% even, (NFFT+1)/2 for NFFT odd, or NFFT if the signal X is complex. If  
% you specify a scalar for WIN, a Hanning window of that length is 
% used.

n = length(x);           % Number of data points
nwind = length(win);     % length of window
if(nwind==1),
  window = hanning(win);
  nwind = win;
else
  window = win;
end
if n < nwind            
   x(nwind)=0;
   n=nwind;
end
% Make x a column vector; do this AFTER the zero-padding in case x is a scalar.
x = x(:);

% Number of windows; (k = fix(n/nwind) for noverlap=0)
k = fix((n-noverlap)/(nwind-noverlap)); 

index = 1:nwind;
KMU = k*norm(window)^2; % Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

% Calculate PSD
Spec = zeros(nfft,1);
for i=1:k
   xw = window.*x(index);
   index = index + (nwind - noverlap);
   Xx = abs(fft(xw,nfft)).^2;
   Spec = Spec + Xx;
end

% Select first half
if rem(nfft,2),         % nfft odd
  select = (2:(nfft+1)/2-1)';  % don't include DC or Nyquist components
  nyq    = (nfft+1)/2;         % they're included below
else
  select = (2:nfft/2)';
  nyq    = nfft/2+1; 
end
% Calculate the single-sided spectrum which includes the full power
Spec = [Spec(1); 2*Spec(select); Spec(nyq)];
freq = [0; Fs*(1:nfft/2)'/nfft ]; 

% Normalize, and scale by 1/Fs to get 
% Power per unit of frequency 
Pxx = Spec / (KMU * Fs);
