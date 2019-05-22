function rss=reynolds_stress( au,av,aw,bu,bv,bw,fs, iverbose )
% REYNOLDS_STRESS - Calculate <u'w'> and <v'w'> using Trowbridge method
% rss=reynolds_stress( au,av,aw,bu,bv,bw,fs, iverbose  )

if(exist('iverbose')~=1),iverbose=0,end;
DirA = 0; not sure what this means
ialign = 0;


[sa da]=pcoord(mean(au),mean(av))
[sb db]=pcoord(mean(bu),mean(bv))
[sab dab]=pcoord(mean([au;bu]),mean([av;bv]))
% rotate so u = mean directio of A
if(ialign)
   da_r = da-DirA+90;
   % rotate to +u = mean flow direction    db_r = db-DirA+90;
   [u_ar,v_ar]=xycoord(sa,da_r);
   [u_br,v_br]=xycoord(sb,db_r);
   xd = u_ar(:)-u_br(:);
   yd = v_ar(:)-v_br(:);
   zd = aw(:)-bw(:);
else
   xd = 
xd = xd(:);
yd = yd(:);
zd = zd(:);
n = length(xd);
n4 = floor(n/4);

cxw = cov([xd zd],1);
rss.rsu = -0.5*cxw(1,2);
cyw = cov([yd zd],1);
rss.rsv = -0.5*cyw(1,2);
end


if(iverbose)
   fprintf(1,'sa, da: %g %g\n',sa,da)

   if(iverbose) % do extra stuff
      n = length(xd);
      n4 = floor(n/4)
      [xxz,lags]=xcorr(xd,zd,n4,'coeff');
      [xyz,lags]=xcorr(yd,zd,n4,'coeff');
      [ax,lags]=xcorr(xd,n4,'coeff');
      [ay,lags]=xcorr(yd,n4,'coeff');
      [az,lags]=xcorr(zd,n4,'coeff');
      lag0 = find(lags==0);
      figure(1)
      subplot(211); clf
      plot(lags/fs,xxz,'-m')
      hold on
      plot(lags/fs,xyz,'-c')
      plot(lags(1:lag0)/fs,ax(1:lag0),'-r')
      plot(lags(lag0:end)/fs,ay(lag0:end),'-b')
      drawnow
      shg

      Tx = cumsum(ax(lag0:end,1))./ax(lag0);
      Ty = cumsum(ay(lag0:end,1))./ay(lag0);
      Tw = cumsum(az(lag0:end,1))./az(lag0);


      % COSPECPLOT -Plot cumulative cospectra
      nseg = 16;
      nfft = length(xd)/nseg;
      [Pxz,ff]=csd(xd,zd,nfft,fss,hanning(nfft),'mean');
      [Pyz,ff]=csd(yd,zd,nfft,fss,hanning(nfft),'mean');
      dff = ff(3)-ff(1);
      %figure(1)
      %clf
      subplot(211)
      semilogx(ff,real(Pxz),'-r','linewidth',2)
      hold on
      semilogx(ff,real(Pyz),'-b','linewidth',2)
      plot([.02;10],[0 0],'--k')
      axis([.02 10 -.05 .05])
      set(gca,'xtick',[.02 .1 1 10]);
      axis('square')
      % unsmoothed estimates for cumulative plots
      nseg = 1;
      nfft = length(xd)/nseg;
      [Pxz,ff]=csd(xd,zd,nfft,fss,boxcar(nfft),'mean');
      [Pyz,ff]=csd(yd,zd,nfft,fss,boxcar(nfft),'mean');
      dff = ff(3)-ff(1);
      ylabel('Co-spectra (m^2/s^2/Hz)')
      subplot(212)
      semilogx(ff,cumsum(real(Pxz)*dff),'-r','linewidth',2)
      hold on
      semilogx(ff,cumsum(real(Pyz)*dff),'-b','linewidth',2)
      plot([.02;10],[0 0],'--k')
      axis([.02 10 -.02 .02])
      axis('square')
      set(gca,'xtick',[.02 .1 1 10])
      xlabel('Frequency (Hz)')
      ylabel('Cumulative Co-spectra (m^2/s^2)')


   end