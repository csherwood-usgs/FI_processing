% SPEEDPROFILE - Fit log profiles to PCADP data
% Version includes multiple elevations
% Last revised on old desktop in CMGLIB 9/24/2005
% dnm = 'c:\crs\proj\GH2001\MIA\pcadp\'
% nc1 = netcdf([dnm,'mia1pcb-cal.nc'])
% nc2 = netcdf([dnm,'mia1pcvp-cal.nc'])

dnm = 'd:\crs\data\EuroSTRAT\702-ch10-flow\7022-pc30\'
nc1 = netcdf([dnm,'pc7022b-cal.nc']);
nc2 = netcdf([dnm,'pc7022vp-cal-bfix.nc']);

% dnm = 'g:\data\EuroSTRAT\712-ch10-flow\7122-pc30\'
% nc1 = netcdf([dnm,'pc7122b-cal.nc']);
% nc2 = netcdf([dnm,'pc7122vp-cal-bfix.nc']);

% dnm = 'g:\data\EuroSTRAT\703-ch20-gp\7033-pc40\'
% nc1 = netcdf([dnm,'pc7033b-cal.nc']);
% nc2 = netcdf([dnm,'pc7033vp-cal-bfix.nc']);

% dnm = 'g:\data\EuroSTRAT\713-ch20-gp\7133-pc40\'
% nc1 = netcdf([dnm,'pc7133b-cal.nc']);
% nc2 = netcdf([dnm,'pc7133vp-cal-bfix.nc']);


months = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';
          'Jul';'Aug';'Sep';'Oct';'Nov';'Dec'];

profiles = nc1.ProfilesPerBurst(:);
cells = nc1.PCADPUserSetupNcells(:);
cellsize = (nc1.PCADPUserSetupCellSize(:));
blank = nc1{'depth'}.blanking_distance(:);

halfway = fix(profiles/2);
jday = nc1{'time'}(:, halfway); %just want 1 time column per burst, using middle profile
msec = nc1{'time2'}(:, halfway);
jtime = jday + msec/(1000*60*60*24); %julian decimal days
nbursts = length(jtime);
pitch = nc2{'Ptch_1216'}(:,1,1);
roll = nc2{'Roll_1217'}(:,1,1);
hdg = nc2{'Hdg_1215'}(:,1,1);

if(1),
    brange = nc2{'brange'}(:,:,1,1);
    brange(:,1)= fillnan(1:nbursts,brange(:,1));
    brange(:,2)= fillnan(1:nbursts,brange(:,2));
    brange(:,3)= fillnan(1:nbursts,brange(:,3));
    brange = vert_bd(brange, pitch, roll);
    medr = median(brange')';
    minr = gmin(brange')';
    % next two lines fix dropouts in brange
    uout = thumbfin(minr,.8,1.2,1);
    minr = fillnan(jtime,uout);
    maxr = gmax(brange')';
end
%%
bno = (1:nbursts)';
ut = NaN*ones(nbursts, cells);
vt = NaN*ones(nbursts, cells);
wt = NaN*ones(nbursts, cells);
binzt = NaN*ones(nbursts, cells);

us4 = NaN*ones(nbursts,1);
zo4 = us4;
r24 = us4;
r2a4 = us4;
use4 = us4;
zoe4 = us4;
P4 = NaN*ones(nbursts,2);
P24 = NaN*ones(nbursts,3);

us4_min = NaN*ones(nbursts,1);
zo4_min = us4;
r24_min = us4;
r2a4_min = us4;
use4_min = us4;
zoe4_min = us4;
P4_min = NaN*ones(nbursts,2);
P24_min = NaN*ones(nbursts,3);

us4_max = NaN*ones(nbursts,1);
zo4_max = us4;
r24_max = us4;
r2a4_max = us4;
use4_max = us4;
zoe4_max = us4;
P4_max = NaN*ones(nbursts,2);
P24_max = NaN*ones(nbursts,3);

n9 = NaN*ones(nbursts,1);
us9 = NaN*ones(nbursts,1);
zo9 = us4;
r29 = us4;
r2a9 = us4;
use9 = us4;
zoe9 = us4;
P9 = NaN*ones(nbursts,2);
P29 = NaN*ones(nbursts,3);

bh_min_off = medr-minr;
bh_max_off = medr-maxr;

%for ib = (gb(1):gb(end))
for ib = (1:nbursts),
    fprintf(1,'%d\n',ib)
    binheights = medr(ib)-(blank+cellsize*([0:(cells-1)]')+0.5*cellsize);    
    binzt(ib,:) = binheights';
    u = nc1{'u_1205'}(ib,:,:);  % dimensions 840x11
    v = nc1{'v_1206'}(ib,:,:);
    w = nc1{'w_1204'}(ib,:,:);

    u = 0.01*mean(u); %burst mean, should be dimesnsions 1x11, now m/s
    v = 0.01*mean(v);
    w = 0.01*mean(w);
    sp = sqrt([u.^2 + v.^2]); %speed

    % problems with 1.5 * bin size in part of this record...going to 2
    %lastcell = find( (binheights-(1.5*cellsize))>0 ,1,'last');
    lastcell = find( (binheights-(2*cellsize))>0 ,1,'last');
    ic4 = (lastcell-3:lastcell);
    if(length(ic4)==4),
        [us4(ib),zo4(ib),r24(ib),r2a4(ib),use4(ib),zoe4(ib),res4,P4(ib,:),P24(ib,:)]=...
            mlogfit(sp(ic4)', binheights(ic4));
        if(min(binheights(ic4)-bh_min_off(ib))>0.5*cellsize),
            [us4_min(ib),zo4_min(ib),r24_min(ib),r2a4_min(ib),use4_min(ib),...
                zoe4_min(ib),res4_min,P4_min(ib,:),P24_min(ib,:)]=...
                mlogfit(sp(ic4)', (binheights(ic4)-bh_min_off(ib)));
        end
        [us4_max(ib),zo4_max(ib),r24_max(ib),r2a4_max(ib),use4_max(ib),...
            zoe4_max(ib),res4_max,P4_max(ib,:),P24_max(ib,:)]=...
            mlogfit(sp(ic4)', (binheights(ic4)-bh_max_off(ib)));
    end
    ic9 = (1:lastcell);
    if(length(ic9)>5),
        n9(ib) = length(ic9);
        [us9(ib),zo9(ib),r29(ib),r2a9(ib),use9(ib),zoe9(ib),res9,P9(ib,:),P29(ib,:)]=...
            mlogfit(sp(ic9)', binheights(ic9));
    end

    if(1)
        figure(1)
        clf
        curve_fit = polyval(P29(ib,:),log(binheights(ic9)) );
        semilogy( curve_fit, binheights(ic9), '-c')
        hold on
        u_fit = log(binheights/zo4(ib))*us4(ib)/.408;
        semilogy(u_fit, (binheights), ':r')
        u_fit4 = log(binheights(ic4)/zo4(ib))*us4(ib)/.408;
        semilogy(u_fit4, (binheights(ic4)),'-r','linewidth',2)
        semilogy(sp',binheights, 'bo')
        ylabel('Elevation (m)')
        xlabel('Speed (m/s)')
        axis([0 1.4 .02 1.2])

        gtime=gregorian(jtime(ib));
        ts0 = sprintf('Burst % 4d on %4d-%s-%02d at %02d:%02d'...
            ,ib,gtime(1),months(gtime(2),:),gtime(3:5));
        zoval = 100*zo4(ib);
        r2val = r24(ib);
        r2aval= r2a4(ib);
        if( (us4(ib)<=0) )
            zoval = 9.999;
            r2val = 0.;
            r2aval = 0.;
        end
        ts1=sprintf('n: %2d: u*: %3.1f zo: %6.3f r2: %6.4f r2a: %6.4f',...
            4,100*us4(ib),zoval,r2val,r2aval);
        zoval = 100*zo9(ib);
        r2val = r29(ib);
        r2aval= r2a9(ib);
        if( (us9(ib)<=0) )
            zoval = 9.999;
            r2val = 0.;
            r2aval = 0.;
        end
        ts2=sprintf('n: % 2d: u*: %3.1f zo: %6.3f r2: %6.4f r2a: %6.4f',...
            n9(ib),100*us9(ib),zoval,r2val,r2aval);
        ts3=sprintf('p: %7.3f',P29(ib,1))
        xtxloc = .4;
        text(xtxloc,.1,ts0);
        text(xtxloc,.07,ts1);
        text(xtxloc,.055,ts2);
        h=text(xtxloc,.045,ts3);
        if(P29(ib,1)>0.01),
            set(h,'color',[0 0 1])
        elseif(P29(ib,1)<-0.01),
            set(h,'color',[1 0 0])
        end
        drawnow
        pn = sprintf('p%04d.png',ib)
        %eval(['print -dpng .\profs7133\',pn]);
    end

    %save variable to time series
    ut(ib, :) = u;
    vt(ib, :) = v;
    wt(ib, :) = w;
end
ncclose

save profiles7033 jtime ut vt wt binzt ...
    medr minr maxr ...
    us4 zo4 use4 zoe4 r24 r2a4 P4 P24 ...
    us4_min zo4_min use4_min zoe4_min r24_min r2a4_min P4_min P24_min ...
    us4_max zo4_max use4_max zoe4_max r24_max r2a4_max P4_max P24_max ...
    n9 us9 zo9 use9 zoe9 r29 r2a9 P9 P29

if(0)
    figure(2)
    clf
    set(0,'defaultaxesfontsize',14)
    set(0,'defaulttextfontsize',14)

    set(gcf, 'PaperPositionMode', 'manual');
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [1.25 1.75 6 8]);

    subplot(4,1,1)
    plot(jtime, ut(:,6), 'color', 'b', 'linewidth', 2)
    hold on
    plot(jtime, vt(:,6), 'color', 'g', 'linewidth', 2)
    ylimits = get(gca, 'ylim');
    timeaxDep1(ylimits(1), ylimits(2), ['a) PCADP H50 (mooring 7583): Velocity'], 0)
    ylabel('Velocity (m/s)')
    legend('u', 'v')

    subplot(4,1,2)
    plot(jtime, (us_use(:,1) + us_use(:,2)), 'color', [.7 .7 1], 'linewidth', 2)
    hold on
    plot(jtime, (us_use(:,1) - us_use(:,2)), 'color', [.7 .7 1], 'linewidth', 2)
    plot(jtime, us_use(:,1), 'b')
    ylimits = get(gca, 'ylim');
    timeaxDep1(ylimits(1), ylimits(2), ['b) u_* with error (m/s)'], 0)

    subplot(4,1,3)
    plot(jtime, (zo_zoe(:,1) + zo_zoe(:,2)), 'color', [.7 .7 1], 'linewidth', 2)
    hold on
    plot(jtime, (zo_zoe(:,1) - zo_zoe(:,2)), 'color', [.7 .7 1], 'linewidth', 2)
    plot(jtime, zo_zoe(:,1), 'b')
    ylimits = get(gca, 'ylim');
    timeaxDep1(ylimits(1), ylimits(2), ['c) z_o with error (m)'], 0)

    subplot(4,1,4)
    plot(jtime, r2t, 'b')
    ylimits = get(gca, 'ylim');
    timeaxDep1(ylimits(1), ylimits(2), ['d) r^2'], 1)
end
