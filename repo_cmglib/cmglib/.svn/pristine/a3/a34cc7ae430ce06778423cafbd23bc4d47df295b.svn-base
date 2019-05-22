function dn = j2dn( time,time2 )
% j2dn - Convert Julian time, time2 to Matlab datenum
[nr nc]=size(time);
time = time(:);
time2 = time2(:);
dn = datenum(gregorian(double(time)))+double(time2)./(24*3600*1000);
dn=reshape(dn,nr,nc);
end

