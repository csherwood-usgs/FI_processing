% ALIGNT
% function idx = align( varargin )

c = { jtime, jta };

[nr,nc]=size(c);
lo = zeros(nc,2);
hi = zeros(nc,2);
% find first and last values and their indexes
for i=1:nc,
    [lo(i,1) lo(i,2)] = min( c{i} );
    if(lo(i,2)~=1)
        fprintf(1,'Lowest value of array %d is not first.\n',i)
    end
    [hi(i,1) hi(i,2)] = max( c{i} );
    if(hi(i,2)~=length(c{i})),
        fprintf(1,'Highest value of array %d is not last.\n',i)
    end
end

% last start time
[fs,fi]=max(lo(:,1))
% first end time
[le,li]=min(hi(:,1))
