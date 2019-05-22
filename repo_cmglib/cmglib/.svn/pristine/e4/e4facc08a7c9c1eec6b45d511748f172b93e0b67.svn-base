function reptxt( fnin, txtin, fnout, txtout )
% REPTXT - Replace text txtin in file fnin with text txtout in file fnout
% reptxt( fnin, txtin, fnout, txtout )

% csherwood@usgs.gov
% April, 2005

fprintf('Looking for %s in %s\n',txtin,fnin)
fprintf('Will replace it with %s in %s\n',txtout,fnout)

fidin = fopen( fnin, 'r' );
c = fread( fidin, 'char=>char' )';
fclose(fidin);

ni = strfind( c, txtin );
if( ~any(ni) ),
    fprintf('No occurrences of %s found in %s\n',txtin,fnin);
    rc = c;
elseif (length(ni)>1),
    fprintf('Found more than one (actually %d) occurrences of %s in %s\n',...
        length(ni),txtin,fnin);
    rc = c;
else
    nchar = length(txtin);
    rc = [c(1:ni-1) txtout c(ni+nchar:end)];
end
fidout = fopen( fnout,'w');
fprintf(fidout,'%s',rc);
fclose(fidout);
%fprintf(1,'Wrote %s\n',fnout);