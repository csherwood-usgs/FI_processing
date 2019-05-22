function B=bailard(u,v,fs,fc,slopeaz,tanBeta,ws)
% BAILARD - Energetics sediment transport calculation
% Approach of Bailard, as modified by Wright et al.
% 
% Input: u and v components of velocity (by definition, u positive toward shore)
%        tanBeta = tangent of nearshore slope (rise/run)
% MKS units

% Wright et al., 1991, Marine Geology 96:19-51, eqn. 10.
% See also Bailard(1981) JGR 86(C11):10,938-10,954, eqn. 11.

[s,d]=pcoord(u,v);
ubar = mean(u);
vbar = mean(v);
uf = lpfilt( (u-ubar),1/fs,fc );
vf = lpfilt( (v-vbar),1/fs,fc );
uw = (u-uf)-ubar;
vw = (v-vf)-vbar;

eb = 0.21;  % Bailard, p. 10942
es = 0.025;
tanPhi = tan( 10*pi/180 );
Cf =.001;
rho = 1035;

if(exist('tanBeta')~=1),
  tanBeta = .005;
end
if(exist('ws')~=1),
  ws = .012; %2.75 phi
end

blcoef = rho*Cf*eb/tanPhi;
sscoef = rho*Cf*es/ws;

U2 = s.^2;;
U3 = U2.*s;
B.Ibl_ubar = blcoef*mean(U2)*ubar;
B.Ibl_uw   = blcoef*mean(U2 .* uw);
B.Ibl_uf   = blcoef*mean(U2 .* uf);
B.Ibl_g    =-blcoef*mean(U3)*tanBeta/tanPhi;
B.Ibl_vbar = blcoef*mean(U2)*vbar;
B.Ibl_vw   = blcoef*mean(U2 .* vw);
B.Ibl_vf   = blcoef*mean(U2 .* vf);
B.Iss_ubar = sscoef*mean(U3)*ubar;
B.Iss_uw   = sscoef*mean(U3.*uw);
B.Iss_uf   = sscoef*mean(U3.*uf);
B.Iss_g    =-sscoef*mean(s.^5)*es/ws;
B.Iss_vbar = sscoef*mean(U3)*vbar;
B.Iss_vw   = sscoef*mean(U3.*vw);
B.Iss_vf   = sscoef*mean(U3.*vf);
B.Ibl_totu = B.Ibl_ubar+B.Ibl_uw+B.Ibl_uf+B.Ibl_g;
B.Ibl_totv = B.Ibl_vbar+B.Ibl_vw+B.Ibl_vf;
B.Iss_totu = B.Iss_ubar+B.Iss_uw+B.Iss_uf+B.Iss_g;
B.Iss_totv = B.Iss_vbar+B.Iss_vw+B.Iss_vf;


