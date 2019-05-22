% TEST_HLPRED - Script to test HL_PRED versions

D=logspace(1,3)*1e-6;  % 10-1000 micron
D=D(:);
d0=0.1*ones(size(D));
x=d0./D;
nD=length(x);
for i=1:nD
  r1(i) = hlpred(D(i),d0(i));
  r2(i) = hlpred2(D(i),d0(i));
end
h1 = reshape([r1.H],nD,1);
h2 = reshape([r2.eta],nD,1);
lam1 = reshape([r1.lam],nD,1);
lam2 = reshape([r2.lamda],nD,1);

clf
figure(1)
subplot(121)
h=loglog(x,[lam1(:)./D(:) lam2(:)./D(:)],'-o');
set(h(2),'marker','x');
legend('Implicit','Explicit');
xlabel('d_o/D');
ylabel('\lambda/D')
title('Comparison of Explicit and Implicit methods for Wiberg-Harris Ripples');
grid

subplot(122)
h=loglog(x,[h1(:)./lam1(:) h2(:)./lam2(:)],'-o');
set(h(2),'marker','x');
legend('Implicit','Explicit');
xlabel('d_o/D');
ylabel('\eta/\lambda')
grid
