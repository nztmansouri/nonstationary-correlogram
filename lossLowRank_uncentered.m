
function [f,dx] = lossLowRank_uncentered(b,y,k,offset,nu)

if nargin<5, offset=0; end
if nargin<6, nu=0; end

[n,m]=size(y);

A = reshape(b(1:(n*k)),n,k);
B = reshape(b((n*k+1):end),k,m);

xb = A*B + offset;
lam = exp(xb);

f = -nansum(nansum(y.*xb - lam)) + nu*sum(sum(abs(A))) + nu*sum(sum(abs(B)));
lam_err = lam-y;
lam_err(~isfinite(lam_err))=0;
dA = (lam_err)*B'+nu*sign(A);
dB = A'*(lam_err)+nu*sign(B);
dx = [dA(:); dB(:)];

