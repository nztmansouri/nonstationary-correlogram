
function [A,B,R2] = poisson_decomp(y,rankk)

[n,p]=size(y);
options=[];
options.useMex=0;

A = randn(n,rankk)/10;
B = randn(rankk,p)/10;
A(:,1)=log(mean(y,2)+1)-mean(log(mean(y,2)+1));
B(1,:)=log(mean(y,1)+1);

% optimize shape and coefficients simultaneously
[x,f,exitflag,output] = minFunc(@lossLowRank_uncentered,[A(:); B(:)],options,y,rankk);
A = reshape(x((1):(n*rankk)),n,rankk);
B = reshape(x((n*rankk+1):end),rankk,p);
yhat = exp(A*B);

% deviances
Ls = sum(y(:).*log(y(:)+(y(:)==0))-y(:));
Dm = 2*(Ls+f);
D0 = 2*(Ls-sum(y(:).*log(mean(y(:)))-mean(y(:))));
R2 = (D0-Dm)/D0;