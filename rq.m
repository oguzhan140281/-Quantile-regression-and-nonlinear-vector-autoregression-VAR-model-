function [coef rho resid serr quants cov tvalue pvalue] = rq(X, y, tau, se, covariance, hs, varargin)

% Construct the dual problem of quantile regression
% Solve it with lp_fnm
% 
% Function rq_fnm of Daniel Morillo & Roger Koenker
% Found at: http://www.econ.uiuc.edu/~roger/rqn/rq.ox
% Translated from Ox to Matlab by Paul Eilers 1999
%
% Extended by Maik Wolters based on R code, 2011
%
if nargin < 4;
    se='false';
end
if nargin < 6;
    hs=1;
end


[m n] = size(X);
u = ones(m, 1);
a = (1 - tau) .* u;

%parameters:
coef = -lp_fnm(X', -y', X' * a, u, a)';

%residuals:
resid=y-X*coef;

%value of objective function:
rho = tau*sum(resid(resid>0))-(1-tau)*sum(resid(resid<0));

%standard errors
%	    "nid" standard errors which presumes local (in tau) linearity (in x) of the
%		the conditional quantile functions and computes a Huber
%		sandwich estimate using a local estimate of the sparsity.
%
% other options for standard errors can be translated from the R code in
% summary.rq


if strcmp(se,'nid')
    n = length(resid);
    p = length(coef);
    rdf = n - p;
    h = bandwidthrq(tau, m, hs);
    if tau + h > 1
        disp('tau + h > 1:  error in the computation of ''nid'' standard errors, set tau+h=0.99');
        bhi = rq(X, y, 0.99);  %coef
        blo = rq(X, y, tau - h);  %coef
    elseif tau - h < 0
        disp('tau - h < 0:  error in the computation of ''nid'' standard errors, set tau-h=0.01');
        bhi = rq(X, y, tau + h);  %coef
        blo = rq(X, y, 0.01);  %coef
    else
        bhi = rq(X, y, tau + h);  %coef
        blo = rq(X, y, tau - h);  %coef
    end
    dyhat = X * (bhi - blo);
    if min(dyhat) <= 0
        disp('sum(dyhat) <= 0), ''non-positive fis''');
    end
    temp1= zeros(size(dyhat,1),1);
    temp2=(2 * h)./(dyhat - eps);
    f = max([temp1 temp2],[],2);
    fxxinv = eye(p);
    for i = 1:size(X,2)
        temp1(:,i)=sqrt(f).*X(:,i);
    end
    temp2=qr(temp1);
    temp3=temp2(1:p,1:p);
    for jnr=2:size(temp3,1)
        for inr=1:jnr-1
            temp3(jnr,inr)=0;
        end
    end
    fxxinv = temp3\fxxinv;
    fxxinv = fxxinv * fxxinv';
    cov = tau * (1 - tau) * fxxinv * X'*X *fxxinv;
    scale = mean(f);
    serr = sqrt(diag(cov));

    %disp standard errors
    tvalue=coef./serr;
    if rdf>0
        pvalue = 2 * (1 - tcdf(abs(tvalue), rdf));
    else 
        pvalue = 'NaN';
    end
elseif strcmp(se,'boot')
    nboot = varargin{1};
    bsize = varargin{2};
    [serr quants]=mbbqr(y,X,bsize,nboot,tau,[0.025 .975]);
    cov=[]; tvalue=[]; pvalue=[];
end
    


function y = lp_fnm(A, c, coef, u, x)
% Solve a linear program by the interior point method:
% min(c * u), s.t. A * x = coef and 0 < x < u
% An initial feasible solution has to be provided as x
%
% Function lp_fnm of Daniel Morillo & Roger Koenker
% Found at: http://www.econ.uiuc.edu/~roger/rqn/rq.ox
% Translated from Ox to Matlab by Paul Eilers 1999
% Modified slightly by Roger Koenker, April, 2001.

% Set some constants
  beta = 0.9995;
  small = 1e-5;
  max_it = 50;
  [m n] = size(A);

% Generate inital feasible point
  s = u - x;
  y = (A' \  c')';
  r = c - y * A;
  z = r .* (r > 0);
  w = z - r;
  gap = c * x - y * coef + w * u;

% Start iterations
  it = 0;
  while gap > small & it < max_it
    it = it + 1;

%   Compute affine step
    q = 1 ./ (z' ./ x + w' ./ s);
    r = z - w;
    Q = sparse(1:n,1:n,q);
    AQ = A * Q;
    AQA = AQ * A';
    rhs = AQ * r';
    dy = (AQA \ rhs)';
    dx = q .* (dy * A - r)';
    ds = -dx;
    dz = -z .* (1 + dx ./ x)';
    dw = -w .* (1 + ds ./ s)';
    
%   Compute maximum allowable step lengths
    fx = bound(x, dx);
    fs = bound(s, ds);
    fw = bound(w, dw);
    fz = bound(z, dz);
    fp = min(fx, fs);
    fd = min(fw, fz);
    fp = min(min(beta * fp), 1);
    fd = min(min(beta * fd), 1);

%   If full step is feasible, take it. Otherwise modify it
    if min(fp, fd) < 1
    
%     Update mu
      mu = z * x + w * s;
      g = (z + fd * dz) * (x + fp * dx) + (w + fd * dw) * (s + fp * ds);
      mu = mu * (g / mu) ^3 / ( 2* n);

%     Compute modified step
      dxdz = dx .* dz';
      dsdw = ds .* dw';
      xinv = 1 ./ x;
      sinv = 1 ./ s;
      xi = mu * (xinv - sinv);
      rhs = rhs + A * ( q .* (dxdz - dsdw -xi));
      dy = (AQA \ rhs)';

      dx = q .* (A' * dy' + xi - r' -dxdz + dsdw);
      ds = -dx;
      dz = mu * xinv' - z - xinv' .* z .* dx' - dxdz';
      dw = mu * sinv' - w - sinv' .* w .* ds' - dsdw';

%     Compute maximum allowable step lengths
      fx = bound(x, dx);
      fs = bound(s, ds);
      fw = bound(w, dw);
      fz = bound(z, dz);
      fp = min(fx, fs);
      fd = min(fw, fz);
      fp = min(min(beta * fp), 1);
      fd = min(min(beta * fd), 1);
      
    end

%   Take the step
    x = x + fp * dx;
    s = s + fp * ds;
    y = y + fd * dy;
    w = w + fd * dw;
    z = z + fd * dz;

    gap = c * x - y * coef + w * u;
%     disp(gap);

end
    
      
    
  
function b = bound(x, dx)
% Fill vector with allowed step lengths
% Support function for lp_fnm
b = 1e20 + 0 * x;
f = find(dx < 0);
b(f) = -x(f) ./ dx(f);
