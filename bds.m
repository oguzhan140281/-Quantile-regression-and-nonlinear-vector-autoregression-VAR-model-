yy = xlsread('data','Monthly Data','B2:D207');
%% Choose the series %%

x = yy(:,1); %



function c = bds(x,eps,m,format)
% usage: c = bds(x,eps,m,format)
%
% Detailed usage:
% bds(x,eps,m,format)
%
% x   = vector time series (n x 1)
% eps = epsilon in units of std of x
% m   = maxdim:  The maximum dimension, m, to estimate the statistic for.
%       BDS estimates will be done over the range 2-m
%
% format = 1 for normalized bds
% returns standardized stat distributed asymptotically N(0,1)
% c(2)/s(2), c(3)/s(3), c(4)s(4), .. c(m)/s(5)   where s(m) is the asymtotic std
%
% format = 2 for raw c(1), c(2), c(3)..., k
% raw c numbers with no normalization
%       (This is useful for bootstrapping distributions.)
%
%  Further information on BDS can be found in:

%
eps = eps*std(x);  % bds routines want raw epsilon
if nargin==4
   c = bdsc(x,eps,m,format);
else
   'Bad argument list'
end