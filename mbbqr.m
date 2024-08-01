function[se qs bvar]=mbbqr(y,x,bsize,nboot,tau,inq)

% y: dependent variable
% x: independent variables including constant
% bsize: Block size for resampling
% nboot: Number of resamples
% tau: Quantile

% programmed by Maik Wolters (Estimating monetary policy reaction
% functions using quantile regression)


warning off;

% Moving block bootstrapping
rand('state',2222); % seed for random number generation
nobs = size(x,1); % number of observations
nx = size(x,2); % number of regressors
bnum = nobs - bsize + 1; % number of possible blocks
m = int8(nobs/bsize); % numbers of blocks to be drawn 
bmean = zeros(nx,1); % initialize vector to save estimates from resampling
bprod = zeros(nx,nx);
betas=zeros(nboot,nx);
for iboot = 1 : nboot % nboot: number of resamples to be drawn
     l = 1;
     for i = 1 : m % draw m blocks for the resample
          f = ceil(bnum.*rand(1,1)); % choose block number 
          xres(l:l+bsize-1,:) = x(f:f+bsize-1,:); % save resampled data
          yres(l:l+bsize-1,:) = y(f:f+bsize-1,:);
          l = l + bsize;
     end;
     savest = rq(xres(1:l-1,:),yres(1:l-1,:),tau); 
     bmean = savest + bmean;
     betas(iboot,:) = savest';
     bprod = bprod + savest*savest';
end
bmean = bmean/nboot;
bvar = bprod/nboot-bmean*bmean'; % centering around sample estimate
se = sqrt(diag(bvar));
qs=quantile(betas,inq)';