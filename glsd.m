% ERS detrending
function [yt,ssr]=glsd(y,z,cbar)
nt=size(y,1);
abar=1+cbar/nt;
ya=zeros(nt,1);
za=zeros(nt,size(z,2));
ya(1,1)=y(1,1); za(1,:)=z(1,:);
ya(2:nt,1)=y(2:nt,1)-abar*y(1:nt-1,1);
za(2:nt,:)=z(2:nt,:)-abar*z(1:nt-1,:);
bhat=inv(za'*za)*za'*ya;
yt=y-z*bhat;
ssr=(ya-za*bhat)'*(ya-za*bhat);