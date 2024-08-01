function mic1=mic1(y,p,kmax,plot_opt)
nt=size(y,1);
z=ones(nt,1); cbar=-7.0; % for p=0
if p==1
    tr_p=[1:nt]'; z=[z tr_p];
    cbar=-13.5;
end
[yt,ssra]=glsd(y,z,cbar); % transforming the data
dyt=diff(yt);
for k=1:kmax
    DY=[];
    effT=length(y)-k-1;
    e=nt-1; b=k+1;
    for j=0:k
        DY=[DY dyt(b-j:e-j,:)];
    end
    Y=DY(:,1); X=[y(k+1:nt-1,1) DY(:,2:end)];
    beta=inv(X'*X)*X'*Y;
    beta0=beta(1,1);
    err=Y-X*beta;
    s2=sum(err(kmax+1:end).^2)/(nt-kmax);
    tau_T=sum(yt(kmax:nt-1).^2)*beta0^2/s2;
    MAIC(k,1)=log(s2)+2*(tau_T+k)/(nt-kmax);
    MBIC(k,1)=log(s2)+log(nt-kmax)*(tau_T+k)/(nt-kmax);
end
[v,IB]=min(MBIC); [v,IA]=min(MAIC);
if plot_opt==1
    plot(MBIC,'*r')
    hold on
    plot(MAIC,'*b')
    hold off
    legend('MBIC','MAIC')
    grid on
end
mic1=[IA IB];