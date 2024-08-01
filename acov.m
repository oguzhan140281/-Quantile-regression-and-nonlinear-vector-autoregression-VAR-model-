function acov=acov(x,y,l)
y=y(:); x=x(:); T=size(x,1);
acov=zeros(2*l+1,1);
c=1;
for h=-l:l
    effT=T-abs(h);
    bx=(abs(h)+1)*(h<0)+(h>=0); by=(abs(h)+1)*(h>0)+(h<=0);
    ex=T*(h<0)+effT*(h>=0); ey=T*(h>0)+effT*(h<=0);
    acov(c)=x(bx:ex,1)'*y(by:ey,1)/T;
    c=c+1;
end