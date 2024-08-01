clc; clear all; close all;
%I modified the code of Nusair, S. A. & Olson, D. (2019). The effects of oil price shocks on Asian exchange rates: Evidence from quantile regression analysis. Energy Economics, 78, 44–63. 
Model=2; % 1 for constant mean; 2 for trend
yy = xlsread('data2','Monthly Data','B2:D207');


series = yy(:,1); %russia
%%series = yy(:,2); %korea
%%series = yy(:,3); %fsi


logseries_o = series;

nt = size(logseries_o,1);
z = ones(nt,1); tr_p=[1:nt]';   if Model==2; z=[z tr_p]; end;
dt = inv(z'*z)*z'*logseries_o;
St = logseries_o - z*dt;
logseries = logseries_o;
kmax = 1; plot_opt = 0;
q = mic1(logseries,Model,kmax,plot_opt);
nlag = max(q);
nobs = size(series,1);
vectau = 0.05:0.05:0.95;
Ntau = size(vectau,2);
bsize = 4;
nboot = 1000;

% Construct data matrices as in the augmented Dickey Fuller form
    y_o = logseries_o;               % dependent variable
    y = logseries(nlag+1:nobs);      % dependent variable
    x = zeros(nobs - nlag, nlag); % independent variable
    x(:,1) = logseries(nlag:nobs-1); % AR(1) term

    for i = 2:nlag
        x(:,i) = logseries(nlag+2-i:nobs+1-i)-logseries(nlag+1-i:nobs-i); %lags of delta(logseries)
    end

S = St(nlag+1:nobs);           % dependent variable
st = zeros(nobs - nlag, nlag); % independent variable
st(:,1) = St(nlag:nobs-1);     % AR(1) term
    for i = 2:nlag
        st(:,i) = St(nlag+2-i:nobs+1-i)-St(nlag+1-i:nobs-i); %lags of delta(logseries)
    end
Ts = size(S,1);
beta_s = inv(st'*st)*st'*S;
e_res = S-st*beta_s;
dom = linspace(-0.04,.04);

% Run Quantile Regression

Decision=[];
for i=1:Ntau
    vtau  = vectau(i);
    tr_p = [1:max(size(x))]';
    if Model==1; tr_p = []; end;
    coeffs = rq([ones(size(x,1),1) x tr_p], y, vtau);
    res = y-[ones(size(x,1),1) x tr_p]*coeffs;
    psi_tau = vtau-(res<0);
    delta0 = corr(e_res,psi_tau);
    M = floor(nt^(1/3));
    
    ac_e_res = acov(e_res,e_res,M); 
    ac_e_resf = acov(e_res,psi_tau,M); 
    ac_p = acov(psi_tau,psi_tau,M);
    h = -M:M;
    kh = 25*M^2*(sin((6*pi*h)/(5*M))./((6*pi*h)/(5*M))-cos(6*pi*h/(5*M)))./(12*pi^2*h.^2);
    kh(M+1) = 1;
    delta_2 = (kh*ac_e_resf/(sqrt(vtau*(1-vtau))*sqrt(kh*ac_e_res)))^2;
    DD(i) = delta_2;
        if delta_2<.1
            delta_2 = .1;
        end
         if delta_2>1
            delta_2 = 1;
        end

            switch Model
                case 1
                    table1 = [1:-.1:.1; -2.86 -2.81 -2.75 -2.72 -2.64 -2.58 -2.51 -2.40 -2.28 -2.12]';
                case 2
                    table1 = [1:-.1:.1; -3.41 -3.33 -3.27 -3.18 -3.10 -2.99 -2.87 -2.73 -2.55 -2.31]';
            end
    CV05(i) = interp1(table1(:,1),table1(:,2),delta_2);
    % we need the t-stat now.
    check = 1;
    za = norminv(1-0.025); fFt=normpdf(norminv(vtau)); nn=size(x,1);
    h_n1(i) = nn^(-1/3)*za^(2/3)*(1.5*fFt^2/(4*norminv(vtau)^2+1))^(1/3);
    h_n2(i) = nn^(-1/5)*(4.5*fFt^4/(4*norminv(vtau)^2+1))^(1/5);
    h_n(i) = h_n2(i);
    coeffs_u = rq([ones(size(x,1),1) x tr_p], y, vtau+h_n(i));
    coeffs_l = rq([ones(size(x,1),1) x tr_p], y, vtau-h_n(i));
    Y_1 = x(:,1); xX = [ones(size(x,1),1) x(:,2:end) tr_p]; 
    PX = eye(size(x,1))-xX*inv(xX'*xX)*xX';
    fF1(i) = (2*h_n(i)/(mean([ones(size(x,1),1) x tr_p])*(coeffs_u-coeffs_l)));
    t_stat(i) = fF1(i)*sqrt(Y_1'*PX'*Y_1)*(coeffs(2)-1)/sqrt(vtau*(1-vtau));
    PERsist(i) = coeffs(2);
    Decision = [Decision; t_stat(i)>CV05(i)];
end
disp(' ');disp(' ');disp('t-statistics'); disp(' ');
disp(t_stat);  % t-statistic for different quantiles
disp('critical values'); disp(' ');
disp(CV05);    % critical values for different quantiles
disp('persistence parameters'); disp(' ');
disp(PERsist); % persistence parameter for different quantiles