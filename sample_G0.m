function pi0=sample_G0(alpha0,m_count,K,KMAX)
%SAMPLE_GI Perform Gibbs update for the top level DP in the hierarchy.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

m=sum(m_count);
alp=[m alpha0];
alp=reshape(alp,K+1,1);
pi0=dirichrnd(alp);
w0=1;
p0tmp=zeros(KMAX,1);
for i=1:KMAX
    b0=betarnd(1,alpha0);
    p0new=w0*b0;
    p0tmp(i)=p0new;
    w0=w0*(1-b0);
end

p0star=pi0(K+1)*w0;
p0tmp=p0tmp*pi0(K+1);
pi0=[pi0(1:K);p0tmp;p0star];