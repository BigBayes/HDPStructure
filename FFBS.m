function [Znew,Snew]=FFBS(weig,lambda,THETA,x)
%FFBS Forwards filtering-backwards sampling.
%   [Znew,Snew]=FFBS(weig,lambda,THETA,x), Propose new row of population 
%       assignments for a given individual using the forwards filtering-
%       backwards sampling algorithm.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

K=length(weig);
pesi=weig';
T=length(x);
PROB=zeros(T,K);
M0=zeros(T,K);
M1=zeros(T,K);
tmp=THETA(1,:);

%%%% Snap to zero.
l=tmp<10^-300;
tmp(l)=10^-300;
l=tmp>.99999;
tmp(l)=.99999;
PROB(1,:)=exp(log(pesi)+x(1)*log(tmp)+(1-x(1))*log(1-tmp));

%%%% Forwards filtering pass.
for t=2:T
    l=PROB((t-1),:)<10^-300;
    PROB((t-1),l)=10^-300;
    tmp=THETA(t,:);
    l=tmp<10^-300;
    tmp(l)=10^-300;
    l=tmp>0.99999;
    tmp(l)=0.99999;
    lm1=-lambda(t)+log(PROB((t-1),:))+x(t)*log(tmp)+(1-x(t))*log(1-tmp);
    m1=exp(lm1);
    lp=log(PROB((t-1),:));
    mp=max(lp);
    nc=mp+log(sum(exp(lp-mp)));
    lm0=nc+log((1-exp(-lambda(t)))) + ...
        log(pesi)+x(t)*log(tmp)+(1-x(t))*log(1-tmp);
    
    m0=exp(lm0);
    mp=max(lm0,lm1);
    temp=mp+log(exp(lm0-mp)+exp(lm1-mp));
    PROB(t,:)=exp(temp);
    M1(t,:)=m1;
    M0(t,:)=m0;
end

%%%% Backwards sampling pass.
Znew=zeros(T,1);
Snew=zeros(T,1);
k=discrnd(PROB(T,:)/sum(PROB(T,:)));
Znew(T)=k;
lprob_s=log(M1(T,k))-log(PROB(T,k));
prob_s=exp(lprob_s);
s=rand< prob_s;
Snew(T)=s;
for i=1:(T-1)
    ind=T-i;
    if s==1
        k=Znew(ind+1);
    else
        k=discrnd(PROB(ind,:)/sum(PROB(ind,:)));
    end
    Znew(ind)=k;
    if ind>1
        lprob_s=log(M1(ind,k))-log(PROB(ind,k));
        prob_s=exp(lprob_s);
        s=rand<prob_s;
    end
    Snew(ind)=s;
end
Snew(1)=0;
