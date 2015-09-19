function [Znew,Snew]=FFBS(weig,lambda,LTHETA1,LTHETA0,x)
%FFBS Forwards filtering-backwards sampling.
%   [Znew,Snew]=FFBS(weig,lambda,THETA1,THETA0,x), Propose new row of population 
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
PROB(1,:)=exp(log(pesi)+x(1)*LTHETA1(1,:)+(1-x(1))*LTHETA0(1,:));

%%%% Forwards filtering pass.
for t=2:T
    PROB(t-1,PROB(t-1,:)<10^-300)=10^-300;
    LPROB=log(PROB(t-1,:));
    lm1=-lambda(t)+LPROB + ...
        x(t)*LTHETA1(t,:)+(1-x(t))*LTHETA0(t,:);
    
    m1=exp(lm1);
    nc=logsumexp1(LPROB);
    lm0=nc+log((1-exp(-lambda(t)))) + ...
        log(pesi)+x(t)*LTHETA1(t,:)+(1-x(t))*LTHETA0(t,:);
    
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
s=rand<prob_s;
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
