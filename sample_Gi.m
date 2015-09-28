function [Znew,Snew,ninew,minew]=sample_Gi(xi,alpha,Zi_old,Si_old,pi0,ni,LTHETA1,LTHETA0,lambda)
%SAMPLE_GI Perform HM sampling for a bottom level DP in the hierarchy.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

T=length(xi);
KTOT=length(pi0);
K=KTOT-1;
nn=ni+alpha*pi0(1:K)';
alp=[nn (alpha*pi0(K+1))];
alp=reshape(alp,K+1,1);
pi=dirichrnd(alp);
Pold=unique(Zi_old);
piminold=min(pi(Pold));
Ci=unifrnd(0,piminold);
pitmp=pi(1:K);
l=pitmp>Ci;
pesi=pitmp(l);
[Ztmp,Snew]=FFBS(pesi,lambda,LTHETA1(:,l),LTHETA0(:,l),xi);

%%%% Change the k in Znew to the original numbering of the populations.
Znew=zeros(T,1);
populations=find(l>0);
for j=1:length(populations)
    l=Ztmp==j;
    Znew(l)=populations(j);
end
pnew2=unique(Znew);
piminnew=min(pi(pnew2));
lacp=log(piminold)-log(piminnew);
if rand>exp(lacp)
    Znew=Zi_old;
    Snew=Si_old;
end

%%%% Update the counts for the number of individuals in each population.
[ninew,minew]=update_counts(Znew,Snew,alpha,pi0,K);