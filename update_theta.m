function THETA=update_theta(X,Z,a0,b0)
%UPDATE_THETA Perform Gibbs update on the allele emission probabilities.
%   THETA=UPDATE_THETA(X,Z,a0,b0), returns a new vector of allele emission
%       probabilities given the current population assignments of the
%       individuals in the population, and also the parameters for the
%       prior on the allele emission probabilities.
%
% Inputs:
%           X NxT. Allele matrix for the population.
%
%           Z NxT. Population assignments for the individuals in the
%                   population. Here, Z(i,j) is the population assignment
%                   of the ith individual at the jth locus.
%
%           a0, b0 1xT. Prior pseudo observations for the allele emission
%                   probabilities. The prior allele emission probabilities
%                   at locus j follow the law of a Beta random variable
%                   with parameters a0(j) and b0(j).
%
% Outputs:
%           THETA. A draw from the conditional distriution of THETA
%                   conditioned on X, Z, a0 and b0.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.


[~,T]=size(X);
z=unique(Z);
K=length(z);
THETA=zeros(T,K);
for k=1:K
    apost=a0;
    bpost=b0;
    for i=1:T
        L=Z(:,i)==z(k);
        nsnp=sum(L);
        x=X(L,i);
        s=sum(x);
        apost(i)=apost(i)+s;
        bpost(i)=bpost(i)+nsnp-s;
    end
    THETA(:,k)=betarnd(apost,bpost);
end