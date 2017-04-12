function [initclus,C,mm] = clusinitial(X)
%CLUSINITIAL Create initial population assignment based on linkage.
%   [initclus,C] = CLUSINITIAL(alpha), returns an initial population
%       assignment of the individuals in the population based on the
%       genetic distance between them. This population assignment (which is
%       not location specific) is used to initialize the MCMC. The
%       population assignment is chosen by looping through MAXCLUS possible
%       total numbers of populations and then forming the population
%       assignment for that number of populations using the 'linkage'
%       function. Then, the population assignment corresponding to the
%       most likely total number of populations is returned.
%
% Inputs:
%           X NxD. Matrix of alleles for N individuals and D loci. Value of
%                   X(i,j) indicates the dosage of the major allele for 
%                   individual i at location j.
%
% Outputs:
%           initclus Nx1. Matrix of population assignments, one for each
%                   individual.
%
%           C scalar. The number of populations appearing in the
%                   assignment.
%
%           mm scalar. The likelihood of the population assignment.

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

[N,M]=size(X);
MAXCLUS=4;
ll=zeros(MAXCLUS,1);
clus=ones(N,MAXCLUS);
for i=2:MAXCLUS
    d = pdist(X);
    z = linkage(d, 'complete');
    c = cluster(z,'maxclust',i);
    
    %%%% Save the population assignment.
    clus(:,i)=c;
end

%%%% Evaluate likelihood.
freq=mean(X,1);
s=freq<0.000000001;
freq(s)=0.000000001;
s=freq==1;
freq(s)=0.999999;
for j=1:M
    ll(1)=ll(1)+sum(X(:,j)*log(freq(j))+(1-X(:,j))*log(1-freq(j)));
end

for k=2:MAXCLUS
    for i=1:k
        ind=clus(:,k)==i;
        Xtmp=X(ind,:);
        freq=mean(Xtmp,1);
        s=freq<0.000000001;
        freq(s)=0.000000001;
        s=freq==1;
        freq(s)=0.999999;
        for j=1:M
            ll(k)=ll(k)+sum(Xtmp(:,j)*log(freq(j))+(1-Xtmp(:,j))*log(1-freq(j)));
        end
    end
end

[mm,C]=max(ll);
initclus=clus(:,C);
