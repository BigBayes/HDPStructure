function [ni,mi]=update_counts(Znew,Snew,alpha,weigh,K)

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

ni=zeros(1,K);
mi=zeros(1,K);
for k =1:K
    l=(Znew==k)&(Snew==0);
    tot=sum(l);
    ni(k)=tot;
    if tot>0
        w=alpha*weigh(k);
        c=1:tot;
        pp=w./(w+c-1);
        bb=rand(1,length(pp));
        mi(k)=sum(bb<pp);
    end
end

