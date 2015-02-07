function alpha=update_alpha(alphaold,a,b,Ncounts,Mcounts)

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

m=sum(sum(Mcounts));
[N,~]=size(Ncounts);
n=sum(Ncounts,2);
alp=(alphaold+1)*ones(N,1);
w=betarnd(alp,n);
prob=n./(alphaold+n);
bb=rand(N,1);
s=(bb<prob);
A=a+m-sum(s);
B=b-sum(log(w));
alpha=gamrnd(A,1/B);