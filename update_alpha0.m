function alpha0 = update_alpha0(Mcounts,a,b,alp0old)

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

n=sum(sum(Mcounts));
K=size(Mcounts,2);
x=betarnd(alp0old+1,n);
c1=a+K-1;
c2=n*(b-log(x));
w=c1/(c1+c2);
B=b-log(x);
if rand<w
    A=a+K;
else
    A=a+K-1;
end
alpha0=gamrnd(A,1/B);