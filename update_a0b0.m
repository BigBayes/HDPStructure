function [a0,b0,mb]=update_a0b0(THETA,mbold)

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

[nSNP,~]=size(THETA);
mb=zeros(nSNP,1);
c=.01;
a0=zeros(nSNP,1);
b0=zeros(nSNP,1);
v=0.01;
for i=1:nSNP
    mb(i)=truncnorm(mbold(i),v,0,1);
    lpostnew=sum(log(betapdf(THETA(i,:),c*mb(i),c*(1-mb(i)))))+log(betapdf(mb(i),0.01,0.01));
    lpostold=sum(log(betapdf(THETA(i,:),c*mbold(i),c*(1-mbold(i)))))+log(betapdf(mbold(i),0.01,0.01));
    lprop=log(trunc_norm_ab_pdf( mbold(i), mb(i), v, 0,1 )) - log(trunc_norm_ab_pdf( mb(i),mbold(i) , v,0,1));
    acpt=exp(lpostnew-lpostold+lprop);
    if rand>acpt
        mb(i)=mbold(i);
    end
    a0(i)=c*mb(i);
    b0(i)=c*(1-mb(i));
end