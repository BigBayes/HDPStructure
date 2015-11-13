function lam=update_recomb(dist,S,lamold,lowL,upL)

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

v=0.1;

%%%% Prior on log(lam) is uniform on the range [lowL, upL]. Posterior\
%%%% density of lam is 1/(y*(upL-lowL))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% ADAPTATION
%u=rand;
%if u < 0.05
 %   ll=truncnorm(log(lamold),v,lowL,upL);
%else
 %  ll=truncnorm(log(lamold),varll*(2.38^2) ,lowL,upL);
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ll=truncnorm(log(lamold),v ,lowL,upL);
lam=exp(ll);
lpriorold=-log(lamold);
lpriornew=-log(lam);
lprop=log(trunc_norm_ab_pdf( log(lamold), ll, v, lowL,upL )) - log(trunc_norm_ab_pdf( ll,log(lamold) , v, lowL,upL ));
[~,T]=size(S);
llikenew=0;
llikeold=0;
for i=2:T
    pold=exp(-lamold*dist(i));
    l=pold<10^-300;
    pold(l)=10^-300;
    pnew=exp(-lam*dist(i));
    l=pnew<10^-300;
    pnew(l)=10^-300;
    temp=S(:,i).*log(pnew)+(1-S(:,i)).*log((1-pnew));
    llikenew=llikenew+sum(temp);
    temp=S(:,i).*log(pold)+(1-S(:,i)).*log((1-pold));
    llikeold=llikeold+sum(temp);
end
lpostnew=lpriornew+llikenew;
lpostold=lpriorold+llikeold;
lacpt=lpostnew-lpostold+lprop;
if rand>exp(lacpt)
    lam=lamold;
end






