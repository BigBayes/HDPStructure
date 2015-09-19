function [POPULATION,RESULTS]=HDPStructure(X,dist,varargin)
%HDPSTRUCTURE Inference for the Bayesian nonparametric structure.
%   POPULATIONS = HDPSTRUCTURE(X, dist, ...), simulates the HDPSTRUCTURE
%       posterior conditioned on X, and returns population assignments of
%       the individuals in the population at each locus of interest.
%
% Inputs:
%           X NxT. Allele matrix for the population. Here, X(i,j) denotes
%                   the allelic value of individual i at locus j.
%
%           dist Tx1. Vector indicating the distance between alleles, in
%                   centimorgans/Mb. The value of dist(t), for t>1
%                   indicates the distance between allele t-1 and allele t.
%                   The value of dist(1) is unused.
%
%           'PARAMETER', VALUE, ... Additional optional arguments must be
%                   provided as additional arguments to the call to the
%                   HDPStructure function in 'PARAMETER', VALUE, ... pairs.
%                   In each pair, 'PARAMETER' specifies which parameter to
%                   set, and VALUE specifies the value of that parameter.
%                   The parameters recognised by HDPStructure are listed as
%                   follows:
%
%           initpop NxT. The initialization for the population assignment.
%                   If you wish to start the MCMC chain at a population
%                   assignment which is already close to a mode, then that
%                   population assignment can be provided by the 'initpop'
%                   parameter. If 'initpop' is not provided, then the
%                   initialization is found using a linkage based method.
%
%           alp0 scalar. We place a gamma prior on the concentration of the
%                   top level Dirichlet process in the hierarchy. The
%                   hyperparameters for that gamma prior are given by
%                   'alp0' and 'bet0'. The default value of 'alp0' is 10.0.
%
%           bet0 scalar. The default value of 'bet0' is 1.0.
%
%           alp scalar. We also place a gamma prior on the concentrations
%                   of each of the bottom level Dirichlet processes in the
%                   hierarchy. The hyperparameters for this gamma prior are
%                   given by 'alp' and 'bet'. The default value of 'alp' is
%                   10.0. 
%
%           bet scalar. The default value of 'bet' is 5.0.
%
%           a0init Tx1. The base measure of the top level Dirichlet process
%                   of the hierarchy is parameterised by the locations of
%                   interest. At each location, the base measure specifies
%                   allele emission probabilities. We place beta priors on
%                   these emission probabilities. In particular, at each
%                   location l, the prior on the emission probabilities is
%                   given by a bet a random variable with parameters a0(l)
%                   and b0(l). We place a prior on each pair of variables
%                   a0(l) and b0(l) as follows. First, we introduce
%                   variables c and m, and we place a symmetric beta prior
%                   on m with parameters 0.01, 0.01, and then we set
%                   a0(l) = c * m, b0(l) = c * (1-m). We fix c = .01. By
%                   default, the initial values of a0 and b0 are chosen by
%                   matching the empirical allele distributions of the
%                   data. In order to override this initialization,
%                   'a0init' and 'b0init' may be provided. Note that these
%                   are resampled using MH updates during each MCMC step
%                   according to the priors specified above. 
%
%           b0init Tx1. The parameter 'b0init' is described above. By
%                   default, 'b0init' is chosen to match the empirical
%                   allele distribution of the data.
%
%           lowL scalar. We place a log uniform prior on the latent
%                   recombination rate. The range of this log uniform prior
%                   is given by ['lowL', 'upL']. The default value of
%                   'lowL' is 0.0.
%
%           upL scalar. The default value of 'upL' is 50.0.
%
%           iters scalar. The number of MCMC iterations to perform is given
%                   by 'iters'. In each iteration, MCMC updates are
%                   performed in turn for all of the latent variables of
%                   the model. In particular, the latent variables are
%                   updated in the following order, using the following
%                   methods:
%
%                       top level DP --- Gibbs
%                       bottom level DPs --- MH
%                       theta --- Gibbs
%                       alpha0 --- MH
%                       alpha --- Gibbs
%                       recombination rate --- MH
%                       a0 and b0 --- MH
%
%                   The default value for 'iters' is 10000.
%
%           burnin scalar. The number of MCMC iterations to discard at the
%                   beginning of the chain before collecting the state of
%                   the MCMC chain is given by 'burnin'. The default value
%                   for 'burnin' is 5000.
%
%           skip scalar. The number of MCMC iterations to thin after the
%                   burnin is completed is given by 'skip'. In particular,
%                   only every 'skip'-th iteration will be collected after
%                   burnin. The default value for 'skip' is 5. The value of
%                   skip must divide iters-burnin.  The value of burnin
%                   must be less than the value of iters --- the total
%                   number of samples used for the posterior estimates is
%                   (iters-burnin)/skip.
%
%           seed scalar. For reproducibility, a seed for the random number
%                   generator driving the MCMC may be provided. If the seed
%                   is not provided, or if the seed provided is <=0, then
%                   the state of matlab's random number generator is not
%                   modified prior to the commencement of the MCMC.
%
%           CORES scalar. Number of cores to use in parallel. Defaults to
%                   one core. The number of cores should be less than or
%                   equal to N.
%
% Outputs:
%           POPULATION NxTxD. This array contains the population assignment
%                   of the individuals in the population for each state in
%                   the MCMC chain that is collected. The number of
%                   collected states is given by D = (iters-burnin)/skip.
%                   Therefore, the population for individual i at location
%                   j specified by the state of the m-th iteration of MCMC
%                   is given by POPULATION(i, j, m).
%
%           RESULTS struct. Additional collected MCMC variables are given
%                   in the structure RESULTS. These variables are as
%                   follows:
%
%                       RESULTS.INDSAMP
%                       RESULTS.KSAMP
%                       RESULTS.ALPHA0SAMP
%                       RESULTS.ALPHASAMP
%                       RESULTS.RECOMSAMP
%                       RESULTS.MEANBASE
%                       RESULTS.THETASAMP
%                       RESULTS.POPNOMI
%
% Examples:
%
%           > % Run HDPStructure on X, a binary array of alleles.
%           > POPULATION=HDPStructure(X,dist);
%           ...
%           > % Display the population assignments of the last state in the
%           > %     MCMC chain.
%           > disp(POPULATION(:,:,end));

%   Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
%       and Yee Whye Teh.

[n,T]=size(X);
parser=inputParser;
default_initpop=[];
default_alp0=10.0;
default_bet0=1.0;
default_alp=10.0;
default_bet=5.0;
default_a0init=[];
default_b0init=[];
default_lowL=0.0;
default_upL=50.0;
default_iters=10000;
default_burnin=5000;
default_skip=5;
default_seed=0;
default_cores=1;
addRequired(parser,'X',@ismatrix);
addRequired(parser,'dist',@isvector);
addOptional(parser,'initpop',default_initpop);
addOptional(parser,'alp0',default_alp0,@isnumeric);
addOptional(parser,'bet0',default_bet0,@isnumeric);
addOptional(parser,'alp',default_alp,@isnumeric);
addOptional(parser,'bet',default_bet,@isnumeric);
addOptional(parser,'a0init',default_a0init,@ismatrix);
addOptional(parser,'b0init',default_b0init,@ismatrix);
addOptional(parser,'lowL',default_lowL,@isnumeric);
addOptional(parser,'upL',default_upL,@isnumeric);
addOptional(parser,'iters',default_iters,@isnumeric);
addOptional(parser,'burnin',default_burnin,@isnumeric);
addOptional(parser,'skip',default_skip,@isnumeric);
addOptional(parser,'seed',default_seed,@isnumeric);
addOptional(parser,'cores',default_cores,@isnumeric);
parse(parser,X,dist,varargin{:});
initpop=parser.Results.initpop;
alp0=parser.Results.alp0;
bet0=parser.Results.bet0;
alp=parser.Results.alp;
bet=parser.Results.bet;
a0init=parser.Results.a0init;
b0init=parser.Results.b0init;

try
    logsumexp1([1 2 3 4]);
catch ME %#ok
    cwd = pwd();
    cd(mfilename('fullpath'));
    mex logsumexp1.c;
    cd(cwd);
end

% Initialize a0 and b0 from empirical measures, if they aren't provided.
if isempty(a0init)
    pm=mean(X);
    a0=pm.*(pm.*(1-pm)./0.1-1);
    l=find(a0<=0);
    a0(l)=.1; %#ok
    a0init=reshape(a0,T,1);
end

if isempty(b0init)
    pm=mean(X);
    a0=pm.*(pm.*(1-pm)./0.1-1);
    l=find(a0<=0);
    b0=a0.*(1-pm)./pm;
    b0(l)=.1; %#ok
    b0init=reshape(b0,T,1);
end

lowL=parser.Results.lowL;
upL=parser.Results.upL;
iters=parser.Results.iters;
burnin=parser.Results.burnin;
skip=parser.Results.skip;


if burnin >= iters
    error('Number of burnin samples should be <= number of iterations.');
end

if mod(iters, skip) ~= 0
    error('Number of skipped (thinned) samples should divide number of iterations.');
end

seed=parser.Results.seed;
if seed>0
    s=RandStream('mt19937ar','Seed',seed);
    RandStream.setGlobalStream(s);
end

cores=parser.Results.cores;

if cores > 1
    try
        parpool(cores);
    catch ME %#ok
        delete(gcp);
        parpool(cores);
    end
end

%%%% Initialize first sample --- if initpop is provided, copy it.
if ~isempty(initpop)
    initclus=initpop;
    K=length(unique(initpop(:)));
else
    [initclus,K]=clusinitial(X);
end

Zold=repmat(initclus,[1,T]);
Ztmp=Zold;

%%%% Initialize parameters.
alpha0=1.5;
alpha=.2;

%%%% Initialization of linkage indicators.
Sold=zeros(n,T);
for i=2:T
    for j=1:n
        if Ztmp(j,i)==Ztmp(j,i-1)
            Sold(j,i)=1;
        end
    end
end

KMAX=10;

%%%% Initialization of table counts.
n_counts=nan(n,K);
m_counts=nan(n,K);
for i =1:n
    [ni,mi]=update_counts(Zold(i,:),Sold(i,:),alpha,ones(K,1)./K,K);
    n_counts(i,:)=ni;
    m_counts(i,:)=mi;
end

a0=a0init;
b0=b0init;
mbold=a0./(a0+b0);
THETA=zeros(T,K);
for i =1:K
    for j=1:T
        l=Ztmp(:,j)==i;
        THETA(j,i)=mean(X(l,j));
    end
end

rr=.1;
lambda=rr*dist;
KTOTAL=K;
POPLABEL=1:K;

%%%% Initialization of output variables.
ALPHA0SAMP=zeros((iters-burnin)/skip,1);
ALPHASAMP=zeros((iters-burnin)/skip,1);
KSAMP=zeros((iters-burnin)/skip,1);
RECOMSAMP=zeros((iters-burnin)/skip,1);
POPULATION=zeros(n,T,(iters-burnin)/skip);
INDSAMP=zeros(n,T,(iters-burnin)/skip);
MEANBASE=zeros(T,(iters-burnin)/skip);
POPNOMI=-ones((iters-burnin+1),50);
indice=0;

%%%% Main loop for MCMC.
for iter = 1:iters
    if mod(iter,100)==0
        fprintf(1, 'Iteration %d/%d\n', iter, iters);
    end
    
    POPLABEL=[POPLABEL (KTOTAL+1):(KTOTAL+KMAX)]; %#ok
    KTOTAL=KTOTAL+KMAX;
    pi0=sample_G0(alpha0,m_counts,K,KMAX);
    for i=1:KMAX
        THETA =[THETA betarnd(a0,b0)]; %#ok
    end
    
    if iter>(burnin+1)
        THETASAMP=[THETASAMP zeros(T,KMAX)]; %#ok
        THETACOUNT=[THETACOUNT zeros(1,KMAX)]; %#ok
    end
    
    NTOT=[n_counts zeros(n,KMAX)];
    MTOT=[m_counts zeros(n,KMAX)];
    if length(POPLABEL)~=(KMAX+K)
        error('Mismatch between actual and expected number of populations after padding.');
    end
    
    % Moved outside FFBS
    tmp=THETA;
    l=tmp<10^-300;
    tmp(l)=10^-300;
    l=tmp>.99999;
    tmp(l)=.99999;
    tmp1 = log(tmp);
    tmp0 = log(1-tmp);
    if cores == 1
        for i=1:n
            [Z,S,ninew,minew]=sample_Gi(X(i,:),alpha,Ztmp(i,:),Sold(i,:),pi0,NTOT(i,:),tmp1,tmp0,lambda);
            Ztmp(i,:)=Z';
            Sold(i,:)=S';
            NTOT(i,:)=ninew;
            MTOT(i,:)=minew;
        end
    else
        parfor i=1:n
            [Z,S,ninew,minew]=sample_Gi(X(i,:),alpha,Ztmp(i,:),Sold(i,:),pi0,NTOT(i,:),tmp1,tmp0,lambda);
            Ztmp(i,:)=Z';
            Sold(i,:)=S';
            NTOT(i,:)=ninew;
            MTOT(i,:)=minew;
        end
    end
    
    Zold=Ztmp;
    for i=1:size(MTOT,2)
        l=Ztmp==i;
        
        %%%% Qui e quello giusto.
        Zold(l)=POPLABEL(i);
    end
    
    l=find(sum(MTOT)==0);
    sn=sum(NTOT);
    if any(sn(l(:))>0)
        error('Unpruned empty populations detected.');
    end
    j=find(sum(MTOT)>0);
    K=length(j);
    m_counts=MTOT(:,j);
    n_counts=NTOT(:,j);
    POPLABEL=POPLABEL(j);
    
    if iter>(burnin+1)
        THETASAMP=THETASAMP(:,j);
        THETACOUNT=THETACOUNT(:,j);
        POPNOMI((iter-burnin+1),1:K)=POPLABEL;
    end
    
    if size(m_counts,2)~=K
        error('Mismatch between actual and expected number of populations after pruning.');
    end
    
    for i=1:K
        l=Zold==POPLABEL(i);
        Ztmp(l)=i;
    end
    
    %%%% MCMC update for theta.
    THETA=update_theta(X,Zold,a0,b0);
    if iter==(burnin+1)
        THETASAMP=THETA;
        pktmp=length(POPLABEL);
        POPNOMI(1,1:pktmp)=POPLABEL;
        THETACOUNT=ones(1,pktmp);
    end
    
    if iter>(burnin+1)
        THETASAMP=THETASAMP+THETA;
        THETACOUNT=THETACOUNT+1;
    end
    
    alp0old=alpha0;
    alpha0=update_alpha0(m_counts,alp0,bet0,alp0old);
    alphaold=alpha;
    alpha=update_alpha(alphaold,alp,bet,n_counts,m_counts);
    
    %%%% Update recombination rate.
    rrold=rr;
    rr=update_recomb(dist,Sold,rrold,lowL,upL);
    lambda=rr*dist;
    
    %%%% Update hyperparameter for HDP base measure.
    [a0,b0,mb]=update_a0b0(THETA,mbold);
    mbold=mb;
    
    %%%% Save output.
    if (mod(iter,skip)==0)&&(iter>burnin)
        indice=indice+1;
        ALPHA0SAMP(indice)=alpha0;
        ALPHASAMP(indice)=alpha;
        KSAMP(indice)=K;
        RECOMSAMP(indice)=rr;
        POPULATION(:,:,indice)=Zold;
        INDSAMP(:,:,indice)=Sold;
        MEANBASE(:,indice)=mbold;
    end
end

for i=1:K
    THETASAMP(:,i)=THETASAMP(:,i)/THETACOUNT(i);
end

%%%% Save results of MCMC sampling in a structure.
RESULTS=struct( ...
    'INDSAMP', INDSAMP, ...
    'KSAMP', KSAMP, ...
    'ALPHA0SAMP', ALPHA0SAMP, ...
    'ALPHASAMP', ALPHASAMP, ...
    'RECOMSAMP', RECOMSAMP, ...
    'MEANBASE', MEANBASE, ...
    'THETASAMP', THETASAMP, ...
    'POPNOMI', POPNOMI ...
);



if cores > 1
    delete(gcp);
end