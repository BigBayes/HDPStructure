# HDPStructure
    HDPStructure(X, dist, ...), simulates the HDPStructure
       posterior conditioned on X, and returns population assignments of
       the individuals in the population at each locus of interest.

       Inputs:
          X NxT. Allele matrix for the population. Here, X(i,j) denotes
                  the allelic value of individual i at locus j.

          dist Tx1. Vector indicating the distance between alleles, in
                  centimorgans/Mb. The value of dist(t), for t>1
                  indicates the distance between allele t-1 and allele t.
                  The value of dist(1) is unused.

          'PARAMETER', VALUE, ... Additional optional arguments must be
                  provided as additional arguments to the call to the
                  HDPStructure function in 'PARAMETER', VALUE, ... pairs.
                  In each pair, 'PARAMETER' specifies which parameter to
                  set, and VALUE specifies the value of that parameter.
                  The parameters recognised by HDPStructure are listed as
                  follows:

          initpop NxT. The initialization for the population assignment.
                  If you wish to start the MCMC chain at a population
                  assignment which is already close to a mode, then that
                  population assignment can be provided by the 'initpop'
                  parameter. If 'initpop' is not provided, then the
                  initialization is found using a linkage based method.

          alp0 scalar. We place a gamma prior on the concentration of the
                  top level Dirichlet process in the hierarchy. The
                  hyperparameters for that gamma prior are given by
                  'alp0' and 'bet0'. The default value of 'alp0' is 10.0.

          bet0 scalar. The default value of 'bet0' is 1.0.

          alp scalar. We also place a gamma prior on the concentrations
                  of each of the bottom level Dirichlet processes in the
                  hierarchy. The hyperparameters for this gamma prior are
                  given by 'alp' and 'bet'. The default value of 'alp' is
                  10.0.

          bet scalar. The default value of 'bet' is 5.0.

          a0init Tx1. The base measure of the top level Dirichlet process
                  of the hierarchy is parameterised by the locations of
                  interest. At each location, the base measure specifies
                  allele emission probabilities. We place beta priors on
                  these emission probabilities. In particular, at each
                  location l, the prior on the emission probabilities is
                  given by a bet a random variable with parameters a0(l)
                  and b0(l). We place a prior on each pair of variables
                  a0(l) and b0(l) as follows. First, we introduce
                  variables c and m, and we place a symmetric beta prior
                  on m with parameters 0.01, 0.01, and then we set
                  a0(l) = c * m, b0(l) = c * (1-m). We fix c = .01. By
                  default, the initial values of a0 and b0 are chosen by
                  matching the empirical allele distributions of the
                  data. In order to override this initialization,
                  'a0init' and 'b0init' may be provided. Note that these
                  are resampled using MH updates during each MCMC step
                  according to the priors specified above.

          b0init Tx1. The parameter 'b0init' is described above. By
                  default, 'b0init' is chosen to match the empirical
                  allele distribution of the data.

          lowL scalar. We place a log uniform prior on the latent
                  recombination rate. The range of this log uniform prior
                  is given by ['lowL', 'upL']. The default value of
                  'lowL' is 0.0.

          upL scalar. The default value of 'upL' is 50.0.

          iters scalar. The number of MCMC iterations to perform is given
                  by 'iters'. In each iteration, MCMC updates are
                  performed in turn for all of the latent variables of
                  the model. In particular, the latent variables are
                  updated in the following order, using the following
                  methods:

                      top level DP --- Gibbs
                      bottom level DPs --- MH
                      theta --- Gibbs
                      alpha0 --- MH
                      alpha --- Gibbs
                      recombination rate --- MH
                      a0 and b0 --- MH

                  The default value for 'iters' is 10000.

          burnin scalar. The number of MCMC iterations to discard at the
                  beginning of the chain before collecting the state of
                  the MCMC chain is given by 'burnin'. The default value
                  for 'burnin' is 5000.

          skip scalar. The number of MCMC iterations to thin after the
                  burnin is completed is given by 'skip'. In particular,
                  only every 'skip'-th iteration will be collected after
                  burnin. The default value for 'skip' is 5. The value of
                  skip must divide iters-burnin.

          seed scalar. For reproducibility, a seed for the random number
                  generator driving the MCMC may be provided. If the seed
                  is not provided, or if the seed provided is <=0, then
                  the state of matlab's random number generator is not
                  modified prior to the commencement of the MCMC.

       Outputs:
          POPULATION NxTxD. This array contains the population assignment
                  of the individuals in the population for each state in
                  the MCMC chain that is collected. The number of
                  collected states is given by D = (iters-burnin)/skip.
                  Therefore, the population for individual i at location
                  j specified by the state of the m-th iteration of MCMC
                  is given by POPULATION(i, j, m).

          RESULTS struct. Additional collected MCMC variables are given
                  in the structure RESULTS. These variables are as
                  follows:

                      RESULTS.INDSAMP
                      RESULTS.KSAMP
                      RESULTS.ALPHA0SAMP
                      RESULTS.ALPHASAMP
                      RESULTS.RECOMSAMP
                      RESULTS.MEANBASE
                      RESULTS.THETASAMP
                      RESULTS.POPNOMI

        Examples:

          > % Run HDPStructure on X, a binary array of alleles.
          > POPULATION=HDPStructure(X,dist);
          ...
          > % Display the population assignments of the last state in the
          > %     MCMC chain.
          > disp(POPULATION(:,:,end));

    Copyright (c) 2015, Maria De Iorio, Lloyd T. Elliott, Stefano Favaro
        and Yee Whye Teh.
