function pdf = trunc_norm_ab_pdf( x, mu, s2, a, b )

%*****************************************************************************
%
%% 
%
%  Parameters:
%
%    Input, real X, the argument of the PDF.
%
%    Input, real MU, SIGMA, the mean and standard deviation of the
%    parent Normal distribution.
%
%    Input, real A, B, the lower and upper truncation limits.
%
%    Output, real PDF, the value of the PDF.
%
  sigma=sqrt(s2);  
  alpha = ( a - mu ) / sigma;
  beta = ( b - mu ) / sigma;
  xi = ( x - mu ) / sigma;

  alpha_cdf = normcdf( alpha );
  beta_cdf = normcdf( beta );
  xi_pdf = normpdf( xi );

  pdf = xi_pdf / ( beta_cdf - alpha_cdf ) / sigma;

  return
end
