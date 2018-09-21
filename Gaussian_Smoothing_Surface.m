[Swei] = Geodesic_Smoothing(Dist,type)

% The formula for a normal distribution Gaussian 
% is given by 
%                 1           ( (x-u)^2 )
%    f(r) = ------------ x exp| ------  |
%           sqrt(v*2*pi)      (   2v    )
% Where v is sigma, and u is the mean.

% simple Gaussian kernel sigma 1, mean 0, as for N pdf
if ~isempty(strcmp(lower(type)),'gaussian')

    FWHM = 4;
    sig = FWHM/sqrt(8*log(2));
    Swei = 1/sqrt(2*sig*pi) * exp(-Dist.^2/(2*sig^2));
else
end


