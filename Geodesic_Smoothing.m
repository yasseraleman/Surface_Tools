function [Swei] = Geodesic_Smoothing(Dist,type)
%
% The formula for a normal distribution Gaussian 
% is given by 
%                 1           ( (x-u)^2 )
%    f(r) = ------------ x exp| ------  |
%           sqrt(v*2*pi)      (   2v    )
% Where v is sigma, and u is the mean.

% simple Gaussian kernel sigma 1, mean 0, as for N pdf
switch type
    case 'gaussian'
        FWHM = 4;
        sig = FWHM/sqrt(8*log(2));
        Swei =  exp(-Dist.^2/(2*sig^2));%Swei = Swei./normm(Swei(:));
    case 'nearest'
        ind = find(Dist == min(Dist));Swei = zeros(size(Dist));
        Swei(ind) = 1;
    case 'lineal'
       Swei= Dist;
end


function norma = normm(M)
norma = sqrt(sum((M').^2))';
return