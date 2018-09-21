function [At] = Tri_Area_Comp(fv);
%
% Syntax :
% [At] = Area_Comp(fv);
%
% This function computes the surface area for patch fv.
%
% Input Parameters:
%   fv       : Surface Patch.
%
% Output Parameters:
%   At       : Surfaces Area.
%
% Related references:
%
%
% See also: Red_Surf Smooth_Surf Plot_Surf Surf_Comp Atlas_Surf
% Plot_oversurf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% February 17th 2007
% Version $1.0

%=========================Main program====================================%
At = 0;
N = size(fv.faces,1);
for i = 1:N;
    di = dist(fv.vertices(fv.faces(i,:),:)');
    A = abs(di(1,2));
    B = abs(di(1,3));
    C = abs(di(2,3));
    p = (A+B+C)/2;
    Ar = sqrt(p*(p-A)*(p-B)*(p-C));
    %Ar = (A*B/2)*(sqrt(1-((A^2+B^2-C^2)^2)/((2*A*B)^2)));
    At(i) = Ar;
end
At = real(At)/100;% cm^2
%========================End of main program==============================%
return