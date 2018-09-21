function varargout = Relabel_by_Size(varargin);
%
% Syntax :
%    outVector = Relabel_by_Size(inVector);
%
% This scripts relabels a vector assigning labels according to how many
% times an element in the input vector is repeated. More 
% repeated labels will relabel with lower label Id.
%
% Input Parameters:
%        inVector                   : Input Labels vector.
%
% Output Parameters:
%        outVector                  : Output Labels.
%
% See also: 
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ================================= Input parameters  ================================= %
if nargin < 1
    error('One Inputs is mandatory');
    return
elseif nargin >= 1
    inVector = varargin{1};
end

if nargout > 1
    error('To Many Output Parameters');
    return;
end

%% ========================== End of Input parameters  ================================= %
%% ================================ Main Program  ====================================== %

clustSizes = accumarray(inVector,inVector*0+1);
indmw = find(inVector == 4000); % Medial Wall Index
if indmw
    clustSizes(4000) = [];
end
[sortSizes,it] = sort(nonzeros(clustSizes),'descend');

[a,reLabels] = ismember(inVector,it);
if indmw
    reLabels(indmw) = 4000;
end

%% =============================== End of Main Program  ================================ %
% Outputs
varargout{1} = reLabels;
return;