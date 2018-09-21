function [varargout] = Surf_Color(varargin);
% Syntax :
% [Colors] = Surf_Color(Surf,cl, N, boolflip);
%
% This function plots the surfaces contained in the cell file Surfa.
%
% Input Parameters:
%   Surf        : Surface variable.
%   cl          : Colormap used to see the results.
%   N           : Useful percent in the colormap
%
% Output Parameters:
%  Colors       : Output colormap matrix.
%
% Related references:
%
%
% See also: Smooth_Surf Surf_Comp Plot_Surf Plot_oversurf Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez
% Neuroimaging Department
% Cuban Neuroscience Center
% September 29th 2007
% Version $1.0


if nargin == 0
    error('Few input parameters');
    return;
elseif nargin == 1
    Surf = varargin{1};
    cl = 'jet';
    N = [1 100];
    boolflip = 0;
elseif nargin == 2
    Surf = varargin{1};
    cl = varargin{2};
    N = [1 100];
    boolflip = 0;
elseif nargin == 3
    Surf = varargin{1};
    cl = varargin{2};
    N = varargin{3};
    boolflip = 0;
elseif nargin == 4
    Surf = varargin{1};
    cl = varargin{2};
    N = varargin{3};
    boolflip = varargin{4};
elseif nargin > 4
    error('To Many Input Parameters');
    return;
end
if nargout > 1
    error('To Many Output Parameters');
    return;
end
if length(N) == 1
    N = [1 N];
end

%% ================= End of checking Input parameters ================== %%

%% ========================= Main program================================ %%
tempVar = whos('cl');
txt=Surf.Is;
if sum(txt-floor(txt))
    [ut,i1t,it2] = unique(txt);
    cent_col=(0-min(ut(:)))/(max(ut(:))-min(ut(:)));
    
    switch tempVar.class;
        case {'double','single'};
            col = cl;
        case 'char';
            [col] = colormaps_colors(cl,size(ut,1),cent_col);
    end
    
    if length(ut) > 1
        Number = [floor(size(col,1)*N(1)/100) floor(size(col,1)*N(2)/100)];
        if Number(1) == 0
            Number(1) = 1;
        end
        tempc = col(Number(1):Number(2),:);
        if boolflip
            tempc = flipdim(tempc,1);
        end
        
        m = size(tempc,1);
        X0 = linspace (1, size(tempc,1), size(ut,1));
        tempc = interp1(1:m,tempc,X0);
        
        Y0 = cumsum([0;diff(ut)/sum(diff(ut))]);
        X1 = linspace(0, 1, size(tempc,1));
        a = 1;b = -1;c = 0; % Line parameters
        Xo = X1(:);Yo = Y0(:);
        x = (b*(b*Xo-a*Yo) -a*c)./(a.^2+b.^2);
        y = (a*(-b*Xo+a*Yo) -b*c)./(a.^2+b.^2);
        
        ind = find(x >max(Xo));
        x(ind) = max(Xo);
        ind = find(x <min(Xo));
        x(ind) = min(Xo);
        
        ind = find(y >max(Yo));
        y(ind) = max(Yo);
        ind = find(y <min(Yo));
        y(ind) = min(Yo);
        col = interp1(X1,tempc,x);
        Colors = abs(col(it2,:));
    else
        Colors = repmat(col, [length(it2) 1]);
    end
elseif sum(txt-floor(txt)) ==0
    %col = [213 221 227;175 206 227;149 196 228; 120 186 230;87 172 231;24 146 232;6 73 119;244 207 154;244 192 117;244 179 86;244 161 43; 212 133 20;158 101 19; 113 71 12]/255;
    
    col = [[1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];[240,163,255;0,117,220;153,63,0;76,0,92;0,92,49;43,206,72;...
        255,204,153;148,255,181;143,124,0;157,204,0;194,0,136;...
        0,51,128;255,164,5;255,168,187;66,102,0;255,0,16;94,241,242;0,153,143;...
        224,255,102;116,10,255;153,0,0;255,255,128;255,255,0;255,80,5]/255];
    
%     col = [1 0 0;0 1 0; 0 0 1; 1 1 0;0 1 1;1 0 1;1 0.5 0;0 0.5 1;1 0 0.5;0.5 1 0; 0.25 0.5 1;0.6 0.3 0.52;0.7 0.5 0.9];
    ind = find(txt == 0);
    txt(ind) = max(txt)+1;
    Ncolor = size(col,1);
    ut = sort(unique(txt));
    re = floor(length(ut)/Ncolor); col = repmat(col,[re+1 1]);
    for j = 1:size(ut,1)
        indpos = find(txt == ut(j)); Colors(indpos,:) = repmat(col(j,:),[size(indpos,1) 1]);
    end
    
    %Creating boundary lines
    %     Is=txt;Is(Is==0)=max(Is)+1;
    %     St=unique(Is);
    %     Ns =size(St,1);
    %     neibb = size(Surf.SurfData.vertices,1)+1;
    %     Ist = 0*Is;
    %     for i =1:Ns
    %         ind = find(Is == St(i));
    %         for j=1:size(ind,1)
    %             Neigh=Surf.SurfData.faces(nonzeros(Surf.Tri(ind(j),3:end)),:);Neigh=unique(Neigh(:));Neigh(Neigh==ind(j))=[];
    %             c=accumarray(double(Is(Neigh)),ones(size(Is(Neigh),1),1));
    %             if size(nonzeros(c),1)~=1&(~logical(sum(ismember(Neigh,neibb))))
    %                 Ist(ind(j))=Is(ind(j));
    %             end
    %         end
    %         neibb = unique([neibb;find(Ist~=0)]);
    %     end
    %     indl =find(Ist~=0);
    %     Colors(indl,:) =repmat([0    0    0],[size(indl,1) 1]);
    Colors(ind,:) = repmat([1 1 1],[length(ind) 1]);
    
    
    
end
varargout{1} = Colors;
%========================End of main program==============================%
return

function s = spectral(m)
%SPECTRAL Black-purple-blue-green-yellow-red-white color map.
%
%         map = spectral(num_colors)
%
% SPECTRAL(M) returns an M-by-3 matrix containing a "spectral" colormap.
% SPECTRAL, by itself, is the same length as the current colormap.
%
% For example, to reset the colormap of the current figure:
%
%           colormap(spectral)
%
% See also HSV, GRAY, PINK, HOT, COOL, BONE, COPPER, FLAG,
%          COLORMAP, RGBPLOT.

if nargin < 1, m = size(get(gcf,'colormap'),1); end
base = [
    0.0000 0.0000 0.0000
    0.4667 0.0000 0.5333
    0.5333 0.0000 0.6000
    0.0000 0.0000 0.6667
    0.0000 0.0000 0.8667
    0.0000 0.4667 0.8667
    0.0000 0.6000 0.8667
    0.0000 0.6667 0.6667
    0.0000 0.6667 0.5333
    0.0000 0.6000 0.0000
    0.0000 0.7333 0.0000
    0.0000 0.8667 0.0000
    0.0000 1.0000 0.0000
    0.7333 1.0000 0.0000
    0.9333 0.9333 0.0000
    1.0000 0.8000 0.0000
    1.0000 0.6000 0.0000
    1.0000 0.0000 0.0000
    0.8667 0.0000 0.0000
    0.8000 0.0000 0.0000
    0.8000 0.8000 0.8000
    ];
n = length(base);
X0 = linspace (1, n, m);
s = interp1(1:n,base,X0);

return