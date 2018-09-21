function varargout = Surface_Rotation(Surf,varargin)
%
% Syntax :
%     Surfout = Surface_Rotation(Surf, angX, angY, angZ);
%
% This function rotates surface points around any axis in 3Dspace.
%
% Input Parameters:
%        Surf                   : Surface variable (file, struct or cellarray).
%
% Output Parameters:
%        Surfout                : Output Surface variable.
%
% See also: Plot_Surf Surf_Color Exp_Surf
%__________________________________________________
% Authors: Yasser Aleman Gomez 
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%

if nargin<1 % the indispensable input arguments are not provided
    error('At least one input is mandatory');
    return;
else
    % Surface Checking
    Surf = Surface_Checking(Surf);
    
    angX = 0;
    angY = 0;
    angZ = 0;
end


% deal with the input arguments
if nargin<1 % the indispensable input arguments are not provided
    error('Two inputs are mandatory');
    return;
else
    if numel(varargin)>0 % optional input arguments are provided
        while ~isempty(varargin)
            if numel(varargin)<2
                error('You need to provide optional input arguments as ''ParameterName''-''ParameterValue'' pairs.');
            end
            switch varargin{1}

                case 'angX'
                    % Limits in X axis (Left-Right axis)
                    angX=varargin{2};
                case 'angY'
                    % Limits in X axis (Anterior-Posterior axis)
                    angY=varargin{2};
                case 'angZ'
                    angZ=varargin{2};
                case 'rotMat'
                     rotMat=varargin{2};
                otherwise
                    error('Unexpected ''ParameterName'' input: %s\n',varargin{1});
            end
            varargin(1:2)=[]; % this pair of optional input arguments has been dealt with -- remove...
        end
    end
end

if nargout > 1
    error('To Many Output Parameters');
    return;
end


%% ========================= End of Checking Inputs ==================== %%

%% ======================== Main Program =============================== %%
if exist('rotMat','var')
    for i = 1:length(Surf)
        t = rotMat*[Surf(i).SurfData.vertices ones(size(Surf(i).SurfData.vertices,1),1)]';
        Surf(i).SurfData.vertices = t(1:3,:)';
    end
else
    angX = angX*pi/180;
    angY = angY*pi/180;
    angZ = angZ*pi/180;
    Rx = [1 0 0 0; 0 cos(angX) -sin(angX) 0; 0 sin(angX) cos(angX) 0; 0 0 0 1];
    Ry = [cos(angY) 0 sin(angY) 0; 0 1 0 0; -sin(angY) 0 cos(angY) 0; 0 0 0 1];
    Rz = [cos(angZ) -sin(angZ) 0 0; sin(angZ) cos(angZ) 0 0; 0 0 1 0; 0 0 0 1];
    Mat = Rx*Ry*Rz;
    for i = 1:length(Surf)
        
        t = [Surf(i).SurfData.vertices ones(size(Surf(i).SurfData.vertices,1),1)]*Mat';
        Surf(i).SurfData.vertices = t(:,1:3);
    end
end
%% ======================== End of Main Program ======================== %%

% Outputs
varargout{1} = Surf;

return;