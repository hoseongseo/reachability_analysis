function uOpt = optCtrl(obj, t, y, deriv, uMode)
% uOpt = optCtrl(obj, t, y, deriv, uMode)

if ~iscell(deriv)
    deriv = num2cell(deriv);
end

uOpt = cell(obj.nu, 1);

% %% Optimal control
% if strcmp(uMode, 'max')
%     uOpt{1} = (deriv{obj.dims==1}>=0)*obj.uMax(1) + (deriv{obj.dims==1}<0)*(-obj.uMax(1));
%     uOpt{2} = (deriv{obj.dims==2}>=0)*obj.uMax(2) + (deriv{obj.dims==2}<0)*(-obj.uMax(2));
%     uOpt{3} = (deriv{obj.dims==3}>=0)*obj.uMax(3) + (deriv{obj.dims==3}<0)*(-obj.uMax(3));
% elseif strcmp(uMode, 'min')
%     uOpt{1} = (deriv{obj.dims==1}>=0)*(-obj.uMax(1)) + (deriv{obj.dims==1}<0)*obj.uMax(1);
%     uOpt{2} = (deriv{obj.dims==2}>=0)*(-obj.uMax(2)) + (deriv{obj.dims==2}<0)*obj.uMax(2);
%     uOpt{3} = (deriv{obj.dims==3}>=0)*(-obj.uMax(3)) + (deriv{obj.dims==3}<0)*obj.uMax(3);
% else
%     error('Unknown uMode!')
% end

%% Feedback Control

x_ = y{1};
y_ = y{2};
theta_ = y{3};

kx = 1.0;
ky = 1.0;
ktheta = 1.0;

uOpt{1} = -kx * (x_ - obj.xD(1));
uOpt{2} = -ky * (y_ - obj.xD(2));
uOpt{3} = -ktheta * (theta_ - obj.xD(3));
end