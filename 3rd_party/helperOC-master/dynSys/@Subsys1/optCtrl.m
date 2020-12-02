function uOpt = optCtrl(obj, t, y, deriv, uMode)
%%%%%
% x(1) = theta
% x(2) = force
% x(3) = x
% x(4) = dx
% u(1) = theta_d
% u(2) = force_d
% mass : mass
% gains = [ktheta_p, kf_p] : inner controller gain
% cD : drag coefficient
%%%%%
% dx(1) = -ktheta_p * x(1) + ktheta_p * u(1) + D(1)
% dx(2) = -kf_p * x(2) + kf_p * u(2) + D(2)
% dx(3) = x(4)
% dx(4) = -cD/m * x(4) + 9.8*sin(x(1)) + x(2)/m*sin(x(1) + D(3)

if ~iscell(deriv)
    deriv = num2cell(deriv);
end

uOpt = cell(obj.nu, 1);

%% Optimal control
theta_ = y{1};
f_ = y{2};
if strcmp(uMode, 'max')
    if deriv{obj.dims==1} >= 0
        uOpt{1} = obj.uMax(1);
    else
        uOpt{1} = -obj.uMax(1);
    end
    
    if deriv{obj.dims==2} >= 0
        uOpt{2} = obj.uMax(2);
    else
        uOpt{2} = -obj.uMax(2);
    end
elseif strcmp(uMode, 'min')
    if deriv{obj.dims==1} >= 0
        uOpt{1} = -obj.uMax(1);
    else
        uOpt{1} = obj.uMax(1);
    end
    
    if deriv{obj.dims==2} >= 0
        uOpt{2} = -obj.uMax(2);
    else
        uOpt{2} = obj.uMax(2);
    end
else
    error('Unknown uMode!')
end

% %% Feedback Control
% x_ = y{1};
% y_ = y{2};
% theta_ = y{3};
%
% kx = 1.0;
% ky = 1.0;
% ktheta = 1.0;
%
% uOpt{1} = -kx * (x_ - obj.xD(1));
% uOpt{2} = -ky * (y_ - obj.xD(2));
% uOpt{3} = -ktheta * (theta_ - obj.xD(3));
end