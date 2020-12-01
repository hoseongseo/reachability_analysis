function uOpt = optCtrl(obj, t, y, deriv, uMode)
%%%%%
% x(1) = x
% x(2) = dx
% x(3) = theta
% u(1) = theta_d
% mass : mass
% gains = ktheta_p : inner controller gain
% cD : drag coefficient
%%%%%
% dx(1) = x(2)
% dx(2) = -cD/m * x(2) + 9.8*sin(x(3)) + D(1)
% dx(3) = -ktheta_p * x(3) + ktheta_p * u(1) + D(2)


if ~iscell(deriv)
    deriv = num2cell(deriv);
end

uOpt = cell(obj.nu, 1);

%% Optimal control

if strcmp(uMode, 'max')
    if deriv{obj.dims==3} > 0
        uOpt{1} = obj.uMax(1);
    elseif deriv{obj.dims==3} < 0
        uOpt{1} = -obj.uMax(1);
    else
        uOpt{1} = 0.0;
    end

elseif strcmp(uMode, 'min')
    if deriv{obj.dims==3} > 0
        uOpt{1} = -obj.uMax(1);
    elseif deriv{obj.dims==3} < 0
        uOpt{1} = obj.uMax(1);
    else
        uOpt{1} = 0.0;
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