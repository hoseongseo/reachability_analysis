function dOpt = optDstb(obj, t, y, deriv, dMode)
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

dOpt = cell(obj.nd, 1);

%% Optimal Disturbance
deriv_dim = [2,3];

if strcmp(dMode, 'max')
    for i = 1:2
        dOpt{i} = (deriv{obj.dims==deriv_dim(i)}>=0)*obj.dMax(i) +...
                  (deriv{obj.dims==deriv_dim(i)}<0)*(-obj.dMax(i));
    end
elseif strcmp(dMode, 'min')
    for i = 1:2
        dOpt{i} = (deriv{obj.dims==deriv_dim(i)}>=0)*(-obj.dMax(i)) +...
                  (deriv{obj.dims==deriv_dim(i)}<0)*obj.dMax(i);
    end
else
    error('Unknown dMode!')
end

end