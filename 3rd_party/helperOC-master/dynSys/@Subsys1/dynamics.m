function dx = dynamics(obj, ~, x, u, d)
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

if nargin < 5
    d = [0; 0; 0];
end

if iscell(x)
    dx = cell(length(obj.dims), 1);
    
    for i = 1:length(obj.dims)
        dx{i} = dynamics_cell_helper(obj, x, u, d, obj.dims, obj.dims(i));
    end
else
    dx = zeros(obj.nx, 1);
    dx(1) = -obj.gains(1) * (x(1) - u(1)) + d(1);
    dx(2) = -obj.gains(2) * (x(2) - u(2)) + d(2);
    dx(3) = x(4);
    dx(4) = -obj.cD/obj.mass * x(4) + 9.8 * sin(x(1)) + x(2) / obj.mass * sin(x(1)) + d(3);
end
end

function dx = dynamics_cell_helper(obj, x, u, d, dims, dim)

switch dim
    case 1
        dx = -obj.gains(1) * (x{1} - u{1}) + d{1};
    case 2
        dx = -obj.gains(2) * (x{2} - u{2}) + d{2};
    case 3
        dx = x{4};
    case 4
        dx = -obj.cD/obj.mass * x{4} + 9.8 * sin(x{1}) + (x{2} / obj.mass).* sin(x{1}) + d{3};
    otherwise
        error('Only dimension 1-3 are defined for dynamics of Subsys1!')
end
end