function dx = dynamics(obj, ~, x, u, d)
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
    dx(1) = x(2);
    dx(2) = -obj.cD/obj.mass * x(2) + 9.8 * sin(x(3)) + d(1);
    dx(3) = -obj.gains(1) * (x(3) - u(1)) + d(2);
%     dx(3) = u(1) + d(2);
end
end

function dx = dynamics_cell_helper(obj, x, u, d, dims, dim)

switch dim
    case 1
        dx = x{2};
    case 2
        dx = -obj.cD/obj.mass * x{2} + 9.8 * sin(x{3}) + d{1};
    case 3
        dx = -obj.gains(1) * (x{3} - u{1}) + d{2};
%         dx = u{1} + d{2};
    otherwise
        error('Only dimension 1-3 are defined for dynamics of HssQuad1!')
end
end