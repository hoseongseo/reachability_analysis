function dx = dynamics(obj, ~, x, u, d)

if nargin < 5
    d = [0; 0];
end

if iscell(x)
    dx = cell(length(obj.dims), 1);
    
    for i = 1:length(obj.dims)
        dx{i} = dynamics_cell_helper(obj, x, u, d, obj.dims, obj.dims(i));
    end
else
    dx = zeros(obj.nx, 1);
    
    dx(1) = x(2);
%     dx(2) = - obj.k(1) * x(1) - (obj.k(2) + obj.Cd) * x(2) - obj.Cd * u(1) + d(1);
%     dx(2) = -obj.Cd * x(2) + u(1) + d(1);
    dx(2) = u(1) + d(1);
end
end

function dx = dynamics_cell_helper(obj, x, u, d, dims, dim)

switch dim
    case 1
        dx = x{2};
    case 2
%         dx = -obj.k(1) * x{1} - (obj.k(2) + obj.Cd) * x{2} - obj.Cd * u{1} + d{1};
%         dx = -obj.Cd * x{2} + u{1} + d{1};
        dx = u{1} + d{1};
    otherwise
        error('Only dimension 1-2 are defined for dynamics of hss!')
end
end