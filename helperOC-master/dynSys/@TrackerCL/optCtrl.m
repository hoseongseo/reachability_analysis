function uOpt = optCtrl(obj, t, y, deriv, uMode)

b0 = PolyBasis(length(obj.p)-1);
b1 = b0.diff(1);
b2 = b1.diff(1);

% uOpt = cell(obj.nu, 1);

%% Optimal control
if strcmp(obj.tMode, 'forward')
    x_ref = b0.eval(t) * obj.p;
    dx_ref = b1.eval(t) * obj.p;
    d2x_ref = b2.eval(t) * obj.p;
    
    e = y{1};
    de = y{2};
    uOpt{1} = -obj.Cd*(de + dx_ref) - obj.k(1) * e - obj.k(2) * de + (obj.k(3)-1) * d2x_ref;
%     uOpt{1} = -obj.k(1) * (y{1} - x_ref) - obj.k(2) * (y{2} - dx_ref) + obj.k(3) * d2x_ref;
elseif strcmp(obj.tMode,'backward')
    x_ref = b0.eval(obj.tMax - t) * obj.p;
    dx_ref = b1.eval(obj.tMax - t) * obj.p;
    d2x_ref = b2.eval(obj.tMax - t) * obj.p;
    
    e = y{1};
    de = y{2};
    uOpt{1} = -obj.Cd*(de + dx_ref) - obj.k(1) * e - obj.k(2) * de + (obj.k(3)-1) * d2x_ref;
%     uOpt{1} = -obj.k(1) * (y{1} - x_ref) - obj.k(2) * (y{2} - dx_ref) + obj.k(3) * d2x_ref;    
end