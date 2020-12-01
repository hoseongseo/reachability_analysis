function uOpt = optCtrl(obj, t, y, ~, ~)
% uOpt = optCtrl(obj, t, y, deriv, uMode)

xg = y{1};
yg = y{2};

px = obj.px;
py = obj.py;

b0 = PolyBasis(length(px)-1);
b1 = b0.diff(1);

B0 = blkdiag(b0.eval(t), b0.eval(t));
B1 = blkdiag(b1.eval(t), b1.eval(t));

Xcur = B0 * [px; py];
dXcur = B1 * [px; py];

%% Optimal control
uOpt{1} = -obj.k*(xg - Xcur(1)) + dXcur(1);
uOpt{2} = -obj.k*(yg - Xcur(2)) + dXcur(2);

end