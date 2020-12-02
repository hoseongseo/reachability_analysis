function dOpt = optDstb(obj, t, y, deriv, dMode)

if ~iscell(deriv)
    deriv = num2cell(deriv);
end

dOpt = cell(obj.nd, 1);

%% Optimal Disturbance
if strcmp(dMode, 'max')
    for i = 1:3
        dOpt{i} = (deriv{obj.dims==i}>=0)*obj.dMax(i) + (deriv{obj.dims==i}<0)*(-obj.dMax(i));
    end
elseif strcmp(dMode, 'min')
    for i = 1:3
        dOpt{i} = (deriv{obj.dims==i}>=0)*(-obj.dMax(i)) + (deriv{obj.dims==i}<0)*obj.dMax(i);
    end
else
    error('Unknown dMode!')
end

end