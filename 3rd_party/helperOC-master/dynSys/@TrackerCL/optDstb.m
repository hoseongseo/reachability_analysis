function dOpt = optDstb(obj, ~, ~, deriv, dMode)

%% Input processing

if ~iscell(deriv)
    deriv = num2cell(deriv);
end

dOpt = cell(obj.nd, 1);

%% Optimal control
if strcmp(dMode, 'max')
    dOpt{1} = (deriv{obj.dims==2}>=0)*obj.dMax + (deriv{obj.dims==2}<0)*(-obj.dMax);
elseif strcmp(dMode, 'min')
    dOpt{1} = (deriv{obj.dims==2}>=0)*(-obj.dMax) + (deriv{obj.dims==2}<0)*obj.dMax;
else
    error('Unknown dMode!')
end

end