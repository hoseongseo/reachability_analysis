classdef hss < DynSys
  properties
    dMax % Disturbance
    dims
    
    k;
    px;
    py;
  end
  
  methods
    function obj = hss(x, dMax, dims, k, px, py)
      
      if numel(x) ~= obj.nx
          error('Initial state does not have right dimension!');
      end
      
      if ~iscolumn(x)
          x = x';
      end
      
      % Basic vehicle properties
      obj.pdim = [find(dims == 1) find(dims == 2)]; % Position dimensions
      obj.nx = length(dims);
      obj.nu = 2;
      obj.nd = 2;
      
      obj.x = x;
      obj.xhist = obj.x;
      
      obj.dMax = dMax;
      obj.dims = dims;
      obj.k = k;
      obj.px = px;
      obj.py = py;
    end
    
  end % end methods
end % end classdef
