classdef TrackerCL < DynSys
  properties
    dMax % Disturbance
    dims

    k;
    p;
    Cd;
    tMode;
    tMax
  end
  
  methods
    function obj = TrackerCL(x, dMax, dims, k, p, Cd, tMode, tMax)
      
      if numel(x) ~= obj.nx
          error('Initial state does not have right dimension!');
      end
      
      if ~iscolumn(x)
          x = x';
      end
      
      % Basic vehicle properties
      obj.nx = length(dims);
      obj.nu = 1;
      obj.nd = 1;
      
      obj.x = x;
      obj.xhist = obj.x;
      
      obj.dMax = dMax;
      obj.dims = dims;
      obj.k = k;
      obj.p = p;
      obj.Cd = Cd;
      obj.tMode = tMode;
      obj.tMax = tMax;
    end
    
  end % end methods
end % end classdef
