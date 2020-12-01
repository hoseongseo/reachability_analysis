classdef Sphere < handle
    %%%% Points on (n+1) dimensional unit sphere
    %%%% S(n) = {x \in R(n+1) | norm(x) = 1}
    properties (GetAccess = public, SetAccess = public)
        x
        phi
    end
    methods (Static)
        function this = Sphere(n, numPts_est)
            % n : dimension -> S(n)
            % numPts : number of total points
            
            % estimated number of points in each grid
            N = round(exp( (1/n) * (log(numPts_est) - log(2)) )); 
            sz = N * ones(1,n);
            sz(end) = 2*N-2; % since periodic
            numPts = prod(sz);
            this.phi = zeros(n, numPts);
            
            gr = Utils.Grid( sz );
            dpi = (gr.x{1}(2,1) - gr.x{1}(1,1))*pi;
            
            for i = 1:n
                if i < n
                    this.phi(i,:) = gr.x{i}(:) * pi;
                else
                    this.phi(i,:) = gr.x{i}(:) * (2*pi - dpi);
                end
            end
            
            this.x = zeros(n+1, numPts);
            this.x(1,:) = cos(this.phi(1,:));
            for i = 2:n
                this.x(i,:) = this.x(i-1,:) .* tan(this.phi(i-1,:)) .* cos(this.phi(i,:));
            end
            this.x(end,:) = this.x(n,:) .* tan(this.phi(n,:));
            this.x = [this.x, this.x(:,1)];
        end
        function this = rand(n, numPts_est)
            % n : dimension -> S(n)
            % numPts : number of total points
            
            % estimated number of points in each grid
            N = round(exp( (1/n) * (log(numPts_est) - log(2)) )); 
            sz = N * ones(1,n);
            sz(end) = 2*N-2;
            
            numPts = prod(sz);
            this.phi = zeros(n, numPts);
            
            gr = Utils.Grid.rand( sz );
            for i = 1:n
                if i < n
                    this.phi(i,:) = gr.x{i}(:) * pi;
                else
                    this.phi(i,:) = gr.x{i}(:) * 2*pi;
                end
            end
            
            this.x = zeros(n+1, numPts);
            this.x(1,:) = cos(this.phi(1,:));
            for i = 2:n
                this.x(i,:) = this.x(i-1,:) .* tan(this.phi(i-1,:)) .* cos(this.phi(i,:));
            end
            this.x(end,:) = this.x(n,:) .* tan(this.phi(n,:));
            this.x = [this.x, this.x(:,1)];
        end
    end
end