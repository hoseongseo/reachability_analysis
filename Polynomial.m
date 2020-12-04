classdef Polynomial < handle
    properties (GetAccess = public, SetAccess = private)
        coeff;
        order;
    end
    methods
        function obj = Polynomial(cff, ord)
            if nargin == 2
                obj.coeff = cff;
                obj.order = ord;
            end
        end
        function p = diff(obj, idx)
            if nargin > 1
            e_ = zeros(1,size(obj.order,2));
            e_(idx) = 1;
            
            order_ = obj.order - repmat(e_, [size(obj.order,1),1]);
            idx_ = order_(:,idx) >= 0;
            order_ = order_(idx_,:);
            
            A = zeros(sum(idx_), length(idx_));
            A(:,idx_) = eye(sum(idx_));
%             if class(obj.coeff) == "function_handle"
%                 coeff_ = @(param) A * (obj.coeff(param) .* obj.order(:,idx));
%             else
                coeff_ = A * (obj.coeff .* obj.order(:,idx));
%             end
            p = Polynomial(coeff_, order_);
            elseif nargin == 1
                p(1:size(obj.order,2)) = Polynomial;
                for j = 1:size(obj.order,2)
                    p(j) = diff(obj, j);
                end
            end
        end
        function p = collect(obj, idx_sub)
            idx_rem = setdiff( 1:size(obj.order,2), idx_sub );
            
            order_sub = obj.order(:,idx_sub);
            max_order_sub = max(order_sub, [], 1);
            
            id_sub = Polynomial.getId( order_sub, max_order_sub );
            id_sub_unique = unique( id_sub );
            
            cff_collected( 1:length(id_sub_unique), 1) = Polynomial;
            for i = 1:length(id_sub_unique)
                tmp = Polynomial;
                tmp.order = obj.order( id_sub == id_sub_unique(i), idx_rem );
                
                idx_ = (id_sub == id_sub_unique(i));
                A = zeros(sum(idx_), length(idx_));
                A(:,idx_) = eye(sum(idx_));
                if class(obj.coeff) == "function_handle"    
                    tmp.coeff = @(param) A * obj.coeff(param);
                else
                    tmp.coeff = A * obj.coeff;
                end
                cff_collected(i) = tmp;
            end
            
            p = Polynomial(cff_collected, Polynomial.getOrder(id_sub_unique, max_order_sub));
        end
        function V = eval(obj, gr)
            V = zeros(gr.shape);
            for i = 1:length(obj.coeff)
                tmp = obj.coeff(i)*ones(gr.shape);
                for j = 1:size(obj.order,2)
                    tmp = tmp.*(gr.xs{j}.^obj.order(i,j));
                end
                V = V + tmp;
            end
        end
        
        function setCoeff(obj, cff)
            obj.coeff = cff;
        end
        
        %%%%% Operator overloadings
        function p = plus1(p1,p2)
            if class(p2) == "Polynomial"
                max_order_ = max( max(p1.order,[],1), max(p2.order,[],1) );
                id1_ = Polynomial.getId( p1.order, max_order_ );
                id2_ = Polynomial.getId( p2.order, max_order_ );
                id_ = union( id1_, id2_ );
                
                [~, i1_] = intersect(id_, id1_);
                [~, i2_] = intersect(id_, id2_);
                
                A1 = zeros(length(id_), length(id1_));
                A2 = zeros(length(id_), length(id2_));
                A1(i1_,:) = eye(length(id1_));
                A2(i2_,:) = eye(length(id2_));
                
                coeff_ = A1 * p1.coeff + A2 * p2.coeff;
                order_ = Polynomial.getOrder( id_, max_order_ );
            else
                max_order_ = max(p1.order,[],1);
                id1_ = Polynomial.getId( p1.order, max_order_ );
                if id1_(1) > 1
                    coeff_ = [p2; p1.coeff];
                    order_ = [zeros(1,size(p1.order,2)); p1.order];
                else % id1_(1) == 1
                    coeff_ = p1.coeff;
                    coeff_(1) = coeff_(1) + p2;
                    order_ = p1.order;
                end
            end
            p = Polynomial(coeff_, order_);
        end
        function p = plus(p1,p2)
            p = arrayOperation(p1,p2,@plus1);
        end
        function p = sum(p1)
            p = Polynomial(p1(1).coeff, p1(1).order);
            for i = 2:numel(p1)
                p = p + p1(i);
            end
        end
        function p = uminus(p1)
            sz = size(p1);
            p(1:sz(1),1:sz(2)) = Polynomial;
            for i = 1:numel(p1)
                p(i) = Polynomial(-p1(i).coeff, p1(i).order);
            end
        end
        function p = minus(p1,p2)
            p = p1 + (-p2);
        end
        function p = times1(p1,p2)
            if class(p2) == "Polynomial"
                max_order_ = max(p1.order,[],1) + max(p2.order,[],1);
                id1_ = Polynomial.getId( p1.order, max_order_ );
                id2_ = Polynomial.getId( p2.order, max_order_ );
                
                id_all_ = Utils.minksum( id1_, id2_ ) - 1;
                id_ = unique( id_all_ );
                
%                 B = repmat(id_all_', [length(id_),1]) == repmat(id_, [1,length(id_all_)]);
                B = zeros(length(id_), length(id_all_));
                for i = 1:length(id_)
                    B(i,:) = (id_all_ == id_(i));
                end
                
                coeff_ = B * kron(p1.coeff, p2.coeff);
                order_ = Polynomial.getOrder( id_, max_order_ );
            else
                coeff_ = p2*p1.coeff;
                order_ = p1.order;
            end
            p = Polynomial( coeff_, order_ );
        end
        function p = times(p1,p2)
            p = arrayOperation(p1,p2,@times1);
        end
        function p = mtimes(p1,p2)
            s1 = size(p1);
            s2 = size(p2);
            if numel(p1) == 1 || numel(p2) == 1
                p = p1.*p2;
            elseif s1(2) == s2(1)
                p(1:s1(1),1:s2(2)) = Polynomial;
                idx = Polynomial.integer_grid([1,1], [s1(1), s2(2)]);
                for k = 1:size(idx,1)
                    i_ = idx(k,1);
                    j_ = idx(k,2);
                    p(i_, j_) = sum(p1(i_,:).*p2(:,j_)');
                end
            end
        end
        function p = power(p1, n)
            if n > 2
                p = power(p1,n-1).*p1;
            elseif n == 2
                tmp = p1.*p1;
                p = Polynomial(tmp.coeff, tmp.order);
            elseif n == 1
                p = Polynomial(p1.coeff, p1.order);
            elseif n == 0
                p = Polynomial(1, zeros(1,size(p1.order,2)));
            end
        end
        function compress(obj)
            eff_idx_ = obj.coeff ~= 0;
%             obj.max_order = max( obj.order(eff_idx_,:), [], 1 );
            obj.order = obj.order( eff_idx_, : );
            obj.coeff = obj.coeff( eff_idx_ );
        end
    end
    methods (Static)
        function id = getId(order, max_order)
            id = order(:,1)+1;
            for i = 2:size(order,2)
                id = id + order(:,i)*(prod(max_order(1:(i-1))+1));
            end
        end
        function order = getOrder(id, max_order)
            p = max_order+1;
            order = zeros(length(id), length(max_order));
            if size(id,2) > 1
                id = id(:);
            end
            val = id-1;
            for i = 1:length(max_order)-1
                order(:,i) = rem(val, p(i));
                val = (val - order(:,i))/p(i);
            end
            order(:,end) = val;
        end
        function p = fromSym(s, x)
            if nargin < 2
                p = Polynomial.fromSym(s, s);
            else
                p(1:size(s,1), 1:size(s,2)) = Polynomial;
                for k = 1:numel(s)
%                     if class(s(k)) == "sym"
                        [c_poly, x_poly] = coeffs(s(k), x);
                        if ~isempty(c_poly)
                            order = zeros(length(x_poly), length(x));
                            for i = 1:length(x)
                                order(:,i) = polynomialDegree(x_poly, x(i));
                            end
                            
                            id = Polynomial.getId( order, max(order, [], 1) );
                            [~, sorted_] = sort(id);
                            order_sorted = order(sorted_, :);
                            coeff_sorted = c_poly( sorted_ )';
                            if isempty( symvar(coeff_sorted) )
                                coeff_sorted = double(coeff_sorted);
                            end
                            p(k) = Polynomial( coeff_sorted, order_sorted );
                        else
                            p(k) = Polynomial( double(s(k)), zeros(1,length(x)) );
                        end
%                     else
%                         p(k) = Polynomial( s(k), zeros(1,length(x)) );
%                     end
                end
            end
        end
    end
    methods (Hidden)
        function [p1,p2] = swap(p1,p2)
            tmp = p1;
            p1 = p2;
            p2 = tmp;
        end
        function p = arrayOperation(p1,p2,fun)
            if class(p1) == "double"
                [p1,p2] = swap(p1,p2);
            end
            if numel(p1) == 1 && numel(p2) == 1
                p = fun(p1,p2);
            elseif numel(p1) > 1 && numel(p2) == 1
                sz = size(p1);
                p(1:sz(1),1:sz(2)) = Polynomial;
                for i = 1:numel(p1)
                    p(i) = fun(p1(i),p2);
                end
            elseif numel(p1) == 1 && numel(p2) > 1
                sz = size(p2);
                p(1:sz(1),1:sz(2)) = Polynomial;
                for i = 1:numel(p2)
                    p(i) = fun(p1,p2(i));
                end
            elseif prod(size(p1) == size(p2))
                sz = size(p1);
                p(1:sz(1),1:sz(2)) = Polynomial;
                for i = 1:numel(p1)
                    p(i) = fun(p1(i),p2(i));
                end
            end
        end
    end
end