function [Qsum, alpha, alpha_hist, cost_hist] = MVOE(Q, alpha0, rho)
% now alpha0 is act as pivoting point
reg = 1e-5;

N = size(Q,3);
d = size(Q,1);

maxIter = 5;

% initialization
if nargin == 1
    alpha = ones(1,N)/N;
    alpha0 = ones(1,N)/N;
    reg = 0;
elseif nargin == 2 % pivot activated
    alpha = alpha0;
else % nargin == 3
    % pivot activated with large regulation
    % this is usually bad case, so increase iteration
    alpha = alpha0;
    reg = rho;
    maxIter = 50;
end

if nargout > 3
    alpha_hist = zeros(maxIter+1,N);
    alpha_hist(1,:) = alpha;
    Q_ = alpha_sum(Q, alpha);
    cost_hist = zeros(maxIter+1,1);
    cost_hist(1) = log(det(Q_));
end

for iter = 1:maxIter
    Qsum = alpha_sum(Q,alpha);
    residue_all = sum(alpha.*(alpha - alpha0));
    
    alpha_new = zeros(size(alpha));
    success = true;
    for k = 1:N
        Qk = Q(:,:,k);
        residue = residue_all - alpha(k)*(alpha(k)-alpha0(k));
        alpha_squared_ = trace(Qsum\Qk) / (d - reg*residue);
        if alpha_squared_ < 0
            % stop updating alpha
            alpha_new(k) = sqrt(alpha_squared_);
        else
            alpha_new(k) = sqrt(alpha_squared_);
        end
    end
    
    if success 
        % accept current step
        alpha = alpha_new;
        alpha = alpha / sum(alpha);
    else
        % reject current step and stop
        disp('Stopped prematurely');
        break;
    end
    
    if nargout > 3        
        alpha_hist(iter+1,:) = alpha;
        Q_ = alpha_sum(Q, alpha);
        cost_hist(iter+1) = log(det(Q_)) + 0.5*reg*norm(alpha-alpha0)^2;
    end
end
end

function Qsum = alpha_sum(Q, alpha)
Qsum = zeros(size(Q(:,:,1)));
for k = 1:size(Q,3)
    Qsum = Qsum + Q(:,:,k) / alpha(k);
end
end
