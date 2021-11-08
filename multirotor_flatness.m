function [R, w, dw] = multirotor_flatness(a, j, s)

g = [0; 0; 9.8];
e3 = [0; 0; 1];

N = size(a, 2);

R = zeros(3,3,N);
w = zeros(3,N);
dw = zeros(3,N);
for i = 1:N
    a_ = a(:,i);
    j_ = j(:,i);
    s_ = s(:,i);
    
    psi_ = 0.0;
    dpsi_ = 0.0;
    d2psi_ = 0.0;
    
    f_ = norm(a_ + g);

    z3_ = (a_ + g) / norm(a_ + g);
    zc_ = [cos(psi_); sin(psi_); 0.0];
    z2_ = cross( z3_, zc_ );
    z2_ = z2_ / norm(z2_);
    z1_ = cross( z2_, z3_ );
    z1_ = z1_ / norm(z1_);
    R_ = [z1_, z2_, z3_];

    wx_ = -(1/f_)*(z2_'*j_);
    wy_ = (1/f_)*(z1_'*j_);
    wz_ = (dpsi_ - e3'*(z1_*wx_+z2_*wy_))/(e3'*z3_);
    w_ = [wx_; wy_; wz_];
    
%     E(:,i) = E_;
%     B_ = [1.0, sin(E_(1))*tan(E_(2)), cos(E_(1))*tan(E_(2));...
%         0.0, cos(E_(1)), -sin(E_(1));...
%         0.0, sin(E_(1))*sec(E_(2)), cos(E_(1))*sec(E_(2))];
%     dE_ = B_ * w_; % euler angle derivative
    
    dwx_ = (1/f_)*(-z2_'*s_ - 2*wx_*(z3_'*j_)) - wy_*wz_;
    dwy_ = (1/f_)*(z1_'*s_ - 2*wy_*(z3_'*j_)) + wx_*wz_;
    dwz_ = (d2psi_ - e3'*(z1_*dwx_ + z2_*dwy_))/(e3'*z3_);
    dw_ = [dwx_; dwy_; dwz_];

    R(:,:,i) = R_;
    w(:,i) = w_;    
    dw(:,i) = dw_;
end