function y = t_nu_data(n,nu,theta_0)

m = 1;
% Generation of the t-distributed data
w = randn(m,n);
R = gamrnd(nu/2,2/nu,1,n);
y = theta_0 + sqrt(1./(repmat(R,m,1))).*w;

end