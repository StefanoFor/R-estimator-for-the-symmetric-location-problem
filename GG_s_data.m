function y = GG_s_data(n,s,b,theta_0)
m = 1;

% Generation of the GG-distributed data
w = randn(m,n);
w_norm = abs(w);
w_n = w./w_norm;
Q = gamrnd(1/(2*s),b*2^s,1,n);
y = theta_0 + Q.^(1/(2*s)).*w_n;
end