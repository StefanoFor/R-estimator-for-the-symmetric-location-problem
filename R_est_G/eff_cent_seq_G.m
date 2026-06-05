function [Delta_G, psi_G_R] = eff_cent_seq_G(x, theta_star)

n = length(x);

% Residuals
R = x - theta_star;

[~, sort_idx] = sort(abs(R));  % indices that would sort |R|
rank = zeros(1,n);
rank(sort_idx) = 1:n;           % ranks of |R|
U = rank / (n+1);               % empirical CDF of |R|

% Gaussian-based score
q = norminv( (1 + U)/2 ); 

Delta_G = sum(sign(R) .* q)/sqrt(n);
psi_G_R = sum(q.^2)/n;
end