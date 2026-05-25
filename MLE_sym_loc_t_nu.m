function theta_R = MLE_sym_loc_t_nu(x, theta_star, nu)

n = length(x);

% Residuals
R = x - theta_star;

[~, sort_idx] = sort(abs(R));  % indices that would sort |R|
rank = zeros(1,n);
rank(sort_idx) = 1:n;           % ranks of |R|
U = rank / (n+1);               % empirical CDF of |R|

% t_nu score
h = tinv( (1 + U)/2 , nu );

K = (nu+1)*h./(nu + h.^2); 

Delta_t = sum(sign(R) .* K);
psi_t = sum(K.^2);

theta_R = theta_star + Delta_t / psi_t;
end