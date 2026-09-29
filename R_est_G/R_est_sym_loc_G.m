function [theta_R_R, hat_psi_R] = R_est_sym_loc_G(x, theta_star, alpha_psi)

n = length(x);

[Delta_G, ~] = eff_cent_seq_G(x, theta_star);

% Estimator of Psi
hat_psi_R = 1/est_psi_hat_vdw(x, theta_star, alpha_psi);

theta_R_R = theta_star + (1/sqrt(n)) * Delta_G / hat_psi_R;

end