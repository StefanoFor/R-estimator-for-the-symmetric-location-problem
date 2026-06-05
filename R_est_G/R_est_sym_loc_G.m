function [theta_R_C, theta_R_R, hat_psi_R] = R_est_sym_loc_G(x, theta_star, pert, alpha_psi)

[Delta_G, psi_G_R] = eff_cent_seq_G(x, theta_star);

h = pert*randn(1);
n = length(x);

[Delta_G_h, ~]= eff_cent_seq_G(x, theta_star + h/sqrt(n));

hat_psi_G = abs(Delta_G_h - Delta_G)/abs(h); 

if hat_psi_G < 0.0001
    hat_psi_G = psi_G_R;
end

% Consistent estimator of Psi
theta_R_C = theta_star + (1/sqrt(n)) * Delta_G / hat_psi_G;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Robust estimator of Psi
hat_psi_R = 1/est_psi_hat_vdw(x, theta_star, alpha_psi);
theta_R_R = theta_star + (1/sqrt(n)) * Delta_G / hat_psi_R;

end