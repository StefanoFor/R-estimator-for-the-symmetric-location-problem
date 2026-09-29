function hat_psi_H = est_psi_hat_vdw(x, theta_star, alpha)
    % Compute the estimator hat_psi_H using van der Waerden scores (phi(u) = norminv(u))
    % as per equation (2.7) in McKean & Hettmansperger (1978).
    %
    % Inputs:
    %   x: vector of observations
    %   theta_star: scalar or vector of model parameters (theta_star)
    %   alpha: significance level for z_alpha (default: 0.10)
    %
    % Output:
    %   hat_psi_H: estimated scale parameter (tau_hat)

    if nargin < 3
        alpha = 0.10; % Default alpha
    end

    residuals = x - theta_star; % Residuals from the fitted model
    n = length(residuals);
    z_alpha = norminv(1 - alpha); % Upper alpha quantile of standard normal

    % Van der Waerden score function: phi(u) = norminv(u)
    % One-sample scores: a^+(i) = phi^+ (i/(n+1)) = norminv( (i/(n+1) + 1)/2 )
    a_plus = @(i, n) norminv((i / (n + 1) + 1) / 2);

    % One-sample residual process S^+(t)
    S_plus = @(t) sum(a_plus(rank_abs(residuals, t), n) .* sign(residuals - t));

    % Helper function to compute ranks of |residuals - t|
    function ranks = rank_abs(e, t)
        abs_diff = abs(e - t);
        [~, idx] = sort(abs_diff);
        [~, ranks] = sort(idx);
    end

    % Target value for S^+(t)
    target = sqrt(n) * z_alpha;

    % Solve S^+(L) = target and S^+(U) = -target using fzero
    options = optimset('Display', 'off');
    L = fzero(@(t) S_plus(t) - target, [min(residuals), max(residuals)], options);
    U = fzero(@(t) S_plus(t) + target, [min(residuals), max(residuals)], options);

    % Compute tau_hat
    hat_psi_H = (sqrt(n) * (U - L)) / (2 * z_alpha);
end