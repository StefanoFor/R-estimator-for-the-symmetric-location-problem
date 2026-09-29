function mu_hat = Huber_M_est_loc(x, c, max_iter)

    x = x(:);  % make x a column vector

    % Initial estimate
    mu_hat = median(x);

    % Iterative reweighted least squares
    tol = 1e-6;

    for k = 1:max_iter

        r = x - mu_hat;

        % Huber weights
        w = ones(size(r));
        idx = abs(r) > c;
        w(idx) = c ./ abs(r(idx));

        % Update
        mu_new = sum(w .* x) / sum(w);

        % Check convergence
        if abs(mu_new - mu_hat) < tol
            mu_hat = mu_new;
            break;
        end

        mu_hat = mu_new;
    end
end