clear all
close all
clc

% Monte carlo runs
Ns=10^5;

% Number of observations
n = 100;

% Center of symmetry
theta_0 = 6;

% Quantile for the estimation of Psi
alpha_psi = 0.1; % Default value 0.1
% Smaller alpha: More robust to outliers (captures more extreme residuals).
% Bugger alpha: More efficient for clean data but less robust.

% Tuning parameter for the Huber M-estimator
alpha_H = 0.9; % Default value 1.345
% Smaller alpha: More robust to outliers (alpha_H -> 0,  Huber M-est = sample median).
% Bugger alpha: More efficient for clean data but less robust (alpha_H -> 0,  Huber M-est = sample mean).

% Maximum number of iteration for the Huber M-est
max_iter_H = 5;

% Shape parameter of the GG distribution
svect=[0.4:0.2:4];
Nl=length(svect);

for il=1:Nl

    s = svect(il)
    
    % Scale parameter of the GG distribution such that, for s = 1, g is a
    % Gaussian distribution with varaince equal to 1
    % b =( (1/2)*gamma(1/(2*s))/ gamma( (1+2)/(2*s) ) )^s;
    b= 1;

    % Fisher Information
    I_0(il) = 4*s^2 * gamma( 2 - 1/(2*s) ) / gamma(1/(2*s)) / (b*2^s)^(1/s);

    MSE_mean = 0;
    MSE_med = 0;
    MSE_H = 0;
    MSE_R_G_Psi_Rob = 0;
    MSE_MLE_t_nu = 0;
    hat_I_G_R = 0;

    parfor ins=1:Ns

        % Generation of the GG-distributed data
        y = GG_s_data(n,s,b,theta_0);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample mean
        t_mean = mean(y);

        err_mean = t_mean-theta_0;
        MSE_mean = MSE_mean + err_mean^2/Ns;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Sample median
        t_med = median(y);

        err_med = t_med-theta_0;
        MSE_med = MSE_med + err_med^2/Ns;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Huber M-estimator
        t_H = Huber_M_est_loc(y, alpha_H, max_iter_H);

        err_H = t_H-theta_0;
        MSE_H = MSE_H + err_H^2/Ns;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Semiparametric One-Step R-estimator with the consistent and the robust estimator of Psi, Gaussian score
        [t_R_G_R, hat_psi_G_R]  = R_est_sym_loc_G(y, t_med, alpha_psi);

        err_R_G_R = t_R_G_R-theta_0;
        MSE_R_G_Psi_Rob = MSE_R_G_Psi_Rob + err_R_G_R^2/Ns;

        hat_I_G_R = hat_I_G_R + hat_psi_G_R/Ns;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    Fro_MSE_mean(il) = MSE_mean;
    Fro_MSE_med(il) = MSE_med;
    Fro_MSE_H(il) = MSE_H;
    Fro_MSE_R_G_Psi_Rob(il) = MSE_R_G_Psi_Rob;
    CRB(il) = 1/I_0(il)/n;
    hat_I_0_G_R(il) = hat_I_G_R;
end


color_matrix(1,:)=[0 0 1]; % Blue
color_matrix(2,:)=[1 0 0]; % Red
color_matrix(3,:)=[0 0.5 0]; % Dark Green
color_matrix(4,:)=[0 0 0]; % Black
color_matrix(5,:)=[0 0.5 1]; % Light Blue
color_matrix(6,:)=[1 0.3 0.6]; % Pink
color_matrix(7,:)=[0 0.9 0]; % Light Green

line_marker{1}='-s';
line_marker{2}='--d';
line_marker{3}=':^';
line_marker{4}='-.p';
line_marker{5}='-o';
line_marker{6}='--h';
line_marker{7}='-.*';

figure(1)
semilogy(svect,Fro_MSE_mean,line_marker{3},'LineWidth',1,'Color',color_matrix(3,:),'MarkerEdgeColor',color_matrix(3,:),'MarkerFaceColor',color_matrix(3,:),'MarkerSize',8);
hold on
semilogy(svect,Fro_MSE_med,line_marker{4},'LineWidth',1,'Color',color_matrix(4,:),'MarkerEdgeColor',color_matrix(4,:),'MarkerFaceColor',color_matrix(4,:),'MarkerSize',8);
hold on
semilogy(svect,Fro_MSE_H,line_marker{6},'LineWidth',1,'Color',color_matrix(6,:),'MarkerEdgeColor',color_matrix(6,:),'MarkerFaceColor',color_matrix(6,:),'MarkerSize',8);
hold on
semilogy(svect,Fro_MSE_R_G_Psi_Rob,line_marker{7},'LineWidth',1,'Color',color_matrix(7,:),'MarkerEdgeColor',color_matrix(7,:),'MarkerFaceColor',color_matrix(7,:),'MarkerSize',8);
grid on
semilogy(svect,CRB,line_marker{1},'LineWidth',1,'Color',color_matrix(1,:),'MarkerEdgeColor',color_matrix(1,:),'MarkerFaceColor',color_matrix(1,:),'MarkerSize',8);
grid on
axis([min(svect) max(svect) 0 2*max(Fro_MSE_med)])
xlabel('Shape parameter: $s$','interpreter','latex');ylabel('MSE \& Lower Bound','interpreter','latex');
legend('Sample mean','Sample median','Huber $M$-est','$R$-est with Gausian score$','Lower bound','interpreter','latex')
lgd = legend;
set(lgd,'FontSize',16);
set(gca, 'FontSize', 16);

figure(2)
semilogy(svect,I_0)
hold on
semilogy(svect,hat_I_0_G_R)
grid on
legend('$I_0$','$\hat{I}_{R}$','interpreter','latex')
lgd = legend;
set(lgd,'FontSize',16);
set(gca, 'FontSize', 16);