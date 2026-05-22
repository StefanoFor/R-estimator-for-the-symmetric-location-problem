clear all
close all
clc

% Monte carlo runs
Ns=10^6;

% Number of observations
n = 100;

% Center of symmetry
theta_0 = 6;

% "Small perturabtion" parameter
pert = 5;

% Parameter of the t_nu score
nu_par = 4;

% Degree of freedom of the t-distribution
nuvect=[1:0.5:10];
Nl=length(nuvect);

for il=1:Nl

    nu = nuvect(il)

    % Fisher Information
    I_0(il) = (nu+1)/(nu+3);

    MSE_mean = 0;
    MSE_med = 0;
    MSE_R_G_Psi_Con = 0;
    MSE_R_G_Psi_Rob = 0;
    MSE_MLE_t_nu = 0;
    hat_I_G = 0;

    parfor ins=1:Ns

        % Generation of the t-distributed data
        y = t_nu_data(n,nu,theta_0)

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
        % Semiparametric One-Step R-estimator with the consistent and the robust estimator of Psi, Gaussian score
        [t_R_G_C, t_R_G_R, hat_psi_G]  = R_est_sym_loc_G(y, t_med, pert);

        err_R_G = t_R_G_C-theta_0;
        MSE_R_G_Psi_Con = MSE_R_G_Psi_Con + err_R_G^2/Ns;

        err_R_G_R = t_R_G_R-theta_0;
        MSE_R_G_Psi_Rob = MSE_R_G_Psi_Rob + err_R_G_R^2/Ns;

        hat_I_G = hat_I_G + hat_psi_G/Ns;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Oracle ML estimator for t-distributed data
        MLE_t = MLE_sym_loc_t_nu(y, t_med, nu);

        err_MLE_t_nu = MLE_t-theta_0;
        MSE_MLE_t_nu = MSE_MLE_t_nu + err_MLE_t_nu^2/Ns;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    Fro_MSE_mean(il) = MSE_mean;
    Fro_MSE_med(il) = MSE_med;
    Fro_MSE_R_G_Psi_Con(il) = MSE_R_G_Psi_Con;
    Fro_MSE_R_G_Psi_Rob(il) = MSE_R_G_Psi_Rob;
    Fro_MSE_MLE_t_nu(il) = MSE_MLE_t_nu;
    CRB(il) = 1/I_0(il)/n;
    hat_I_0_G(il) = hat_I_G;
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
semilogy(nuvect,Fro_MSE_mean,line_marker{3},'LineWidth',1,'Color',color_matrix(3,:),'MarkerEdgeColor',color_matrix(3,:),'MarkerFaceColor',color_matrix(3,:),'MarkerSize',8);
hold on
semilogy(nuvect,Fro_MSE_med,line_marker{4},'LineWidth',1,'Color',color_matrix(4,:),'MarkerEdgeColor',color_matrix(4,:),'MarkerFaceColor',color_matrix(4,:),'MarkerSize',8);
hold on
semilogy(nuvect,Fro_MSE_R_G_Psi_Con,line_marker{5},'LineWidth',1,'Color',color_matrix(5,:),'MarkerEdgeColor',color_matrix(5,:),'MarkerFaceColor',color_matrix(5,:),'MarkerSize',8);
hold on
semilogy(nuvect,Fro_MSE_R_G_Psi_Rob,line_marker{7},'LineWidth',1,'Color',color_matrix(7,:),'MarkerEdgeColor',color_matrix(7,:),'MarkerFaceColor',color_matrix(7,:),'MarkerSize',8);
hold on
semilogy(nuvect,Fro_MSE_MLE_t_nu,line_marker{6},'LineWidth',1,'Color',color_matrix(6,:),'MarkerEdgeColor',color_matrix(6,:),'MarkerFaceColor',color_matrix(6,:),'MarkerSize',8);
grid on
semilogy(nuvect,CRB,line_marker{1},'LineWidth',1,'Color',color_matrix(1,:),'MarkerEdgeColor',color_matrix(1,:),'MarkerFaceColor',color_matrix(1,:),'MarkerSize',8);
grid on
axis([1 10 0 2*max(Fro_MSE_R_G_Psi_Con)])
xlabel('Degrees of freedom: $\nu$','interpreter','latex');ylabel('MSE \& Lower Bound','interpreter','latex');
legend('Sample mean','Sample median','$R$-est with G-score and Consistent $\widehat{\Psi}$','$R$-est with G-score and Robust $\widehat{\Psi}$','Oracle MLE','Lower bound','interpreter','latex')
lgd = legend;
set(lgd,'FontSize',16);
set(gca, 'FontSize', 16);

figure(2)
semilogy(nuvect,I_0)
hold on
semilogy(nuvect,hat_I_0_G)
grid on
legend('$I_0$','$\hat{I}_0$ G score','interpreter','latex')
lgd = legend;
set(lgd,'FontSize',16);
set(gca, 'FontSize', 16);


%save dataPlot Fro_MSE_mean Fro_MSE_med Fro_MSE_R_G_Psi_Con Fro_MSE_R_G_Psi_Rob Fro_MSE_MLE_t_nu CRB nuvect