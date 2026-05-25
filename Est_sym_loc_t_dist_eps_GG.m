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
pert = 2;

% Degree of freedom of the t-distribution
nu = 4;

% Shape and scale parameters of the GG contaminating distribution
s = 0.9; % if s = 1, the GG density is actually the Gaussian density 
%b =( (1/2)*gamma(1/(2*s))/ gamma( (1+2)/(2*s) ) )^s; %for s = 1, b is the variance of the resulting Gaussian distribution
b = 10; 

% Contamination parameter in (0,0.5)
eps_vect=[0:0.02:0.3];
Nl=length(eps_vect);

for il=1:Nl

    eps_cond = eps_vect(il)

    MSE_mean = 0;
    MSE_med = 0;
    MSE_R_G_Psi_Con = 0;
    MSE_R_G_Psi_Rob = 0;
    MSE_MLE_t_nu = 0;
    hat_I_G = 0;

    parfor ins=1:Ns

        % Generation of the t-distributed GG_contaminated data
        r_cond = rand(1,n)>=eps_cond;
        GG_dist = GG_s_data(n,s,b,theta_0);
        t_dist = t_nu_data(n,nu,theta_0);
        y = r_cond.*t_dist + (1-r_cond).*GG_dist;

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

    end
    Fro_MSE_mean(il) = MSE_mean;
    Fro_MSE_med(il) = MSE_med;
    Fro_MSE_R_G_Psi_Con(il) = MSE_R_G_Psi_Con;
    Fro_MSE_R_G_Psi_Rob(il) = MSE_R_G_Psi_Rob;
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
semilogy(eps_vect,Fro_MSE_mean,line_marker{3},'LineWidth',1,'Color',color_matrix(3,:),'MarkerEdgeColor',color_matrix(3,:),'MarkerFaceColor',color_matrix(3,:),'MarkerSize',8);
hold on
semilogy(eps_vect,Fro_MSE_med,line_marker{4},'LineWidth',1,'Color',color_matrix(4,:),'MarkerEdgeColor',color_matrix(4,:),'MarkerFaceColor',color_matrix(4,:),'MarkerSize',8);
hold on
semilogy(eps_vect,Fro_MSE_R_G_Psi_Con,line_marker{5},'LineWidth',1,'Color',color_matrix(5,:),'MarkerEdgeColor',color_matrix(5,:),'MarkerFaceColor',color_matrix(5,:),'MarkerSize',8);
hold on
semilogy(eps_vect,Fro_MSE_R_G_Psi_Rob,line_marker{7},'LineWidth',1,'Color',color_matrix(7,:),'MarkerEdgeColor',color_matrix(7,:),'MarkerFaceColor',color_matrix(7,:),'MarkerSize',8);
grid on
axis([min(eps_vect) max(eps_vect) 0 2*max(Fro_MSE_R_G_Psi_Con)])
xlabel('contamination parameter: $\epsilon$','interpreter','latex');ylabel('MSE','interpreter','latex');
legend('Sample mean','Sample median','$R$-est with G-score and Consistent $\widehat{\Psi}$','$R$-est with G-score and Robust $\widehat{\Psi}$','interpreter','latex')
lgd = legend;
set(lgd,'FontSize',16);
set(gca, 'FontSize', 16);

save dataPlot Fro_MSE_mean Fro_MSE_med Fro_MSE_R_G_Psi_Con Fro_MSE_R_G_Psi_Rob eps_vect