function c_n = L1_L2_halo_c_n(mu, varargin)

    % Computes c_n coefficients for halo orbit around L1, L2 Lagrange
    % points
    
    if nargin == 1
        n = 4;
    elseif nargin == 2
        n = varargin;
    end

    c_n_EM = NaN(2,n); % First row is L1 coefficients, second is L2
    c_n_SE = NaN(2,n); % First row is L1 coefficients, second is L2

    for i = 1:2
        for j = 1:n
            if i==1
                c_n_EM(i,j) = 1/gammas_EM(i)^3*((1)^j*mu_EM + (-1)^j*((1-mu_EM)*gammas_EM(i)^(j+1))/((1 - gammas_EM(1))^(j+1)));
                c_n_SE(i,j) = 1/gammas_SE(i)^3*((1)^j*mu_SE + (-1)^j*((1-mu_SE)*gammas_SE(i)^(j+1))/((1 - gammas_SE(1))^(j+1)));
            elseif i==2
                c_n_EM(i,j) = 1/gammas_EM(i)^3*((-1)^j*mu_EM + (-1)^j*((1-mu_EM)*gammas_EM(i)^(j+1))/((1 + gammas_EM(1))^(j+1)));
                c_n_SE(i,j) = 1/gammas_SE(i)^3*((-1)^j*mu_SE + (-1)^j*((1-mu_SE)*gammas_SE(i)^(j+1))/((1 + gammas_SE(1))^(j+1)));
            end
        end
    end