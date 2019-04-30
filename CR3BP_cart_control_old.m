function Xdot = CR3BP_cart_control(t,X,mu,c)
    
    %% Unpack
    x = X(1);
    y = X(2);
    z = X(3);
    xdot = X(4);
    ydot = X(5);
    zdot = X(6);
    m = X(7);
    ux = X(8);
    uy = X(9);
    uz = X(10);
    
    STM = reshape(X(11:end),10,[]);
    
    %% Dynamics
    
    r1 = sqrt((x+mu)^2 + y^2 + z^2); 
    r2 = sqrt((x-1+mu)^2 + y^2 + z^2);

    xddot = 2*ydot + x -(1-mu)*((x+mu)/(r1^3)) - mu*(x-1+mu)/(r2^3) + (ux/m);%*200;%/0.002730736489750;
    yddot = -2*xdot + y - (1-mu)*(y/(r1^3)) - mu*(y)/(r2^3) + (uy/m);%*200;%/0.002730736489750;
    zddot = -(1-mu)*((z)/(r1^3)) - mu*(z)/(r2^3) + (uz/m);%*200;%/0.002730736489750;
    u_epsilon = [1e-8; 1e-10; 1e-12];
    mdot = -sqrt(ux^2 + uy^2 + uz^2 + u_epsilon(3))/c;
%     mdot = -(ux + uy + uz)/c;
    
    %% STM Dynamics
    fXmat = fX(X,c,mu);
    STMdot = fXmat*STM;
    
    %% Repack
    Xdot = [xdot; ydot; zdot; xddot; yddot; zddot; mdot; 0; 0; 0; reshape(STMdot,[],1)];
    
end