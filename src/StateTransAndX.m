function [t,PHItf,x,xf,PHI] = StateTransAndX(X0,tend,MU)
    % Variational Equation initial values (phi_0 and x_0)
    N = length(X0);
    Phi_t0_t0 = reshape(eye(N),N^2,1);
    VarEq_InitVal = zeros(N^2+N,1);
    VarEq_InitVal(1:N^2) = Phi_t0_t0;
    VarEq_InitVal(N^2+1:N^2+N) = X0;
    % Integrate
    options  = odeset('Reltol',1e-14,'Abstol',1e-16); %G_var.IntFunc.ODEoptions;
    fun      = @(t,x)VarEqAndSTMDOT(t,x,MU);          %G_var.IntFunc.VarEqAndSTMdot;
    tspan    = [0 tend];
    [t,xInteg] = ode45(fun,tspan,VarEq_InitVal,options);
    % Final output values
    PHI = xInteg;
    PHItf = reshape(xInteg(end,1:36),6,6);
    x = xInteg(:,37:42);
    xf = xInteg(end,37:42);
end