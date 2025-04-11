function [P_Static,mu_s, mu, gain, second_P] = all_equal_gain(A, B, d)
    
    n = size(B{1}, 1);
    m = size(B{1}, 2);
    
    Control_Bound = 4;
    State_Bound = 2;
    
    F_x = [-eye(n, n);
           eye(n, n)];
    
    G_x = [State_Bound;
           State_Bound;
           State_Bound;
           State_Bound];
    F_u = [-eye(m, m);
           eye(m, m)];
    G_u = [Control_Bound;
           Control_Bound];
    cx = size(F_x, 1);
    cu = size(F_u, 1);


    Nd = size(B, 2);
    mu = sdpvar(Nd, 1);
%     d = rand(n, Nd);
    P = sdpvar(n, n);
    Y = sdpvar(m, n);
    Z = sdpvar(cx, cx);
    G = sdpvar(cu, cu);
    Constraints = [];


    for i = 1 : Nd
        LMI_Invariance = [P  (A{i} * P + B{i} * Y);
                          (A{i} * P + B{i} * Y)' P];
        Constraints = [Constraints;
                       LMI_Invariance >= 0;];
        for k = 1 : cx
            Constraints = [Constraints; 
                           Z(k, k) <= G_x(k)^2;];
        end
    
        LMI_Maximum_Invariance = [1 mu(i) * d(:, i)';
                                  mu(i) * d(:, i) P];
        Constraints = [Constraints;
                       LMI_Maximum_Invariance >= 0;
                       mu(i) >= 0];
        for k = 1 : cu
            Constraints = [Constraints; 
                           G(k, k) <= G_u(k)^2;];
        end
    end

    
    LMI_Constraints_x = [Z F_x * P;
                        (F_x * P)' P];
    LMI_Constraints_u = [G F_u * Y;
                        (F_u * Y)' P];

    Constraints = [Constraints; 
               LMI_Constraints_x >= 0;
               LMI_Constraints_u >= 0;];
    Objective = -ones(1, Nd) * mu;
    options = sdpsettings('verbose',1,'solver','mosek');
    sol = optimize(Constraints,Objective,options);
    
    if sol.problem ~= 0
        error("Infeasible!")
    end

    P_Static = double(P);
    mu_s = double(ones(1, Nd) * mu);
    Y = double(Y);
    mu = double(mu);
    gain = Y * inv(P_Static);





    Q1 = 100 * eye(n, n);
    R1 = 0.001;
    second_P = sdpvar(n, n);
    second_P_LMI = {};
    Constraints = [];
    for i = 1 : Nd
        second_P_LMI{i} = second_P - (A{i} + B{i} * gain)' * second_P * (A{i} + B{i} * gain) - Q1 - gain' * R1 * gain;
        Constraints = [Constraints; second_P_LMI{i} >= 0];

    end
    

    options = sdpsettings('verbose',1,'solver','mosek');
    sol = optimize(Constraints,1,options);
    
    if sol.problem ~= 0
        error("Infeasible!")
    end

    second_P = double(second_P);
    
