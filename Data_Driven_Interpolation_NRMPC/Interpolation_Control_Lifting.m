clc
clear
close all
% rng(10)
% rng(88)
r = randi(20000)
rng(r)
% A = {};
% B = {};
% 
% %% System Description
% 
% 

% Duff
f_u =  @(t,x,u)([x(2,:) ; -x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );
% f_u =  @(t,x,u)([ x(2,:) ; -0.5*x(2, :) + 1 * x(1, :) - 4 * x(1, :).^3 + u] );
% 
% % VDP
% % f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 3*x(1,:).^2.*x(2,:) - 1*x(2,:) - u; ] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
% f_u =  @(t,x,u)(-[ -x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
% f_u =  @(t,x,u)([ x(2,:) ; x(1,:) - 0.5 * (1 - x(1, :).^2) .* x(2, :) + (x(1, :) + 1).* u; ] );
% f_u =  @(t,x,u)([ x(2,:) ; x(1,:) + 0.5 * (1 - x(1, :).^2) .* x(2, :) + (abs(x(2, :)) + 1).* u; ] );
% f_u =  @(t,x,u)([ x(2,:) ; x(1,:) - 0.5 * (1 - x(1, :).^2) .* x(2, :) + (x(1, :) ).* u; ] );
% f_u =  @(t,x,u)([ x(2,:) ; x(1,:) - 0.5 * (1 - x(1, :).^2) .* x(2, :) + (x(1, :) ).* u; ] );


% % Inverted Pendulm
% MM = 0.0762;
% g = 9.80;
% J = 2.44e-4;
% Km = 10.51;
% l = 0.041;
% tau = 0.398;
% f_u =  @(t,x,u)([ x(2,:) ;  MM * g * l / J * sin(x(1,:)) - 1 / tau * x(2, :) + Km / tau * u] );

% % RotationalRobotic Manipulator System
% MM = 0.0292;
% Tc = 0.416;
% Ts = 0.4657;
% kappa = 1000;
% vs = 0.2;
% sigma = 0.0135;
% BB = 16;
% Cc = @(x) Tc * 2 / pi * atan(kappa * x);
% Cs = @(x) (Ts - Tc) * exp(- (x / vs)^2) * 2 / pi * atan(kappa * x);
% Cv = @(x) sigma * x;
% f_u =  @(t,x,u)([ x(2,:) ;  -inv(MM) * (Cc(x(2, :)) + Cs(x(2, :)) + Cv(x(2, :))) + inv(MM) * BB * u] );



deltaT = 0.05;
n = 2;
m = 1;
Nd = 20;
% Runge-Kutta 4

k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );

Phi = @(x) [x; x(1) / 100 + x(2)  ];
% Phi = @(x) [x; exp(-abs(sin(x(2)))) - 1];

NK = size(Phi(zeros(n, 1)), 1);
m = 1;
%% Koopman informed Initial Value of the Model
Total_Times = 1000;
Koopman_AB = {};
Continuous_Koopman_A = {};
Continuous_Koopman_B = {};
Start_points = [];
for times = 1 : Total_Times
    Koopman_X = [];
    Koopman_U = [];
    Koopman_Y = [];
%     Ns = randi([10, 20]);
    Ns = 4;
    Nt = 20;
    segment = 4 / Total_Times;
    start_points = -2 + (times - 1) * segment;
    Start_points = [Start_points start_points];
    for k = 1 : Nt
        x = start_points + 2 * segment * rand(n, 1) - segment;
        x_lift = Phi(x);
        for i = 1 : Ns
            u = 4 * rand - 2;
            x_next = f_ud(0, x, u);
            x_next_lift = Phi(x_next);
            Koopman_X = [Koopman_X x_lift];
            Koopman_Y = [Koopman_Y x_next_lift];
            Koopman_U = [Koopman_U u];
            x = x_next;
            x_lift = Phi(x);
        end
    end
    Koopman_AB{times} = Koopman_Y * [Koopman_X; Koopman_U]' * pinv([Koopman_X; Koopman_U] * [Koopman_X; Koopman_U]')
%     Continuous_Koopman_A{times} = (Koopman_AB{times}(1 : n, 1 : n) - eye(n, n)) / deltaT;
%     Continuous_Koopman_B{times} = Koopman_AB{times}(1 : n, n + 1 : n + m) / deltaT;
end

% minus_Koopman_AB = {};
% 
% for i = 1 : Total_Times - 1
%     for j = i + 1 : Total_Times
%         minus_Koopman_AB{i} = Koopman_AB{i} - Koopman_AB{j};
%     end
% end

LDI_Koopman_Max = zeros(size(Koopman_AB{1}));
LDI_Koopman_Min = zeros(size(Koopman_AB{1}));

for j = 1 : NK
    for k = 1 : NK + m
        max_element = -1e9;
        min_element = 1e9;
        for i = 1 : Total_Times
            if max_element < Koopman_AB{i}(j, k)
                max_element = Koopman_AB{i}(j, k);
            end
            if min_element > Koopman_AB{i}(j, k)
                min_element = Koopman_AB{i}(j, k);
            end
        end
        LDI_Koopman_Max(j, k) = max_element;
        LDI_Koopman_Min(j, k) = min_element;
    end
end

LDI_Koopman_Max
LDI_Koopman_Min

control_epsilon = 0.01;

rank(ctrb(LDI_Koopman_Min(1:NK, 1:NK), LDI_Koopman_Min(1:NK, NK + 1: NK + m)))



A = {};
B = {};
% for k = 1 : 2
%     for i = 1 : n
%         for j = 1 : n + m
%             if j <= n
%                 if k == 1 
%                     A{i}(i, j) = LDI_Koopman_Max(i, j);
%                 end
%                 if k == 2
%                     A{i}(i, j) = LDI_Koopman_Min(i, j);
%                 end
%             else
%                 if k == 1 
%                     B{i}(i, j) = LDI_Koopman_Max(i, j);
%                 end
%                 if k == 2 
%                     B{i}(i, j) = LDI_Koopman_Min(i, j);
%                 end
%             end
%         end
%     end
% end


A = {};
for k = 1 : Nd
    for i = 1 : NK
        for j = 1 : NK
            A{k}(i, j) = LDI_Koopman_Min(i, j) + (LDI_Koopman_Max(i, j) - LDI_Koopman_Min(i, j)) * randi([0, 1]);
        end
    end
end


B = {};
for k = 1 : Nd
    for i = 1 : NK
        for j = NK + 1 : NK + m
            B{k}(i, j - NK) = LDI_Koopman_Min(i, j) + (LDI_Koopman_Max(i, j) - LDI_Koopman_Min(i, j)) * randi([0, 1]);
        end
    end
end
% A{Nd + 1} = diag(0 * ones(n, 1));
% B{Nd + 1} = 0 * ones(n, 1);

Nd_temp = 2 * Nd;

for i = Nd + 1 : Nd_temp
    A{i} = -A{i - Nd};
    B{i} = -B{i - Nd};
end


A{Nd_temp + 1} = diag(0 * ones(NK, 1));
B{Nd_temp + 1} = 0 * ones(NK, 1);
Nd = Nd_temp + 1;

AB = {};
for i = 1 : Nd
    AB{i} = [A{i} B{i}];
end

%% Constraints Description
% load CDC2024-2.mat
C = [eye(n, n), zeros(n, NK - n)];
n = size(B{1}, 1);
m = size(B{1}, 2);
Control_Bound = 4;
State_Bound = 2;

F_x = [-eye(n, n);
       eye(n, n)];

G_x = [State_Bound;
       State_Bound;
       State_Bound;
       State_Bound;
       State_Bound;
       State_Bound];
F_u = [-eye(m, m);
       eye(m, m)];
G_u = [Control_Bound;
       Control_Bound];
cx = size(F_x, 1);
cu = size(F_u, 1);

%% Controller Design Offline
P = {};
Y = {};
Z = {};

mu = sdpvar(Nd, 1);
d = 4 * rand(n, Nd) - 2;
% d = kron(ones(1, Nd), [1.67; -1.41]);
% d = kron(ones(1, Nd), [-1; 0]);

for i = 1 : Nd
    P{i} = sdpvar(n, n);
    Y{i} = sdpvar(m, n);
    Z{i} = sdpvar(cx, cx);
    G{i} = sdpvar(cu, cu);
end

Constraints = {};
for i = 1 : Nd
    for j = 1 : Nd
        LMI_Invariance = [P{i}  (A{i} * P{j} + B{i} * Y{j});
                          (A{i} * P{j} + B{i} * Y{j})' P{j}];
        Constraints = [Constraints;
                       LMI_Invariance >= 0;];
    end
    LMI_Constraints_x = [Z{i} F_x * P{i};
                        (F_x * P{i})' P{i}];
    LMI_Constraints_u = [G{i} F_u * Y{i};
                        (F_u * Y{i})' P{i}];
    LMI_Maximum_Invariance = [1 mu(i) * d(:, i)';
                              mu(i) * d(:, i) P{i}];
    Constraints = [Constraints;
                   LMI_Constraints_x >= 0;
                   LMI_Constraints_u >= 0;
                   LMI_Maximum_Invariance >= 0];
    for k = 1 : cx
        Constraints = [Constraints; 
                       Z{i}(k, k) <= G_x(k)^2;];
    end

    for k = 1 : cu
        Constraints = [Constraints; 
                       G{i}(k, k) <= G_u(k)^2;];
    end
end

Objective = -ones(1, Nd) * mu;
options = sdpsettings('verbose',1,'solver','mosek');
sol = optimize(Constraints,Objective,options)

if sol.problem ~= 0
    error("Infeasible!")
end
mu_s_interpolation = double(ones(1, Nd) * mu);
mu_interpolation = double(mu);

for i = 1 : Nd
    P{i} = double(P{i});
    Y{i} = double(Y{i});
    K{i} = Y{i} * inv(P{i});
end

Upper_Bound_P = sdpvar(n * Nd, n * Nd);
K_z = [];
K_e = [];
for i = 1 : Nd
    if i == 1
        K_z = [K_z K{i}];
    else
        K_e = [K_e K{i} - K{1}];
        K_z = [K_z K{i} - K{1}];
    end
end
Q = 100 * eye(n, n);
Q1 = [eye(n, n) zeros(n, n * Nd - n)]' * Q * [eye(n, n) zeros(n, n * Nd - n)];
R = 0.001;
R1 = K_z' * R * K_z;
LMI_UpperBound = {};
Phi_z = {};
for i = 1 : Nd
    Temp = [];
    for j = 1 : Nd
        if j == 1
            RowVector = [A{1} + B{1} * K{1}, B{1} * K_e];
            Temp = [Temp; RowVector];
        else
            ej = eye(Nd, Nd);
            ej = ej(j, :);
            RowVector = kron(ej, A{j} + B{j} * K{j});
            Temp = [Temp; RowVector];
        end
    end
    Phi_z{i} = Temp;
end
Constraints = [];
for i = 1 : Nd
    LMI_UpperBound{i} = Upper_Bound_P - Q1 - R1 - Phi_z{i}'*Upper_Bound_P*Phi_z{i};
    Constraints = [Constraints LMI_UpperBound{i} >= 0];
end

sol = optimize(Constraints,1,options);

if sol.problem ~= 0
    error("Infeasible!")
end
Upper_Bound_P = double(Upper_Bound_P);

InvariantF = figure;
for i = 1 : Nd
    R = chol(inv(C*P{i}*C'));
    t = 0 : 0.001 : 2 * pi + 0.3;
    XX = inv(R) * [cos(t); sin(t)];
    plot(XX(1, :), XX(2, :));
    hold on
    grid on
end

% CP = [eye(n, n) eye(n, n * Nd - n)];
% R_Static = chol(inv(CP * Upper_Bound_P * CP'));
% t = 0 : 0.001 : 2 * pi + 0.3;
% XX = inv(R_Static) * [cos(t); sin(t)];
% plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');



[P_all_equal, mu_s_all_equal, mu_all_equal, K_equal_gain, second_P] = all_equal_gain(A(1: Nd), B(1: Nd), d);
R_Static = chol(inv(C*P_all_equal*C'));
t = 0 : 0.001 : 2 * pi + 0.3;
XX = inv(R_Static) * [cos(t); sin(t)];
plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');

% R_Static = chol(second_P);
% t = 0 : 0.001 : 2 * pi + 0.3;
% XX = inv(R_Static) * [cos(t); sin(t)];
% plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');

while 1
    xx = 2 * rand(2, 1) -1;
    if xx'* inv(C*P_all_equal*C') * xx <= 1 && Phi(xx)' * inv(P_all_equal) * Phi(xx) <= 1 
        break;
    end
end
% x = [1; -1];
[K_Static, nosolution_flag, Q, gamma, P_Static, mu_static] = getOptimalGain(Phi([0.6; -0.8]), A(1: Nd), B(1: Nd), d);
if nosolution_flag ~= 0
    error("Infeasible!")
end
mu_s_static = ones(1, Nd) * mu_static


Terminal_Invariant = inv(Q);
R_Static = chol(C*inv(Q)*C');
t = 0 : 0.001 : 2 * pi + 0.3;
XX = inv(R_Static) * [cos(t); sin(t)];
plot(XX(1, :), XX(2, :), 'LineWidth', 2.0, 'LineStyle','-.');





%% Simlulation
N = 5000;
% x_init = [4.3;
%           -3];
% x_init = [-.7;
%            3];
x_init = [.8;
          -1];
% x_init = [-1.3; 1.8];
% x_init = 2 * rand(2, 1) - 1;
x = Phi(x_init);
X = [C*x];
U = [];
x_rec = x;
X_recon = [C*x_rec];
Lambda = [];
P_i = [];
M_Set = {};
V_Set = [];
for i = 1 : N
%     A_k = zeros(n, n);
%     B_k = zeros(n, m);
%     p(1) = rand;
%     p(2) = 1 - p(1);
%     for j = 1 : Nd
%         A_k = A_k + p(j) * A{j};
%         B_k = B_k + p(j) * B{j};
%     end

%     if x'*inv(P_all_equal)*x <= 1
%         u = K_equal_gain * x;
%         V = x'*P_all_equal*x
%     gamma - x'*second_P*x
%     if 0
    if x'*Terminal_Invariant*x <= 1
        u = K_Static * x;
        V = x'*P_Static*x
    else
        [lambda, xi, cost, vi] = Get_Online_Control_CVX(x, A, P, Upper_Bound_P);
        u_Decomposition = zeros(m, Nd);
        z = [];
        for j = 1 : Nd
            u_Decomposition(:, j) = K{j} * xi(:, j);
            if j == 1
                z = [z x'];
            else
                z = [z xi(:, j)'];  
            end
        end
        V = cost + gamma
        u = u_Decomposition * ones(Nd, 1);
    end
    V_Set = [V_Set V];
%     u = 2 * rand - 1;
%     u = K{3} * x;
%     X(:, i + 1) = A_k * X(:, i) + B_k * u;
    for j = 1 : Nd
        x_candidate{j} = A{j} * x_rec + B{j} * u;
        AB{j} = vec([A{j} B{j}]);
    end
    X(:, i + 1) = f_ud(0, X(:, i), u);
    M = cell2mat(x_candidate);
    M_Set{i} = M;
    p_i = lsqlin(M, Phi(X(:, i + 1)), -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
    P_i = [P_i p_i];
    A_Nd = cell2mat(A);
    B_Nd = cell2mat(B);
    A_Rec = A_Nd * kron(p_i, eye(n));
    B_Rec = B_Nd * p_i;
    ERROR = M * p_i - Phi(X(:, i + 1));
    x_rec = A_Rec * x_rec + B_Rec * u ;
    X_recon = [X_recon C*x_rec];
    x = Phi(X(:, i + 1))
    U = [U u];
%     Lambda = [Lambda lambda];
end

tspan = deltaT * (1 : size(X_recon, 2));

figure
plot(tspan, X(1, (1 : size(X_recon, 2))), 'LineStyle','--','LineWidth', 2.0);
hold on 
plot(tspan, X_recon(1, (1 : size(X_recon, 2))), 'LineStyle','-','LineWidth', 1.0);
legend('x_1', 'reconstruncted x_1')


figure
plot(tspan, X(2, (1 : size(X_recon, 2))), 'LineStyle','--','LineWidth', 2.0);
hold on 
plot(tspan, X_recon(2, (1 : size(X_recon, 2))), 'LineStyle','-','LineWidth', 1.0);

legend('x_2', 'reconstruncted x_2')

% figure
% for i = 1 : size(P_i, 2)
%     if mod(i, 1) == 0
%         PP = M_Set{i}';
%         vertices = convhull(PP);
%         plot3(i * ones(size(PP(vertices, 1), 1)), PP(vertices, 1), PP(vertices, 2), 'LineWidth', 2)
%     end
%     hold on
%     scatter3(i, X(1, i + 1), X(2, i + 1), 'black', 'filled')
%     xlabel('Steps');
%     ylabel('$x_1$', 'Interpreter','latex');
%     zlabel('$x_2$', 'Interpreter','latex');
%     pause(0.1)
% end

figure
plot(U);

figure(InvariantF)
h1 = plot(X(1, :), X(2, :), 'LineWidth', 3.0);
legend([h1], 'Phase Trajectory')


figure
plot(V_Set, 'LineWidth', 3.0);

% %% 
% N = 1000;
% x_init = [0.6;
%           -1.0];
% x = x_init;
% X = [x];
% U = [];
% Lambda = [];
% for i = 1 : N
%     if x' * Terminal_Invariant * x <= 1
%         u = K_Static * x;
%     else 
%         [lambda, vi] = Get_Online_Control_CVX(x, A, P);
%         u_Decomposition = zeros(m, Nd);
%         for j = 1 : Nd
%             u_Decomposition(:, j) = K{j} * vi(:, j);
%         end
%         u = u_Decomposition * lambda;
%     end
%     X(:, i + 1) = f_ud(0, X(:, i), u);
%     x = X(:, i + 1)
%     U = [U u];
%     Lambda = [Lambda lambda];
% end
% 
% figure
% plot(X(1, :));
% hold on 
% plot(X(2, :));
% 
% figure
% plot(U);
% 
% figure
% S_Lambda = sum(Lambda, 1);
% plot(S_Lambda);




%% Controller Design Online
function [lambda, xi, cost, vi] = Get_Online_Control_CVX(x, A, P, UpperBound_P)
    scatter(x(1, :), x(2, :))

    Nd = size(A, 2);
    n = size(A{1}, 1);
    cvx_clear  
    cvx_begin sdp
    
        cvx_solver('Mosek_5')
%         cvx_solver('SeDuMi_3')
        variable lambda(Nd, 1)
        variable xi(n, Nd)
        subject to
        for i = 1 : Nd
            LMI_vi = [lambda(i) xi(:, i)'
                      xi(:, i)  lambda(i) * P{i}];
            LMI_vi >= 0;
            lambda(i) >= 0;
        end
        xi * ones(Nd, 1) == x;
        lambda' * ones(Nd, 1) <= 1;
        minimize([x; vec(xi(:, 2 : Nd))]' * UpperBound_P * [x; vec(xi(:, 2 : Nd))] + lambda'*lambda)
    cvx_end


    if cvx_status ~= "Solved"
        UpperBound_P = double(UpperBound_P);
        norm(UpperBound_P)
        error('Infeasible!')
    end
    pause(0.1)
    
    lambda = double(lambda);
    xi = double(xi);
    cost = [x; vec(xi(:, 2 : Nd))]' * UpperBound_P * [x; vec(xi(:, 2 : Nd))] + lambda'*lambda;
    vi = zeros(n, Nd);
    for i = 1 : Nd
        if lambda(i) >= 1e-6
          vi(:, i) = xi(:, i) ./ lambda(i);
        end
    end
    scatter(vi(1, :), vi(2, :))
end
























