clc
clear
close all
r = randi(2000)
% 1147
% 748 good 1 step
% 368 good 3 steps
% rng(823)
% rng(1332)
rng(1243)
n = 2;
m = 1;
Nd = 10;
NK = n;
Det = [];
Eig = [];
P_i = [];


% f_u =  @(t,x,u)([ x(2,:) ; 4 * 9.8 * sin(x(1,:)) + 3 * u * cos(x(1, :))] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 3*x(1,:).^2.*x(2,:) - 1*x(2,:) - u] );
% f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 2*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);

% % Inverted Pendulm
% MM = 0.0762;
% g = 9.80;
% J = 2.44e-4;
% Km = 10.51;
% l = 0.041;
% tau = 0.398;
% f_u =  @(t,x,u)([ x(2,:) ;  MM * g * l / J * sin(x(1,:)) - 1 / tau * x(2, :) + Km / tau * u] );
% f_u =  @(t,x,u)(-[ -x(2,:) ; 1*x(1,:) + 1*x(1,:).^2.*x(2,:) - 1 * x(2,:) - u]);
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);

% RotationalRobotic Manipulator System


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




% Duff
f_u =  @(t,x,u)([ x(2,:) ; -1*x(2, :) + 2 * x(1, :) - 2 * x(1, :).^3 + u] );
% VDP
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 0.8*x(1,:) + 3*x(1,:).^2.*x(2,:) - 1*x(2,:) - u] );
% f_u =  @(t,x,u)(-[ -2*x(2,:) ; 1*x(1,:) + 3*x(1,:).^2.*x(2,:) - 0.8*x(2,:) - u]);
% f_u =  @(t,x,u)(-[ -1*x(2,:) ; 1*x(1,:) + 2*x(1,:).^2.*x(2,:) - 2*x(2,:) - u]);

deltaT = 0.05;
%Runge-Kutta 4
k1 = @(t,x,u) (  f_u(t,x,u) );
k2 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT/2,u) );
k3 = @(t,x,u) ( f_u(t,x + k2(t,x,u)*deltaT/2,u) );
k4 = @(t,x,u) ( f_u(t,x + k1(t,x,u)*deltaT,u) );
f_ud = @(t,x,u) ( x + (deltaT/6) * ( k1(t,x,u) + 2*k2(t,x,u) + 2*k3(t,x,u) + k4(t,x,u)  )   );


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
    Ns = 3;
    Nt = 20;
    segment = 4 / Total_Times;
    start_points = -2 + (times - 1) * segment;
    Start_points = [Start_points start_points];
    for k = 1 : Nt
        x = start_points + 2 * segment * rand(n, 1) - segment;
        for i = 1 : Ns
            u = 40 * rand - 20;
            x_next = f_ud(0, x, u);
            Koopman_X = [Koopman_X x];
            Koopman_Y = [Koopman_Y x_next];
            Koopman_U = [Koopman_U u];
            x = x_next;
        end
    end
    Koopman_AB{times} = Koopman_Y * [Koopman_X; Koopman_U]' * pinv([Koopman_X; Koopman_U] * [Koopman_X; Koopman_U]')
    Continuous_Koopman_A{times} = (Koopman_AB{times}(1 : n, 1 : n) - eye(n, n)) / deltaT;
    Continuous_Koopman_B{times} = Koopman_AB{times}(1 : n, n + 1 : n + m) / deltaT;
end

minus_Koopman_AB = {};

for i = 1 : Total_Times - 1
    for j = i + 1 : Total_Times
        minus_Koopman_AB{i} = Koopman_AB{i} - Koopman_AB{j};
    end
end

LDI_Koopman_Max = zeros(size(Koopman_AB{1}));
LDI_Koopman_Min = zeros(size(Koopman_AB{1}));


for j = 1 : n
    for k = 1 : n + m
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
A = {};
B = {};


A = {};
for k = 1 : Nd
    for i = 1 : n
        for j = 1 : n
            A{k}(i, j) = LDI_Koopman_Min(i, j) + (LDI_Koopman_Max(i, j) - LDI_Koopman_Min(i, j)) * randi([0, 1]);
        end
    end
end


B = {};
for k = 1 : Nd
    for i = 1 : n
        for j = n + 1 : n + m
            B{k}(i, j - n) = LDI_Koopman_Min(i, j) + (LDI_Koopman_Max(i, j) - LDI_Koopman_Min(i, j)) * randi([0, 1]);
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


A{Nd_temp + 1} = diag(0 * ones(n, 1));
B{Nd_temp + 1} = 0 * ones(n, 1);
Nd = Nd_temp + 1;

AB = {};
for i = 1 : Nd
    AB{i} = [A{i} B{i}];
end


%% Polytopic Modelling Validation

% Nd = size(A, 2);
e = [];
index_identification = 1;

epsilon = 0;
eta = 1.00;
cnt = 0;
cnt_OK = 0;
N = 50;
states = 50;
Epsilon = [];
X_Set = 4 * rand(n, states) - 2;
fprintf('Progress:\n');

while cnt_OK < 3
    epsilon = 0;
    e = [];
    for k = 1 : states
        Tspan = (1 : N) * deltaT;
        
%         x = 4 * rand(2, 1) - 2;
        x = X_Set(:, k);
        x_rec = x;
        X_test = [x];
        X_recon_test = [x_rec];
        U = [];


%         figure
        for i = 1 : N
            w = 0;
            u = 4 * rand - 2;
            for j = 1 : Nd
                x_candidate{j} = A{j} * x_rec + B{j} * u;
            end
            w = 0;
            x_next = f_ud(0, x, u) + w;
        %     M = cell2mat(AB);
        %     M = Plot_Mapping;
            M = cell2mat(x_candidate);
            p_i = lsqlin(M, x_next, -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
            epsilon = epsilon + M * p_i - x_next;
        %     p_i = lsqlin(M, vec(A_B), -eye(size(M, 2), size(M, 2)), zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
            
%             if mod(i, 10) == 0
%                 P = M';
%                 vertices = convhull(P);
%                 plot3(i * ones(size(P(vertices, 1), 1)), P(vertices, 1), P(vertices, 2), 'LineWidth', 2)
%             end
%             hold on
%             scatter3(i, x_next(1), x_next(2), 'black', 'filled')
%             xlabel('Steps');
%             ylabel('$x_1$', 'Interpreter','latex');
%             zlabel('$x_2$', 'Interpreter','latex');
            
    %         pause(0.05)
        
            P_i = [P_i p_i];
            A_Nd = cell2mat(A);
            B_Nd = cell2mat(B);
            A_Rec = A_Nd * kron(p_i, eye(n));
            B_Rec = B_Nd * p_i;
            e = [e M * p_i - x_next];
        
        %     A_B_Rec = reshape(M * p_i, [n, n + m]);
        %     A_Rec = A_B_Rec(1 : n, 1 : n);
        %     B_Rec = A_B_Rec(1 : n, n + 1 : n + m);
        %     e = [e M * p_i - vec(A_B)]
        
            A_Rec_Set{i} = A_Rec;
            B_Rec_Set{i} = B_Rec;
            x_rec = A_Rec * x_rec + B_Rec * u ;
            X_test = [X_test x_next];
            X_recon_test = [X_recon_test x_rec];
            U = [U u];
            x = x_next;
        end
        percent = k / states * 100;
        fprintf('\r%.0f%% completed', percent);
    end
    total_error = norm(epsilon)
    max_e = norm(e', inf)
    Epsilon = [Epsilon max_e];
    if norm(e', inf) <=  0.1
       cnt_OK = cnt_OK + 1;
    else
       eta = eta + 0.001;
    end

    cnt = cnt + 1;
 
    
    % CDC Approach
    randomized_index = randi((Nd - 1) / 2);
    AB{randomized_index} = eta * AB{randomized_index};
    
    

    for i = (Nd - 1) / 2 + 1 : Nd - 1
        AB{i} = -1 * AB{i - (Nd - 1) / 2};
    end

    for i = 1 : Nd
        A{i} = 1 * AB{i}(1 : n, 1 : n);
        B{i} = 1 * AB{i}(1 : n, n + 1 : n + m);
    end
        



    figure
    mesh(deltaT * (1 : size(P_i, 2)), 1 : size(M, 2), P_i)
    xlabel('Time(sec)')
    ylabel('i-th vertex')
    zlabel('p_i')
    
    
    
    figure
    plot(X_test(1, :))
    hold on 
    plot(X_recon_test(1, :))
    legend('Original', 'LPV')
    ylabel('$x_1$', 'Interpreter','latex')
    
    
    figure
    plot(X_test(2, :))
    hold on 
    plot(X_recon_test(2, :))
    legend('Original', 'LPV')
    ylabel('$x_2$', 'Interpreter','latex')
    

    
    now_error = norm(X_test - X_recon_test) / N

end

figure
plot(Epsilon)
ylabel('$\epsilon$', 'Interpreter','latex')



%% Control - Interpolation-Termination control  %% 

%% Constraints Description
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
Q = 1000 * eye(n, n);
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
figure(InvariantF)
hold on
for i = 1 : Nd
    R = chol(inv(P{i}));
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
figure(InvariantF)
hold on
R_Static = chol(inv(P_all_equal));
t = 0 : 0.001 : 2 * pi + 0.3;
XX = inv(R_Static) * [cos(t); sin(t)];
plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');

% R_Static = chol(second_P);
% t = 0 : 0.001 : 2 * pi + 0.3;
% XX = inv(R_Static) * [cos(t); sin(t)];
% plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');

x_max_index = find(sum(XX.^2, 1) == max(sum(XX.^2, 1)))
% x = [1; -1];
[K_Static, nosolution_flag, Q, gamma, P_Static, mu_static] = getOptimalGain(XX(:, x_max_index), A(1: Nd), B(1: Nd), d);
if nosolution_flag ~= 0
    error("Infeasible!")
end
mu_s_static = ones(1, Nd) * mu_static;


Terminal_Invariant = inv(Q);

figure(InvariantF)
hold on
R_Static = chol(inv(Q));
t = 0 : 0.001 : 2 * pi + 0.3;
XX = inv(R_Static) * [cos(t); sin(t)];
plot(XX(1, :), XX(2, :), 'LineWidth', 2.0, 'LineStyle','-.');





%% Simlulation
N = 2000;
% x_init = [4.3;
%           -3];
% x_init = [-.7;
%            3];
x_init = [1.0;
          -1];
% x_init = [-1.3; 1.8];
% x_init = 2 * rand(2, 1) - 1;
x = x_init;
X = [x];
U = [];
x_rec = x;
X_recon = [x_rec];
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
    if x'*Terminal_Invariant*x <= 1
        u = K_Static * x;
        V = x'*P_Static*x
    else

        [lambda, xi, cost] = Get_Online_Control_CVX(x, A, P, Upper_Bound_P, InvariantF);

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
    V_Set = [V_Set V]
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
    p_i = lsqlin(M, X(:, i + 1), -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
    P_i = [P_i p_i];
    A_Nd = cell2mat(A);
    B_Nd = cell2mat(B);
    A_Rec = A_Nd * kron(p_i, eye(n));
    B_Rec = B_Nd * p_i;
    ERROR = M * p_i - X(:, i + 1);
    x_rec = A_Rec * x_rec + B_Rec * u ;
    X_recon = [X_recon x_rec];
    x = X(:, i + 1)
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
hold on
h1 = plot(X(1, :), X(2, :), 'LineWidth', 3.0);
legend([h1], 'Phase Trajectory')


figure
plot(V_Set, 'LineWidth', 3.0);





%% Controller Design Online
function [lambda, xi, cost] = Get_Online_Control_CVX(x, A, P, UpperBound_P, InvariantF)
    
    figure(InvariantF)
    hold on
    scatter(x(1, :), x(2, :))

    Nd = size(A, 2);
    n = size(A{1}, 1);

    cvx_begin sdp
    
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
        lambda' * ones(Nd, 1) == 1;
        minimize([x; vec(xi(:, 2 : Nd))]' * UpperBound_P * [x; vec(xi(:, 2 : Nd))] + lambda'*lambda)
    cvx_end

    if cvx_status ~= "Solved"
        error('infeasible')
    end
    
    lambda = double(lambda);
    xi = double(xi);
    cost = [x; vec(xi(:, 2 : Nd))]' * UpperBound_P * [x; vec(xi(:, 2 : Nd))] + lambda'*lambda;
%     vi = zeros(n, Nd);
%     for i = 1 : Nd
%         if lambda(i) <= 1e-5
%             vi(:, i) = xi(:, i).*0;
%         else
%             vi(:, i) = xi(:, i) ./ lambda(i);
%         end
%     end
% 
%     figure(InvariantF)
%     hold on
%     scatter(vi(1, :), vi(2, :))
end















