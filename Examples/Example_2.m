clc
clear
close all
r = randi(2000)
% 1147
% 748 good 1 step
% 368 good 3 steps
% rng(823)
% rng(1332)
% rng(1243) %dmd
rng(939)
% rng(r) %EDMD

n = 2;
m = 1;
Nd = 13;
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

Phi = @(x) [x; x(1) / 10 + x(2)  ];
% Phi = @(x) [x; exp(-abs(sin(x(2)))) - 1];
% Phi = @(x) [x  ];

NK = size(Phi(zeros(n, 1)), 1);
m = 1;

C = [eye(n,n ) zeros(n, NK - n)];
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
    Ns = NK + m;
    Nt = 20;
    segment = 4 / Total_Times;
    start_points = -2 + (times - 1) * segment;
    Start_points = [Start_points start_points];
    for k = 1 : Nt
        x = start_points + 2 * segment * rand(n, 1) - segment;
        x = Phi(x);
        for i = 1 : Ns
            u = 40 * rand - 20;
            x_next = f_ud(0, C* x, u);
            Koopman_X = [Koopman_X x];
            Koopman_Y = [Koopman_Y Phi(x_next)];
            Koopman_U = [Koopman_U u];
            x = Phi(x_next);
        end
    end
    Koopman_AB{times} = Koopman_Y * [Koopman_X; Koopman_U]' * pinv([Koopman_X; Koopman_U] * [Koopman_X; Koopman_U]')
    Continuous_Koopman_A{times} = (Koopman_AB{times}(1 : NK, 1 : NK) - eye(NK, NK)) / deltaT;
    Continuous_Koopman_B{times} = Koopman_AB{times}(1 : NK, NK + 1 : NK + m) / deltaT;
end

minus_Koopman_AB = {};

for i = 1 : Total_Times - 1
    for j = i + 1 : Total_Times
        minus_Koopman_AB{i} = Koopman_AB{i} - Koopman_AB{j};
    end
end

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
A = {};
B = {};


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
A{Nd + 1} = diag(0 * ones(n, 1));
B{Nd + 1} = 0 * ones(n, 1);

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


%% Polytopic Modelling Validation
% Nd = size(A, 2);
e = [];
index_identification = 1;

epsilon = 0;
eta = 1.0;
cnt = 0;
cnt_OK = 0;
N = 50;
states = 50;
Epsilon = [];
X_Set = 4 * rand(n, states) - 2;
fprintf('Progress:\n');
each_state_error_Set = [];
U_Set = [];
while cnt_OK <= 3
    epsilon = 0;
    e = [];
    U = [];
    each_state_error_Set = [];
    for k = 1 : states
        Tspan = (1 : N) * deltaT;
        
%         x = 4 * rand(2, 1) - 2;
        x = Phi(X_Set(:, k));
        x_rec = x;
        X_test = [x];
        X_recon_test = [x_rec];
        U = [];
        each_state_error = [];
%         figure
        for i = 1 : N
            w = 0;
            u = 12 * rand - 6;
            for j = 1 : Nd
                x_candidate{j} = A{j} * x_rec + B{j} * u;
            end
            w = 0;
            x_next = f_ud(0, C * x, u) + w;
        %     M = cell2mat(AB);
        %     M = Plot_Mapping;
            M = cell2mat(x_candidate);
            p_i = lsqlin(M, Phi(x_next), -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
            epsilon = epsilon + M * p_i - Phi(x_next);
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
            A_Rec = A_Nd * kron(p_i, eye(NK));
            B_Rec = B_Nd * p_i;
            e = [e M * p_i - Phi(x_next)];
            
            each_state_error = [each_state_error M * p_i - Phi(x_next)];
        %     A_B_Rec = reshape(M * p_i, [n, n + m]);
        %     A_Rec = A_B_Rec(1 : n, 1 : n);
        %     B_Rec = A_B_Rec(1 : n, n + 1 : n + m);
        %     e = [e M * p_i - vec(A_B)]
        
            A_Rec_Set{i} = A_Rec;
            B_Rec_Set{i} = B_Rec;
            x_rec = A_Rec * x_rec + B_Rec * u ;
            X_test = [X_test Phi(x_next)];
            X_recon_test = [X_recon_test x_rec];
            U = [U u];
            x = Phi(x_next);
        end
        percent = k / states * 100;
        fprintf('\r%.0f%% completed', percent);
        now_e = norm(each_state_error', inf);
        each_state_error_Set = [each_state_error_Set now_e];
        U_Set = [U_Set; U];
    end
    total_error = norm(epsilon)
    max_e = norm(e', inf)
    Epsilon = [Epsilon max_e];
    if cnt == 0
        max_index = find(each_state_error_Set == max(each_state_error_Set));
    end

    


    figure
    mesh(deltaT * (1 : size(P_i, 2)), 1 : size(M, 2), P_i)
    xlabel('Time(sec)')
    ylabel('i-th vertex')
    zlabel('p_i')
    




    x = Phi(X_Set(:, max_index));

    U_Test = U_Set(max_index, :);
    save('Test_EDMD.mat', 'x', 'U_Test');
%     x = Phi([0.04; -0.86]);  
%     x = Phi([-1.8091; 1.5093]);

    x_rec = x;
    X_test = [x];
%     load U_Test_EDMD.mat
%     load U_Test.mat
    X_recon_test = [x_rec];
    for i = 1 : N
        w = 0;
        u = U_Set(max_index, i);
%         u = U_Test_EDMD(i);
%         u = U_Test_EDMD(i);
        for j = 1 : Nd
            x_candidate{j} = A{j} * x_rec + B{j} * u;
        end
        w = 0;
        x_next = f_ud(0, C*x, u) + w;
        M = cell2mat(x_candidate);
        p_i = lsqlin(M, Phi(x_next), -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
        epsilon = epsilon + M * p_i - Phi(x_next);
        P_i = [P_i p_i];
        A_Nd = cell2mat(A);
        B_Nd = cell2mat(B);
        A_Rec = A_Nd * kron(p_i, eye(NK));
        B_Rec = B_Nd * p_i;
        e = [e M * p_i - Phi(x_next)];
    
    %     A_B_Rec = reshape(M * p_i, [n, n + m]);
    %     A_Rec = A_B_Rec(1 : n, 1 : n);
    %     B_Rec = A_B_Rec(1 : n, n + 1 : n + m);
    %     e = [e M * p_i - vec(A_B)]
    
        A_Rec_Set{i} = A_Rec;
        B_Rec_Set{i} = B_Rec;
        x_rec = A_Rec * x_rec + B_Rec * u ;
        X_test = [X_test Phi(x_next)];
        X_recon_test = [X_recon_test x_rec];
        x = Phi(x_next);
    end

    
% figure
% plot(tspan, X(1, (1 : size(X_recon, 2))), 'LineStyle','--','LineWidth', 2.0);
% hold on 
% plot(tspan, X_recon(1, (1 : size(X_recon, 2))), 'LineStyle','-','LineWidth', 1.0);
% legend('x_1', 'reconstruncted x_1')

    figure
    subplot(211)
    tspan = deltaT * (1 : size(X_test, 2));
    plot(tspan, X_test(1, :), 'LineStyle','--','LineWidth', 2.0)
    hold on 
    plot(tspan, X_recon_test(1, :), 'LineStyle','-','LineWidth', 1.0)
    legend('$x_1$', 'reconstructed $x_1$') 
    ylabel('$x_1$', 'Interpreter','latex')
    xlabel('$t/s$', 'Interpreter','latex')
    
    
    subplot(212)
    plot(tspan, X_test(2, :), 'LineStyle','--','LineWidth', 2.0)
    hold on 
    plot(tspan, X_recon_test(2, :), 'LineStyle','-','LineWidth', 1.0)
    legend('$x_2$', 'reconstructed $x_2$') 
    ylabel('$x_2$', 'Interpreter','latex')
    xlabel('$t/s$', 'Interpreter','latex')

    figure
    plot(tspan(1:end-1), U_Set(max_index, :), 'LineStyle','-','LineWidth', 2.0);
    ylabel('$u$', 'Interpreter','latex')
    xlabel('$t/s$', 'Interpreter','latex')


    now_error = norm(X_test - X_recon_test) / N




    if max_e <=  0.01
       cnt_OK = cnt_OK + 1;
    else
       % CDC Approach
        randomized_index = randi((Nd - 1) / 2);
        AB{randomized_index} = (eta + 0.01) * AB{randomized_index};
        
        
    
        for i = (Nd - 1) / 2 + 1 : Nd - 1
            AB{i} = -1 * AB{i - (Nd - 1) / 2};
        end
    
        for i = 1 : Nd
            A{i} = 1 * AB{i}(1 : NK, 1 : NK);
            B{i} = 1 * AB{i}(1 : NK, NK + 1 : NK + m);
        end
    end

    cnt = cnt + 1;
 
    
    
        

end

figure
plot(Epsilon, 'LineStyle','-','LineWidth', 2.0)
ylabel('$\Sigma_k \, ||\varepsilon_k||$', 'Interpreter','latex')
xlabel('Iterative step', 'Interpreter','latex')
grid on


%% Control - Interpolation-Termination control  %% 
n = NK;

%% Constraints Description
n = size(B{1}, 1);
m = size(B{1}, 2);
Control_Bound = 6;
State_Bound = 2;

F_x = [-eye(n, n);
       eye(n, n)];

G_x = [State_Bound;
       State_Bound;
       State_Bound;
       State_Bound;
       State_Bound / 10 + State_Bound;
       State_Bound / 10 + State_Bound];
F_u = [-eye(m, m);
       eye(m, m)];
G_u = [Control_Bound;
       Control_Bound];
cx = size(F_x, 1);
cu = size(F_u, 1);

%% Controller Design Offline
rng(939)
C = [eye(2, 2), zeros(2, NK - 2)];
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
K = {};
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
XX_Set = [];
hold on
for i = 1 : Nd
    R = chol(inv(P{i}));
    
    % 构造单位球面
    [x, y, z] = sphere(50); % 生成50x50的单位球面
    XYZ = [x(:) y(:) z(:)]';
    
    ellipsoid_points = inv(R) * XYZ;
    x_ellip = reshape(ellipsoid_points(1, :), size(x));
    y_ellip = reshape(ellipsoid_points(2, :), size(y));
    z_ellip = reshape(ellipsoid_points(3, :), size(z));

    XX_Set = [XX_Set ellipsoid_points];

%     % 绘图
%     surf(x_ellip, y_ellip, z_ellip,'FaceColor', [0.6, 1, 0.6], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'EdgeAlpha',  0.1);
%     axis equal
%     xlabel('x'); ylabel('y'); zlabel('z');
%     title('Region defined by X^T P X ≤ 1');
%     grid on
end


% 假设你有一个 Nx3 的三维点集
% 随机生成一些示例点
% 计算凸包
XX_Set = XX_Set';
[KK, volume] = convhull(XX_Set);

% 绘制凸包
convex_p = trisurf(KK, XX_Set(:,1), XX_Set(:,2), XX_Set(:,3), ...
    'FaceColor', [0.97, 0.75, 0.75], 'FaceAlpha', 0.2, 'EdgeColor', [0.97, 0.75, 0.75], 'EdgeAlpha',  0.5);
% [x_vals, y_vals] = meshgrid(-2:0.01:2, -2:0.01:2);
% z_vals = x_vals / 10 + y_vals;  % Plane equation z = x/10 + y
% surf(x_vals, y_vals, z_vals, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0, 0.2, 0.4], 'EdgeAlpha',  0.1);  % Plot the surface
valid_indices = abs(XX_Set(:,1) / 10 + XX_Set(:,2) - XX_Set(:,3)) <= 1e-2;

% 美化图形
axis equal
xlabel('$x_1$'); ylabel('$x_2$'); zlabel('$\Phi(x_1, x_2)$', 'Interpreter','latex');
title(sprintf('Convex hull of invariant sets (Volume = %.3f)', volume));
grid on



XX_Set_Indices = [XX_Set(valid_indices, 1) XX_Set(valid_indices, 2) XX_Set(valid_indices, 3)];

[KK, volume] = convhull(XX_Set_Indices);

% 绘制凸包
projection = trisurf(KK, XX_Set_Indices(:,1), XX_Set_Indices(:,2), XX_Set_Indices(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [0, 0.2, 0.4], 'EdgeAlpha',  0.1);



% CP = [eye(n, n) eye(n, n * Nd - n)];
% R_Static = chol(inv(CP * Upper_Bound_P * CP'));
% t = 0 : 0.001 : 2 * pi + 0.3;
% XX = inv(R_Static) * [cos(t); sin(t)];
% plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');



[P_all_equal, mu_s_all_equal, mu_all_equal, K_equal_gain, second_P] = all_equal_gain(A(1: Nd), B(1: Nd), d);
figure(InvariantF)
hold on
R_Static = chol(inv(P_all_equal));
[x, y, z] = sphere(50); % 生成50x50的单位球面
XYZ = [x(:) y(:) z(:)]';
ellipsoid_points = inv(R_Static) * XYZ;
% XX_Set = [XX_Set ellipsoid_points];
x_ellip = reshape(ellipsoid_points(1, :), size(x));
y_ellip = reshape(ellipsoid_points(2, :), size(y));
z_ellip = reshape(ellipsoid_points(3, :), size(z));


% plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');

% R_Static = chol(second_P);
% t = 0 : 0.001 : 2 * pi + 0.3;
% XX = inv(R_Static) * [cos(t); sin(t)];
% plot(XX(1, :), XX(2, :), 'LineWidth', 5.0, 'LineStyle','-.');

% x_max_index = find(sum(XX.^2, 1) == max(sum(XX.^2, 1)))
% Phi(XX(:, x_max_index)
% x = [1; -1];
[K_Static, nosolution_flag, Q, gamma, P_Static, mu_static] = getOptimalGain(Phi([0.7; -1.72]), A(1: Nd), B(1: Nd), d);
if nosolution_flag ~= 0
    error("Infeasible!")
end
mu_s_static = ones(1, Nd) * mu_static;


Terminal_Invariant = inv(Q);

figure(InvariantF)
hold on
R_Static = chol(inv(Q));

[x, y, z] = sphere(50); % 生成50x50的单位球面
XYZ = [x(:) y(:) z(:)]';
ellipsoid_points = inv(R_Static) * XYZ;
% XX_Set = [XX_Set ellipsoid_points];
x_ellip = reshape(ellipsoid_points(1, :), size(x));
y_ellip = reshape(ellipsoid_points(2, :), size(y));
z_ellip = reshape(ellipsoid_points(3, :), size(z));

% ellipsoid_points = ellipsoid_points';
% valid_indices = abs(ellipsoid_points(:,1) / 10 + ellipsoid_points(:,2) - ellipsoid_points(:,3)) <= 1e-3;
% XX_Set_Indices = [ellipsoid_points(valid_indices, 1) ellipsoid_points(valid_indices, 2) ellipsoid_points(valid_indices, 3)];
% 
% [KK, volume] = convhull(XX_Set_Indices);
% 
% % 绘制凸包
% projection_2 = trisurf(KK, XX_Set_Indices(:,1), XX_Set_Indices(:,2), XX_Set_Indices(:,3), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'FaceColor', [1, 0.5, 0.31], 'EdgeAlpha',  0.1);
 

% 绘图
terminal_invariant = surf(x_ellip, y_ellip, z_ellip, 'FaceAlpha', 0.2, 'EdgeColor', 'none'...
    , 'FaceColor', [0.83, 0.68, 0.22], 'EdgeAlpha', 0.5);
axis equal
grid on

%% Simlulation

X_init_Set = [1 -1.5;
              -1.2 0]';



U_Control_Input_Fig = figure;
Lyapunov_Fig = figure;
State_Each_Fig_1 = figure;
State_Each_Fig_2 = figure;
elapsedTime_Total = 0;
for times = 1 : 2
N = 500;
% x_init = [4.3;
%           -3];
% x_init = [-.7;
%            3];
% x_init = [1.42;
%           -.79];
% x_init = [-.94;
%           -.02];
x_init = X_init_Set(:, times);
% x_init = [-1.3; 1.8];
% x_init = 2 * rand(2, 1) - 1;
x = Phi(x_init);
X_3 = [C * x];
U_3 = [];
x_rec = x;
X_recon_3 = [x_rec];
X_Lift_Set_3 = [Phi(x_init)];
Lambda = [];
P_i = [];
M_Set = {};
V_Set_3 = [];
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
% 
    if 0
        u = K_Static * x;
        V = x'*P_Static*x
    else
        tic
        [lambda, xi, cost] = yalmip_online_interpolation(x, A, P, Upper_Bound_P);
        elapsedTime  = toc;
        elapsedTime_Total = elapsedTime_Total + elapsedTime;
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
        V = cost + gamma - gamma
        u = u_Decomposition * ones(Nd, 1);
    end
    V_Set_3 = [V_Set_3 V]
%     u = 2 * rand - 1;
%     u = K{3} * x;
%     X(:, i + 1) = A_k * X(:, i) + B_k * u;



    for j = 1 : Nd
        x_candidate{j} = A{j} * x_rec + B{j} * u;
        AB{j} = vec([A{j} B{j}]);
    end
    X_3(:, i + 1) = f_ud(0, X_3(:, i), u);
    M = cell2mat(x_candidate);
    M_Set{i} = M;
    p_i = lsqlin(M, Phi(X_3(:, i + 1)), -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
    P_i = [P_i p_i];
    A_Nd = cell2mat(A);
    B_Nd = cell2mat(B);
    A_Rec = A_Nd * kron(p_i, eye(n));
    B_Rec = B_Nd * p_i;
    ERROR = M * p_i - Phi(X_3(:, i + 1));
    x_rec = A_Rec * x_rec + B_Rec * u ;
    X_recon_3 = [X_recon_3 x_rec];
    x = Phi(X_3(:, i + 1))
    U_3 = [U_3 u];
    X_Lift_Set_3 = [X_Lift_Set_3 x];

%     Lambda = [Lambda lambda];
end

% tspan = deltaT * (1 : size(X_recon, 2));
% 
% figure
% plot(tspan, X_3(1, (1 : size(X_recon, 2))), 'LineStyle','--','LineWidth', 2.0);
% hold on 
% plot(tspan, X_recon_3(1, (1 : size(X_recon, 2))), 'LineStyle','-','LineWidth', 1.0);
% legend('x_1', 'reconstructed x_1')
% 
% 
% figure
% plot(tspan, X_3(2, (1 : size(X_recon, 2))), 'LineStyle','--','LineWidth', 2.0);
% hold on 
% plot(tspan, X_recon_3(2, (1 : size(X_recon, 2))), 'LineStyle','-','LineWidth', 1.0);
% 
% legend('x_2', 'reconstructed x_2')

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

% figure
% plot(U_3);

% figure(InvariantF)
% hold on
% h1 = plot3(X_Lift_Set_3(1, :), X_Lift_Set_3(2, :), X_Lift_Set_3(3, :), 'LineWidth', 3.0, 'Color','r');
% legend([convex_p terminal_invariant projection h1], {'Convex hull of the invariant sets', 'Interpolation-terminated invariant set', 'Projection of Koopman embedded observable ','Phase trajectory with interpolated control'}...
%     , 'Interpreter','latex');


% figure
% plot(V_Set_3, 'LineWidth', 3.0);
averaged_time = elapsedTime_Total / N

% Simlulation interpolation-terminated control
N = 500;
% x_init = [4.3;
%           -3];
% x_init = [-.7;
%            3];
% x_init = [1.42;
%           -.79];
% x_init = [-.94;
%           -.02];
x_init = X_init_Set(:, times);
% x_init = [-1.3; 1.8];
% x_init = 2 * rand(2, 1) - 1;
x = Phi(x_init);
X_4 = [C * x];
U_4 = [];
x_rec = x;
X_recon_4 = [x_rec];
X_Lift_Set_4 = [Phi(x_init)];
Lambda = [];
P_i = [];
M_Set = {};
V_Set_4 = [];
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
% 
    if x'*P_Static*x <= gamma
        u = K_Static * x;
        V = x'*P_Static*x
    else
        [lambda, xi, cost] = yalmip_online_interpolation(x, A, P, Upper_Bound_P);
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
    V_Set_4 = [V_Set_4 V]
%     u = 2 * rand - 1;
%     u = K{3} * x;
%     X(:, i + 1) = A_k * X(:, i) + B_k * u;



    for j = 1 : Nd
        x_candidate{j} = A{j} * x_rec + B{j} * u;
        AB{j} = vec([A{j} B{j}]);
    end
    X_4(:, i + 1) = f_ud(0, X_4(:, i), u);
    M = cell2mat(x_candidate);
    M_Set{i} = M;
    p_i = lsqlin(M, Phi(X_4(:, i + 1)), -eye(size(M, 2), size(M, 2))  , zeros(size(M, 2), 1),  ones(1, size(M, 2)), 1);
    P_i = [P_i p_i];
    A_Nd = cell2mat(A);
    B_Nd = cell2mat(B);
    A_Rec = A_Nd * kron(p_i, eye(n));
    B_Rec = B_Nd * p_i;
    ERROR = M * p_i - Phi(X_4(:, i + 1));
    x_rec = A_Rec * x_rec + B_Rec * u ;
    X_recon_4 = [X_recon_4 x_rec];
    x = Phi(X_4(:, i + 1))
    U_4 = [U_4 u];
    X_Lift_Set_4 = [X_Lift_Set_4 x];

%     Lambda = [Lambda lambda];
end

tspan = deltaT * (1 : size(X_recon_3, 2));

figure(State_Each_Fig_1)
subplot 211

plot(tspan, X_3(1, (1 : size(X_recon_3, 2))), 'LineStyle','--','LineWidth', 4.0);
hold on 
plot(tspan, X_recon_3(1, (1 : size(X_recon_3, 2))), 'LineStyle','-','LineWidth', 2.0);

subplot 212
plot(tspan, X_4(1, (1 : size(X_recon_4, 2))), 'LineStyle','--','LineWidth', 4.0);
hold on 
plot(tspan, X_recon_4(1, (1 : size(X_recon_4, 2))), 'LineStyle','-','LineWidth', 2.0);



figure(State_Each_Fig_2)
subplot 211
plot(tspan, X_3(2, (1 : size(X_recon_3, 2))), 'LineStyle','--','LineWidth', 4.0);
hold on 
plot(tspan, X_recon_3(2, (1 : size(X_recon_3, 2))), 'LineStyle','-','LineWidth', 2.0);

subplot 212
plot(tspan, X_4(2, (1 : size(X_recon_4, 2))), 'LineStyle','--','LineWidth', 4.0);
hold on 
plot(tspan, X_recon_4(2, (1 : size(X_recon_4, 2))), 'LineStyle','-','LineWidth', 2.0);

% legend('$x_2$ using Algorithm 4', 'reconstructed $x_{2,rec}$ using Algorithm 4')

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


figure(U_Control_Input_Fig)
% subplot 211
xlabel('$t/s$', 'Interpreter','latex')
ylabel('Control input $u$', 'Interpreter','latex')
plot(tspan(1:end - 1), U_3(1:end), 'LineWidth', 3.0)
hold on

% subplot 212
plot(tspan(1:end - 1), U_4(1:end), 'LineWidth', 3.0)
grid on
xlabel('$t/s$', 'Interpreter','latex')
ylabel('Control input $u$', 'Interpreter','latex')
hold on

figure(InvariantF)
hold on
h1 = plot3(X_Lift_Set_3(1, :), X_Lift_Set_3(2, :), X_Lift_Set_3(3, :), 'LineWidth', 3.0, 'Color', 'r');
h2 = plot3(X_Lift_Set_4(1, :), X_Lift_Set_4(2, :), X_Lift_Set_4(3, :), 'LineWidth', 3.0, 'LineStyle','--', 'Color', 'b');
legend([convex_p terminal_invariant projection h1 h2], {'Convex hull of the invariant sets', 'Interpolation-terminated invariant set', 'Projection of Koopman embedded observable ', 'Phase trajectory using Algorithm 3', 'Phase trajectory using Algorithm 4'}...
    , 'Interpreter','latex');


figure(Lyapunov_Fig)
% subplot 211
plot(tspan(1:end - 1), V_Set_3, 'LineWidth', 3.0)
grid on
xlabel('$t/s$', 'Interpreter','latex')
ylabel('Lyapunov function $V_1, \, V_3$', 'Interpreter','latex')
hold on

% subplot 212
plot(tspan(1:end - 1), V_Set_4, 'LineWidth', 3.0);
grid on
xlabel('$t/s$', 'Interpreter','latex')
ylabel('Lyapunov function $V_1, \, V_3$', 'Interpreter','latex')
hold on


end

str1 = {};
str2 = {};
% str3 = {};
% str4 = {};
for i = 1 : times
    str1{i} = ['$x_1$ using Algorithm 3 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
    str2{i} = ['reconstructed $x_{1,rec}$ using Algorithm 3 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']'];
%     str3{i} = ['$x_1$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
%     str4{i} = ['reconstructed $x_{1,rec}$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']'];
end


figure(State_Each_Fig_1)
subplot 211
legend(str1{1},str2{1},str1{2},str2{2});

str1 = {};
str2 = {};
for i = 1 : times
    str1{i} = ['$x_1$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
    str2{i} = ['reconstructed $x_{1,rec}$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']'];
end

subplot 212
legend(str1{1},str2{1},str1{2},str2{2});





figure(State_Each_Fig_2)

subplot 211
str1 = {};
str2 = {};
for i = 1 : times
    str1{i} = ['$x_2$ using Algorithm 3 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
    str2{i} = ['reconstructed $x_{2,rec}$ using Algorithm 3 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']'];
%     str3{i} = ['$x_2$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
%     str4{i} = ['reconstructed $x_{2,rec}$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']'];
end

legend(str1{1},str2{1},str1{2},str2{2});

subplot 212
str3 = {};
str4 = {};

for i = 1 : times
    str3{i} = ['$x_2$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
    str4{i} = ['reconstructed $x_{2,rec}$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']'];
end
legend(str3{1},str4{1},str3{2},str4{2});




stru1 = {};
stru2 = {};
for i = 1 : times
    stru1{i} = ['Control input $u$ using Algorithm 3 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
    stru2{i} = ['Control input $u$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
end
figure(U_Control_Input_Fig)
% subplot 211
legend(stru1{1}, stru2{1}, stru1{2}, stru2{2})


% stru = {};
% for i = 1 : times
%     stru{i} = ['Control input $u$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
% end
% figure(U_Control_Input_Fig)
% subplot 212
% legend(stru{1}, stru{2}, stru{3}, stru{4})



figure(Lyapunov_Fig)
strV1 = {};
strV2 = {};
for i = 1 : times
    strV1{i} = ['Lyapunov function $V_1$ using Algorithm 3 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
    strV2{i} = ['Lyapunov function $V_3$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
end
% subplot 211
legend(strV1{1}, strV2{1}, strV1{2}, strV2{2})


% strV = {};
% for i = 1 : times
%     strV{i} = ['Lyapunov function $V_3$ using Algorithm 4 with intial value [' , num2str(X_init_Set(1, i)), ',', num2str(X_init_Set(2, i)) ,']' ];
% end
% subplot 212
% legend(strV{1}, strV{2}, strV{3}, strV{4})




%% Controller Design Online
function [lambda, xi, cost] = Get_Online_Control_CVX(x, A, P, UpperBound_P, InvariantF)
    
%     figure(InvariantF)
%     hold on
%     scatter(x(1, :), x(2, :))

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
        end
        xi * ones(Nd, 1) == x;
        lambda' * ones(Nd, 1) <= 1;
        lambda >= 0;
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




function [lambda, xi, cost, vi] = yalmip_online_interpolation(x, A, P, UpperBound_P)
    persistent ctrl Nd n  % 保持问题结构 & 加速

    % 初始化一次 YALMIP 优化器
    if isempty(ctrl)
        Nd = size(A, 2);
        n = size(A{1}, 1);

        xi_var = sdpvar(n, Nd, 'full');
        lambda_var = sdpvar(Nd, 1);
        x_param = sdpvar(n, 1);

        Constraints = [];
        for i = 1:Nd
            Constraints = [Constraints, ...
                [lambda_var(i), xi_var(:, i)'; xi_var(:, i), lambda_var(i) * P{i}] >= 0];
        end
        Constraints = [Constraints, ...
            xi_var * ones(Nd, 1) == x_param, ...
            lambda_var >= 0, ...
            sum(lambda_var) <= 1];

        z = [x_param; reshape(xi_var(:, 2:Nd), [], 1)];
        Objective = z' * UpperBound_P * z + lambda_var' * lambda_var;

        options = sdpsettings('solver', 'mosek', 'verbose', 0, ...
                              'debug', 0, ...
                              'cachesolvers', 1);
        ctrl = optimizer(Constraints, Objective, options, x_param, {lambda_var, xi_var, Objective});
    end

    % 在线调用
    tic
        out = ctrl(x);
        lambda = out{1};
        xi = out{2};
        cost = out{3};
    toc
%     x - xi * ones(Nd, 1)
    % 计算 vi（只对 lambda 不为 0 的情况）
    vi = zeros(size(xi));
    for i = 1:Nd
%         if lambda(i) >= 1e-6
%             vi(:, i) = xi(:, i) ./ lambda(i);
%         end
        vi(:, i) = xi(:, i) ./ lambda(i);
    end

%     % 可视化（可选）
%     scatter(x(1), x(2)); hold on;
%     scatter(vi(1, :), vi(2, :));
end













