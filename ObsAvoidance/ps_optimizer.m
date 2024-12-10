classdef ps_optimizer
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %   lengthBuffer is (real length - 1)

    % 可以求解非线性规划问题（NLP）
    % 问题维度不限，但不建议求解大型复杂非线性问题（因为不好用）
    % 代价函数支持状态量、控制量及控制量导数为变量，且能够接受时变参数
    % 约束支持状态量、控制量及控制量导数的边界约束，且开放g函数自由定义其他约束

    properties
        %% collocation parameters
        n_collocation
        method
        root
        D_matrix
        weight
        L
        L_start
        L_final

        %% NLP parameters
        t_s
        t_f
        n_block
        T_s
        T_f
        h
        x_init
        x_0
        u_0
        n_x
        n_u
        ub_x
        lb_x
        ub_u
        lb_u
        ub_du
        lb_du
        w
        cellPointer_w
        w0
        lb_w
        ub_w
        vectorPointer_w
        g
        cellPointer_g
        lb_g
        ub_g
        vectorPointer_g
        J
        w_opt
        X_opt
        U_opt
        grid
        vectorPointer_grid
        
        %% function flags
        flag_set_NLP_horizon = 0;
        flag_gen_empty_NLP = 0;
        flag_gen_grid = 0;
        flag_set_constraint = 0;
        flag_fill_NLP = 0;
        flag_solve_NLP = 0;
        flag_extract_solution = 0;

    end

    methods
        %% initial function
        function obj = ps_optimizer(n_collocation,method)
            % INPUTS
            % n_collocation: the length of legendre polynomial's root
            % method: the type of legendre-gauss integral point
            
            obj.n_collocation = n_collocation;
            obj.method = method;

            obj = legendre_lagrange_param_gen(obj,n_collocation - 1,method);

            obj.flag_set_NLP_horizon = 0;
            obj.flag_gen_empty_NLP = 0;
            obj.flag_set_constraint = 0;
            obj.flag_fill_NLP = 0;
            obj.flag_solve_NLP = 0;
            obj.flag_extract_solution = 0;
            obj.flag_gen_grid = 0;
        end

        %% basic function
        function obj = legendre_lagrange_param_gen(obj,N,type)
            % INPUTS
            % N: the order of legendre polynomial. (the length of root is N + 1)
            % type: the type of legendre-gauss integral point
            
            P_N_plus1 = legendre_coefficient(N + 1);
            
            P_N = legendre_coefficient(N);
            
            P_N_plus1_dot = polyder(P_N_plus1);
            
            P_N_dot = polyder(P_N);
            
            obj.L = zeros(N + 1,N + 1);
            
            switch type
                case 'LG'
            
                    obj.root = sort(roots(P_N_plus1));
            
                    for k = 1:N + 1
                        obj.L(k,:) = 1/(polyval(P_N_plus1_dot,obj.root(k)))*deconv(P_N_plus1,[1 -obj.root(k)]);
                    end
            
                case 'LGR'
            
                    obj.root = sort(roots(poly_add(P_N,P_N_plus1)));
            
                    for k = 1:N + 1
                        obj.L(k,:) = (1 - obj.root(k))/(2*(N + 1)*polyval(P_N,obj.root(k)))*deconv(poly_add(P_N,P_N_plus1),[1 -obj.root(k)]);
                    end
            
                case 'LGL'
            
                    obj.root = sort(roots(conv(P_N_dot,[-1 0 1])));
                    
                    for k = 1:N + 1
                        obj.L(k,:) = 1/(N*(N + 1)*polyval(P_N,obj.root(k)))*deconv(conv(P_N_dot,[1 0 -1]),[1 -obj.root(k)]);
                    end
            
                otherwise
            
            end
            
            obj.L_start = zeros(N + 1,0);
            obj.L_final = zeros(N + 1,0);
            obj.D_matrix = zeros(N + 1,N + 1);
            obj.weight = zeros(N + 1,1);
            for k = 1:N + 1
                coeff = obj.L(k,:);
            
                obj.L_start(k) = polyval(coeff,-1);
                obj.L_final(k) = polyval(coeff,1);
            
                pder = polyder(coeff);
                for r = 1:N + 1
                    obj.D_matrix(r,k) = polyval(pder,obj.root(r));
                end
            
                pint = polyint(coeff);
                obj.weight(k) = polyval(pint,1.0) - polyval(pint,-1.0);
            end
        end

        function obj = set_NLP_horizon(obj,t_s,t_f,n_block)
            % INPUTS
            % t_s: strat time of NLP
            % t_f: final time of NLP
            % n_block: number of blocks in [t_s,t_f]
            
            D_t = (t_f - t_s)/n_block;

            obj.t_s = t_s;
            obj.t_f = t_f;
            obj.n_block = n_block;

            obj.T_s = zeros(n_block,1);
            obj.T_f = zeros(n_block,1);
            
            t = t_s;
            for i = 1:n_block
                obj.T_s(i) = t;
                t = t + D_t;
                obj.T_f(i) = t;
            end
            obj.h = (t_f - t_s)/2/n_block;

            obj.flag_set_NLP_horizon = 1;
            obj.flag_gen_empty_NLP = 0;
            obj.flag_set_constraint = 0;
            obj.flag_fill_NLP = 0;
            obj.flag_solve_NLP = 0;

            obj = obj.gen_grid();
        end

        function obj = gen_grid(obj)
            if obj.flag_set_NLP_horizon == 0
                fprintf('horizon has not been set \n');
                return
            end

            obj.grid = zeros(obj.n_collocation*obj.n_block,1);
            obj.vectorPointer_grid = 1;
            for i = 1:obj.n_block
                obj.grid(obj.vectorPointer_grid:obj.vectorPointer_grid + obj.n_collocation - 1) = ...
                    obj.root*(obj.T_f(i) - obj.T_s(i))/2 + (obj.T_f(i) + obj.T_s(i))/2;
                obj.vectorPointer_grid = obj.vectorPointer_grid + obj.n_collocation;
            end

            obj.flag_gen_grid = 1;
        end

        function obj = gen_empty_NLP(obj,x_0,u_0)
            % INPUTS
            % x_0: warm start value of x
            % u_0: warm start value of u

            if obj.flag_set_NLP_horizon == 0
                fprintf('horizon has not been set \n');
                return
            end

            obj.x_0 = x_0;
            obj.u_0 = u_0;
            
            % length of variables
            obj.n_x = length(x_0);
            obj.n_u = length(u_0);
            
            % Start with an empty NLP
            obj.w = cell(obj.n_collocation*obj.n_block*2,1);
            obj.cellPointer_w = 1;
            obj.w0 = zeros((obj.n_x + obj.n_u)*obj.n_collocation*obj.n_block,1);
            obj.lb_w = zeros((obj.n_x + obj.n_u)*obj.n_collocation*obj.n_block,1);
            obj.ub_w = zeros((obj.n_x + obj.n_u)*obj.n_collocation*obj.n_block,1);
            obj.vectorPointer_w = 1;
            
            obj.g = cell(obj.n_block*(obj.n_collocation + 2) - 1,1);
            obj.cellPointer_g = 1;
            obj.lb_g = zeros(obj.n_block*obj.n_x + (obj.n_block - 1)*obj.n_u + obj.n_collocation*obj.n_block*obj.n_x,1);
            obj.ub_g = zeros(obj.n_block*obj.n_x + (obj.n_block - 1)*obj.n_u + obj.n_collocation*obj.n_block*obj.n_x,1);
            obj.vectorPointer_g = 1;
            
            obj.J = 0;

            obj.flag_gen_empty_NLP = 1;
        end

        function obj = set_constraint(obj,x_cons,u_cons,du_cons)
            % INPUTS
            % x_cons: [lb_x,ub_x]
            % u_cons: [lb_u,ub_u]
            % du_cons: [lb_du,ub_du]
            
            obj.lb_x = x_cons(:,1);
            obj.ub_x = x_cons(:,2);
            obj.lb_u = u_cons(:,1);
            obj.ub_u = u_cons(:,2);
            obj.lb_du = du_cons(:,1);
            obj.ub_du = du_cons(:,2);

            obj.flag_set_constraint = 1;
            obj.flag_fill_NLP = 0;
            obj.flag_solve_NLP = 0;
        end

        function obj = fill_NLP(obj,x_init,fun_derivative)
            % INPUTS
            % x_init: initial value of x
            % fun_derivative: functions like @(X,U,d_U,t)fun
            % , and return [derivative of X, derivative of J]

            if obj.flag_gen_empty_NLP == 0
                fprintf('NLP has not been created \n')
                return
            end

            if obj.flag_set_constraint == 0
                fprintf('constraint has not been set \n')
                return
            end

            if obj.flag_gen_grid == 0
                fprintf('grid has not been created \n')
                return
            end

            % init state value
            obj.x_init = x_init;

            % "Lift" initial conditions
            Xk_final_last = x_init;

            % reset vectorPointer_grid
            obj.vectorPointer_grid = 1;
            
            % fill the NLP
            for k = 1:obj.n_block
                % creat NLP control variables
                Ukj = cell(obj.n_collocation,1);
                for j = 1:obj.n_collocation
                    Ukj{j} = casadi.MX.sym(['U_' num2str(k) '_' num2str(j)],obj.n_u,1);
                
                    obj.w{obj.cellPointer_w} = Ukj{j};
                    obj.cellPointer_w = obj.cellPointer_w + 1;
            
                    obj.lb_w(obj.vectorPointer_w:obj.vectorPointer_w + obj.n_u - 1) = obj.lb_u;
                    obj.ub_w(obj.vectorPointer_w:obj.vectorPointer_w + obj.n_u - 1) = obj.ub_u;
                    obj.w0(obj.vectorPointer_w:obj.vectorPointer_w + obj.n_u - 1) = obj.u_0;
                    obj.vectorPointer_w = obj.vectorPointer_w + obj.n_u;
                end
                
                % crear NLP state variables
                Xkj = cell(obj.n_collocation,1);
                for j = 1:obj.n_collocation
                    Xkj{j} = casadi.MX.sym(['X_' num2str(k) '_' num2str(j)],obj.n_x,1);
                
                    obj.w{obj.cellPointer_w} = Xkj{j};
                    obj.cellPointer_w = obj.cellPointer_w + 1;
            
                    obj.lb_w(obj.vectorPointer_w:obj.vectorPointer_w + obj.n_x - 1) = obj.lb_x;
                    obj.ub_w(obj.vectorPointer_w:obj.vectorPointer_w + obj.n_x - 1) = obj.ub_x;
                    obj.w0(obj.vectorPointer_w:obj.vectorPointer_w + obj.n_x - 1) = obj.x_0;
                    obj.vectorPointer_w = obj.vectorPointer_w + obj.n_x;
                end
                
                % Loop over collocation points
                Xk_start = obj.L_start(1)*Xkj{1};
                Xk_final = obj.L_final(1)*Xkj{1};
                Uk_start = obj.L_start(1)*Ukj{1};
                Uk_final = obj.L_final(1)*Ukj{1};
                for j = 1:obj.n_collocation
                    % Expression for derivative at the collocation point
                    xp = obj.D_matrix(j,1)*Xkj{1};
                    up = obj.D_matrix(j,1)*Ukj{1};
                    for r = 2:obj.n_collocation
                       xp = xp + obj.D_matrix(j,r)*Xkj{r};
                       up = up + obj.D_matrix(j,r)*Ukj{r};
                    end
                    
                    % Append collocation equations
                    [df_j,q_j] = fun_derivative(Xkj{j},Ukj{j},up/obj.h,obj.grid(obj.vectorPointer_grid));
                    obj.vectorPointer_grid = obj.vectorPointer_grid + 1;
            
                    obj.g{obj.cellPointer_g} = obj.h*df_j - xp;
                    obj.cellPointer_g = obj.cellPointer_g + 1;
            
                    obj.lb_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_x - 1) = zeros(obj.n_x,1);
                    obj.ub_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_x - 1) = zeros(obj.n_x,1);
                    obj.vectorPointer_g = obj.vectorPointer_g + obj.n_x;

                    obj.g{obj.cellPointer_g} = up;
                    obj.cellPointer_g = obj.cellPointer_g + 1;

                    obj.lb_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_u - 1) = obj.h*obj.lb_du;
                    obj.ub_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_u - 1) = obj.h*obj.ub_du;
                    obj.vectorPointer_g = obj.vectorPointer_g + obj.n_u;
                    
                    % Add contribution to start and final state
                    Xk_start = Xk_start + obj.L_start(j)*Xkj{j};
                    Xk_final = Xk_final + obj.L_final(j)*Xkj{j};
                    Uk_start = Uk_start + obj.L_start(j)*Ukj{j};
                    Uk_final = Uk_final + obj.L_final(j)*Ukj{j};
                    
                    % Add contribution to quadrature function
                    obj.J = obj.J + obj.weight(j)*q_j*obj.h;
                end
                Xk_start = Xk_start - obj.L_start(1)*Xkj{1};
                Xk_final = Xk_final - obj.L_final(1)*Xkj{1};
                Uk_start = Uk_start - obj.L_start(1)*Ukj{1};
                Uk_final = Uk_final - obj.L_final(1)*Ukj{1};
                
                % Add equality constraint
                if k ~= 1
                    obj.g{obj.cellPointer_g} = Uk_start-Uk_final_last;
                    obj.cellPointer_g = obj.cellPointer_g + 1;
                    
                    obj.lb_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_u - 1) = zeros(obj.n_u,1);
                    obj.ub_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_u - 1) = zeros(obj.n_u,1);
                    obj.vectorPointer_g = obj.vectorPointer_g + obj.n_u;
                end
                
                obj.g{obj.cellPointer_g} = Xk_start-Xk_final_last;
                obj.cellPointer_g = obj.cellPointer_g + 1;
            
                obj.lb_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_x - 1) = zeros(obj.n_x,1);
                obj.ub_g(obj.vectorPointer_g:obj.vectorPointer_g + obj.n_x - 1) = zeros(obj.n_x,1);
                obj.vectorPointer_g = obj.vectorPointer_g + obj.n_x;
                
                % update equality constraint variables
                Xk_final_last = Xk_final;
                Uk_final_last = Uk_final;
            end

            obj.flag_fill_NLP = 1;
            obj.flag_solve_NLP = 0;
        end

        function obj = solve_NLP(obj)
            if obj.flag_fill_NLP == 0
                fprintf('NLP has somthing worng \n')
                return
            end

            % create an NLP solver
            prob = struct('f', obj.J, 'x', vertcat(obj.w{:}), 'g', vertcat(obj.g{:}));

            opts = struct();
            opts.ipopt.max_iter = 1e3;
            opts.ipopt.print_level = 0;
            opts.print_time = 1;
            opts.ipopt.acceptable_tol = 1e-8;
            opts.ipopt.acceptable_obj_change_tol = 1e-6;

            solver = casadi.nlpsol('solver','ipopt',prob,opts);
            
            % solve the NLP
            sol = solver('x0',obj.w0,'lbx',obj.lb_w,'ubx',obj.ub_w,...
                        'lbg',obj.lb_g,'ubg',obj.ub_g);
            obj.w_opt = full(sol.x);

            obj.flag_solve_NLP = 1;
        end

        %% auxiliary function
        function obj = reset_param(obj)
            % reset NLP parameters
            obj.x_init = [];
            obj.cellPointer_w = 1;
            obj.vectorPointer_w = 1;
            obj.vectorPointer_g = 1;
            obj.J = [];
            obj.w_opt = [];
            obj.X_opt = [];
            obj.U_opt = [];
            obj.grid = [];
            obj.vectorPointer_grid = 1;
        end

        function obj = update_horizon(obj,t_s,t_f)
            % INPUTS
            % t_s: strat time of NLP
            % t_f: final time of NLP
            % n_block: number of blocks in [t_s,t_f

            if obj.flag_set_NLP_horizon == 0
                return
            end

            D_t = (t_f - t_s)/obj.n_block;

            obj.t_s = t_s;
            obj.t_f = t_f;

            obj.T_s = zeros(obj.n_block,1);
            obj.T_f = zeros(obj.n_block,1);
            
            t = t_s;
            for i = 1:obj.n_block
                obj.T_s(i) = t;
                t = t + D_t;
                obj.T_f(i) = t;
            end
            obj.h = (t_f - t_s)/2/obj.n_block;

            obj = obj.reset_param();

            % reset function flags
            obj.flag_set_NLP_horizon = 1;
            obj.flag_gen_empty_NLP = 0;
            obj.flag_gen_grid = 0;
            obj.flag_set_constraint = 0;
            obj.flag_fill_NLP = 0;
            obj.flag_solve_NLP = 0;
            obj.flag_extract_solution = 0;

            obj = obj.gen_grid();
        end

        function obj = extract_solution(obj)
            if obj.flag_solve_NLP == 0
                fprintf('NLP has not been solved \n')
                return
            end

            obj.X_opt = zeros(obj.n_collocation*obj.n_block*obj.n_x,1);
            x_lengthBuffer = 1;
            obj.U_opt = zeros(obj.n_collocation*obj.n_block*obj.n_u,1);
            u_lengthBuffer = 1;
            
            for i = 1:obj.n_block
                obj.U_opt((u_lengthBuffer - 1)*obj.n_collocation + 1:(u_lengthBuffer + obj.n_u - 1)*obj.n_collocation) = ...
                    obj.w_opt((i - 1)*(obj.n_x + obj.n_u)*obj.n_collocation ...
                    + 1:(i - 1)*(obj.n_x + obj.n_u)*obj.n_collocation + obj.n_u*obj.n_collocation);
                u_lengthBuffer = u_lengthBuffer + obj.n_u;

                obj.X_opt((x_lengthBuffer - 1)*obj.n_collocation + 1:(x_lengthBuffer + obj.n_x - 1)*obj.n_collocation) = ...
                    obj.w_opt((i - 1)*(obj.n_x + obj.n_u)*obj.n_collocation + obj.n_u*obj.n_collocation ...
                    + 1:(i - 1)*(obj.n_x + obj.n_u)*obj.n_collocation + obj.n_u*obj.n_collocation + obj.n_x*obj.n_collocation);
                x_lengthBuffer = x_lengthBuffer + obj.n_x;
            end

            obj.flag_extract_solution = 1;
        end

        function obj = show_solution(obj,arg1)
            if nargin == 1
                arg1 = 1;
            end

            if obj.flag_extract_solution == 0
                obj = obj.extract_solution();
            end

            figure(arg1);
            clf;
            subplot(2,1,1);
            hold on

            le_txt = {};

            for i = 1:obj.n_x
                plot(obj.grid,obj.X_opt(i:obj.n_x:end),'Marker','o','MarkerSize',4,'LineWidth',1.5)
                le_txt = [le_txt(:)',{['x_',num2str(i)]}];
            end

            for i = 1:obj.n_u
                plot(obj.grid,obj.U_opt(i:obj.n_u:end),'Marker','o','MarkerSize',.1,'LineWidth',1.5)
%                 stairs(obj.grid,obj.U_opt(i:obj.n_u:end),'Marker','o','MarkerSize',.1,'LineWidth',1.5)
                le_txt = [le_txt(:)',{['u_',num2str(i)]}];
            end

            xlabel('t')
            legend(le_txt,"AutoUpdate","off",'Location','best');
            legend('boxoff')
            title('优化结果')

            if obj.flag_gen_grid == 0
                fprintf('grid has not been set \n');
                return
            end

            subplot(2,1,2);
            hold on

            plot(obj.grid,ones(size(obj.grid)),'Marker','|','MarkerSize',4,'LineStyle','none','LineWidth',1)
            title('采样点分布')
        end

        function [L_tau,pointer_block]= get_L_tau(obj,t)
            % 先定位在哪一个block里
            pointer_block = ceil(obj.n_block*(t - obj.t_s)/(obj.t_f - obj.t_s));
            
            % 将位置转化至[-1,1]的相对位置
            tau = 2*(t - obj.T_s(pointer_block))/(obj.T_f(pointer_block) - obj.T_s(pointer_block)) - 1;

            % 求L_tau
            L_tau = zeros(obj.n_collocation,0);
            for i = 1:obj.n_collocation
                coeff = obj.L(i,:);
            
                L_tau(i) = polyval(coeff,tau);
            end
        end

        function U_t = get_U_t(obj,t)
            if obj.flag_fill_NLP == 0
                fprintf('NLP has not been filled \n')
                return
            end

            [L_tau,pointer_block] = obj.get_L_tau(t);

            % 将对应的U_k进行加权得到U_t
            U_t = obj.w{(pointer_block - 1)*obj.n_collocation*2 + 1}*L_tau(1);
            for i = (pointer_block - 1)*obj.n_collocation*2 + 2:...
                    (pointer_block - 1)*obj.n_collocation*2 + obj.n_collocation
                U_t = U_t + obj.w{i}*L_tau(i - (pointer_block - 1)*obj.n_collocation*2);
            end
        end

        function X_t = get_X_t(obj,t)
            if obj.flag_fill_NLP == 0
                fprintf('NLP has not been filled \n')
                return
            end

            [L_tau,pointer_block] = obj.get_L_tau(t);
            
            % 将对应的U_k进行加权得到U_t
            X_t = obj.w{(pointer_block - 1)*obj.n_collocation*2 + obj.n_collocation + 1}*L_tau(1);
            for i = (pointer_block - 1)*obj.n_collocation*2 + obj.n_collocation + 2:...
                    (pointer_block - 1)*obj.n_collocation*2 + obj.n_collocation + obj.n_collocation
                X_t = X_t + obj.w{i}*L_tau(i - (pointer_block - 1)*obj.n_collocation*2 - obj.n_collocation);
            end
        end

        function U = get_U_full(obj)
            if obj.flag_fill_NLP == 0
                fprintf('NLP has not been filled \n')
                return
            end
            
            U = cell(obj.n_collocation*obj.n_block,1);
            cellPointer_w_tmp = 1;
            cellPointer_U = 1;

            for i = 1:obj.n_block
                for j = 1:obj.n_collocation
                    U{cellPointer_U} = obj.w{cellPointer_w_tmp};
                    cellPointer_U = cellPointer_U + 1;
                    cellPointer_w_tmp = cellPointer_w_tmp + 1;
                end

                cellPointer_w_tmp = cellPointer_w_tmp + obj.n_collocation;
            end
        end

        function X = get_X_full(obj)
            if obj.flag_fill_NLP == 0
                fprintf('NLP has not been filled \n')
                return
            end

            X = cell(obj.n_collocation*obj.n_block,1);
            cellPointer_w_tmp = obj.n_collocation + 1;
            cellPointer_X = 1;

            for i = 1:obj.n_block
                for j = 1:obj.n_collocation
                    X{cellPointer_X} = obj.w{cellPointer_w_tmp};
                    cellPointer_X = cellPointer_X + 1;
                    cellPointer_w_tmp = cellPointer_w_tmp + 1;
                end

                cellPointer_w_tmp = cellPointer_w_tmp + obj.n_collocation;
            end
        end
    end
end