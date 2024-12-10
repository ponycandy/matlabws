clear

%% 创建空的NLP变量
n_collocation = 5;
type = 'LGL';
obj = ps_optimizer(n_collocation,type);

s_start = 100;
s_final = 110;
n_block = 10;
obj = obj.set_NLP_horizon(s_start,s_final,n_block);
% obj = obj.update_horizon(90,110);

x_0 = [0;0];
u_0 = 0;
obj = obj.gen_empty_NLP(x_0,u_0);

x_cons = [-2,2;-0.785,0.785];
u_cons = [-0.1,0.1];
du_cons = inf*[-1,1];
obj = obj.set_constraint(x_cons,u_cons,du_cons);

%% 假设收到了一个路径
map = figure(1);

% 取数据
load("data/m_track.mat");

% 坐标系变换
Track = xy2sc(X,Y);

% 画出来
clf
plot(Track.x,Track.y);
axis equal

%% 建立NLP问题
x_init = [0;0];
obj = obj.fill_NLP(x_init,@(X,U,D_U,s)s_curv_model(X,U,D_U,s,Track));

%% 增加障碍物
[x,y] = sce2xy(Track,[103.5;105;106.22],[0.2;0.3;-0.4]);
figure(map)
axis([87,96,-55,-47])
hold on
plot([x;x(1)],[y;y(1)],'Marker','none','Color',[1,0,0],'LineStyle','-',LineWidth=2);

%% 一次障碍物处理
% 将点映射到曲率坐标系上
[s,e] = xy2sce(Track,x,y);

% 判断该往哪边跑，这里的逻辑是按照几何中心确定方向
e_O = 0;
for i = 1:length(s)
    e_O = e_O + e(i);
end
e_O = e_O/length(s);
e_O = e_O - x_init(1);

% 检查边缘点是否在路径上
while s(1) <= s_start || s(end) >= s_final
    if s(1) <= s_start
        s = s(2:end);
        e = e(2:end);
    end

    if s(end) >= s_final
        s = s(1:end - 1);
        e = e(1:end - 1);
    end
end

if ~isempty(s)
    n_s = length(s);

    X = cell(n_s,1);
    for i = 1:n_s
        X{i} = obj.get_X_t(s(i));
    end
    
    g_add = cell(n_s,1);
    
    for i = 1:n_s
        g_add{i} = X{i}(1);
    end
    
    obj.g = {obj.g;g_add};
    obj.g = vertcat(obj.g{:});
    obj.cellPointer_g = obj.cellPointer_g + n_s;
    
    if e_O <= 0
        obj.ub_g = [obj.ub_g;inf*ones(n_s,1)];
        obj.lb_g = [obj.lb_g;e];
    else
        obj.ub_g = [obj.ub_g;e];
        obj.lb_g = [obj.lb_g;-inf*ones(n_s,1)];
    end
    obj.vectorPointer_g = obj.vectorPointer_g + n_s;
end

%% 求解
obj = obj.solve_NLP();

%% 查看求解结果
obj = obj.show_solution(2);

X_opt = obj.X_opt;
X_opt = reshape(X_opt,[2,n_collocation*n_block]);
e = X_opt(1,:)';
erro_PHI = X_opt(2,:)';
s = obj.grid;

[X_new,Y_new] = sce2xy(Track,s,e);

figure(map)
hold on
plot(X_new,Y_new,'Marker','none','Color',[0,0,1],'LineStyle','-','LineWidth',2);
plot(X_new(1:n_collocation:end),Y_new(1:n_collocation:end),'LineStyle','none','Marker','o','Color',[0 0 1],'LineWidth',2)