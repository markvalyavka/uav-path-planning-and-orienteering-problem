clear all
close all

g = 9.81

refcpc = table2array(readtable('cpc_final.csv','NumHeaderLines',1));
refcpc_a = refcpc(:,15:17);
refcpc_t = refcpc(:,1);
refcpc_thrust = refcpc_a(:,:) - [0,0,g];


% ref = table2array(readtable('../samples_equidistant_08.csv','NumHeaderLines',1));
ref = table2array(readtable('../samples_pmm.csv','NumHeaderLines',1));
ref = ref(1:3:end,:)
% ref = ref(1:7:end,:)
reft=ref(:,1);
refp=ref(:,2:4);
refv=ref(:,9:11);
refa=ref(:,15:17);
probmin = optimproblem('ObjectiveSense','min');
N_points = 40
t = optimvar('t',1,1,'LowerBound',0,'UpperBound', reft(N_points,:));
p = optimvar('p',N_points,3);
v = optimvar('v',N_points,3);
a = optimvar('a',N_points,3,'LowerBound',-30,'UpperBound',30);
probmin.Objective = t; %sum(t);

probmin.Constraints.constpstart = p(1,:) == refp(1,:)
probmin.Constraints.constpend = p(N_points,:) == refp(N_points,:)
probmin.Constraints.constvstart = v(1,:) == refv(1,:)
%probmin.Constraints.constvend = v(N_points,:) == refv(N_points,:)

max_jerk = [440,440,440]

% for i:2:N_points-1
% x0.t = diff(reft(N_points,:));
x0.t = reft(N_points,:);
x0.p = refp(1:N_points,:);
x0.v = refv(1:N_points,:);
x0.a = refa(1:N_points,:);
% x0.a = zeros(N_points,3);
% end
max_diff_p = [0.1,0.1,0.1]

for i=1:N_points
    const_name_a = ['consta' num2str(i)]
    probmin.Constraints.(const_name_a) = a(i,1)^2 + a(i,2)^2 + (a(i,2)+g)^2 <= 32.94^2;
end

% for i=2:N_points-1
%     const_name_p_plus = ['constpclose_plus' num2str(i)]
%     const_name_p_minus = ['constpclose_minus' num2str(i)]
%     probmin.Constraints.(const_name_p_plus) = p(i,:)<=x0.p(i,:)+max_diff_p;
%     probmin.Constraints.(const_name_p_minus) = p(i,:)>=x0.p(i,:)-max_diff_p;
% end

for i=2:N_points
    const_name_p = ['constp' num2str(i)]
    const_name_v = ['constv' num2str(i)]
    const_name_ajplus = ['constajplus' num2str(i)]
    const_name_ajminus = ['constajminus' num2str(i)]
    t_part = t/(N_points-1);
%     t_part = t(i-1);
    probmin.Constraints.(const_name_p) = p(i,:)==p(i-1,:)+v(i-1,:)*t_part+0.5*a(i-1,:)*t_part^2;
    probmin.Constraints.(const_name_v) = v(i,:)==v(i-1,:)+a(i-1,:)*t_part;
    probmin.Constraints.(const_name_ajplus) = a(i,:)<=a(i-1,:)+max_jerk*t_part;
    probmin.Constraints.(const_name_ajminus) = a(i,:)>=a(i-1,:)-max_jerk*t_part;
end

options_front = optimoptions('fmincon','MaxFunctionEvaluations',20000,'MaxIterations',10000);
show(probmin)
sol = solve(probmin,x0,'options',options_front)
new_time = sum(sol.t)
oldtime = reft(N_points)

vecnorm(refcpc_a,2,2);

diff_t = diff(refcpc_t(1:end-1));
diff_t_stack = repmat(diff_t,1,3);
diff_a = diff(refcpc_thrust(1:end-1,:));
j = diff_a ./ diff_t_stack;

G_vec = zeros(N_points,3);
G_vec(:,3) = -g;
T = sol.a-G_vec;
vecnorm(T,2,2);
plot3(sol.p(:,1),sol.p(:,2),sol.p(:,3),'g')
hold on
plot3(refp(1:N_points,1),refp(1:N_points,2),refp(1:N_points,3),'r')
