function [result,info] = minball_socp(p,q0,w_in,gurobi_params,local_params)

% Solve
% min R s.t. ||c - wj|| <= R, j = 1,...q  
% <<<no samples discarded>>>

if nargin < 4
  gurobi_params.outputflag = 0;
  gurobi_params.Timelimit = 60;
end
if nargin < 5
  local_params.tol = 1e-8;
end

w = w_in(1:q0,:); %normrnd(0,1,q0,p.nx)*sqrtm(p.w_var) + ones(q0,1)*p.w_mean';
b = true(q0,1);
cones(q0).index = 0;

% X = [R; c; xv(1:p.nx*ns_loc)]
AA = [zeros(p.nx*q0,1), repmat(eye(p.nx),q0,1), eye(p.nx*q0)];
bb = reshape(w(b,:)',q0*p.nx,1);
for i = 1:q0
  cones(i).index = [1,1+p.nx+(i-1)*p.nx+(1:p.nx)];
end
model.A = sparse(AA);
model.rhs = bb;
model.cones = cones(1:q0);
model.sense = '=';
model.vtype = 'C';
model.obj = [1; zeros((1+q0)*p.nx,1)];
model.ub = p.bnd*ones(1+p.nx+p.nx*q0,1);
model.lb = -model.ub;

result = gurobi(model,gurobi_params);
R = result.x(1);
c = result.x(1+(1:p.nx));
info.nact = sum(result.qcslack < sqrt(local_params.tol));

% check constraints
ns = size(w_in,1);
rhs_sq = (R+local_params.tol)^2;
viol = sum(sum((w - ones(q0,1)*c').^2,2) > rhs_sq);
info.q0obs = q0 - viol;
viol = viol + sum(sum((w_in(q0+1:end,:) - ones(ns-q0,1)*c').^2,2) > rhs_sq);
info.q = ns - viol;
