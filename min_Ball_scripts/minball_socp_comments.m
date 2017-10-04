
function [result,info] = minball_socp_comments(p,q,w_in,gurobi_params,local_params)

% Solve
%   min R s.t. ||c - w_j|| <= R, j = 1,2,...q  
% over optimization variables R (scalar) and c (vector with nx
% elements) and count the number of samples satisfying
%   ||c - w_j|| <= R for j in {1,2,...ns}
% where ns >= q
%
% Input arguments
%  p.nx: dimension of c and w_j
%  p.bnd: max absolute value of any optimization variable (set this
%    to some arbitrarily large number e.g. 1e8)
%  q: number of samples included in constraint w_j, j = 1,...q  
%  w_in: matrix with jth row equal to w_j'
%  gurobi_params: structure of parameters to be passed to gurobi,
%    e.g. set gurobi_params.outputflag = 0 to suppress gurobi print
%    output
%  local_params.tol: tolerance value for deciding whether a
%  constraint is satisfied (e.g. 1e-8)
%
% Return arguments
%  result.x: optimal solution
%  result.objval: optimal value of objective
%  result.status: gurobi termination status (optimal/infeasible/etc)
%  result.runtime: gurobi execution time
%  info.nact: number of active SOCP constraints
%  info.qobs: number of samples defined by the first q rows of w that
%    satisfy constraints (info.qobs should always equal q)
%  info.q: number of samples specified by all rows of w that
%    satisfy constraints  

if nargin < 4
  gurobi_params.outputflag = 0;
  gurobi_params.Timelimit = 60;
end
if nargin < 5
  local_params.tol = 1e-8;
end

w = w_in(1:q,:); 
% store the first q rows of w_in in w
% e.g. w_in = normrnd(0,1,ns,p.nx)*sqrtm(w_var) + ones(ns,1)*w_mean';
% if samples are normally distributed with mean w_mean and
% covariance w_var

cones(q).index = 0;
% initialize cones as an array (of length q) of structures with a
% field named index

% X = [R; c; xv(1:p.nx*q)] 
% define X as the vector containing the optimization variables:
% scalar R (the minimum radius), vector c with nx elements (the
% centre of the hypersphere), and vectors xv_1, ... xv_q, each
% with nx elements, which will be used to specify the second order
% cone (SOC) constraints 

AA = [zeros(p.nx*q,1), repmat(eye(p.nx),q,1), eye(p.nx*q)];
bb = reshape(w',q*p.nx,1);
% AA and bb define the linear constraints
% AA and bb are defined in blocks consisting of nx rows: 
%   AA = [0 | I | I 0 ... 0   bb = [w_1
%         0 | I | 0 I ... 0         w_2
%                 ...               ...
%         0 | I | 0 0 ... I]        w_q]
% the linear constraints will be implemented as
% AA*X = bb, which (given the definition of X) is equivalent to
% xv_j = w_j - c, j = 1,2,... q

for i = 1:q
  cones(i).index = [1,1+p.nx+(i-1)*p.nx+(1:p.nx)];
end
% define the SOC constraints: cones(j).index = [k, idx] specifies
% the jth cone constraint as
%   X(k)^2 >= sum(X(idx).^2), X(k) >= 0
% here idx points to xv_j, so this is equivalent to 
%   R^2 >= ||xv_j||^2, R >= 0 
% which, given the linear constraint defining xv, is equivalent to
%   R^2 >= ||w_j - c||^2, R >= 0

model.A = sparse(AA);      % gurobi requires linear constraint
                           % matrix to be sparse (not sure why)
model.rhs = bb;            % rhs of linear constraints
model.cones = cones(1:q);  % SOC constraints
model.sense = '=';         % all linear constraints are equality
model.vtype = 'C';         % all variables are continous (i.e. no
                           % integer variables)
model.obj = [1; zeros((1+q)*p.nx,1)];
                           % objective vector, f, defined so that
                           % f'*X = R 
model.ub = p.bnd*ones(1+p.nx+p.nx*q,1);
model.lb = -model.ub;      % absolute bounds on all variables 

result = gurobi(model,gurobi_params);
                           % solve the optimization problem
R = result.x(1);
c = result.x(1+(1:p.nx));
info.nact = sum(result.qcslack < sqrt(local_params.tol));

% count how many samples in w_in satisfy the constraints 
ns = size(w_in,1);
rhs_sq = (R+local_params.tol)^2;
viol = sum(sum((w - ones(q,1)*c').^2,2) > rhs_sq);
info.qobs = q - viol;
viol = viol + sum(sum((w_in(q+1:end,:) - ones(ns-q,1)*c').^2,2) > rhs_sq);
info.q = ns - viol;
