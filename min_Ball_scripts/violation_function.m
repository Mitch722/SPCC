function [ no_violate, output_points, viol_fact ] = violation_function( R, c, global_lite )

% This function finds the points in the global sample which violate the constraints generated from the sub-sample.
%
%    input Args: 
%            1. R the radius calculated from Gurobi interior point method
%
%            2. c, the centre of the hypersphere genereated by Gurobi SOCP 
%            
%            3. sub_samp, the sub-sample of points used to generate the model
%
%            4. global, the global distibution from which the sub-sample was generated
%
%
%    returns:
% 
%            1. no_violate, the number of points which violate the sub-sample model
%
%            2. points, the points which violate the sub-sample model
%
%            3. viol_fact, epsilon the proportion of points that have violated the model.

[rows_glob, cols_glob] = size(global_lite);


% global points minus fromt the centre
centre_points = global_lite - c;
% square all the individual points
centre_points_sq = centre_points.^2;

% Sum the squares of the points along the rows
sum_vects = sum(centre_points_sq, 2);
% R: the radius of the hypersphere squared
R_sqd = R^2;

% a matrix of logics showing where the points which are larger than R
% squared are
logics = sum_vects > R_sqd;
% log_mat = logics;

one_mat = ones( length( logics ), 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_violate = logics' * one_mat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

logic_n = repmat(logics, 1, cols_glob);

points = global_lite;
points(~logic_n) = 0;

pont = points(:);
pont(pont == 0) = [];

output_points = reshape(pont, [], cols_glob);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

viol_fact = no_violate / rows_glob;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
