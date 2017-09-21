% Description: works out the distance between two points
% Args: A = point a described as a 2 x 1 array
%       B = point a described as a 2 x 1 array
%
% Returns: dist, which is the distance between the two points A and B

function dist = find_distance(A,B)

difference = A - B;

dist = norm(difference);

