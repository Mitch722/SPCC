
A = [ 1,1 ];
B = [ 0,0 ];

dist = find_distance(A,B);

c = A';
d = B';

dist2 = find_distance(c,d);

% a test to see if arrays can be handled
assert(dist == dist2, 'The two distances are not equal')

assert(dist == sqrt(2), 'The dist is not correct')