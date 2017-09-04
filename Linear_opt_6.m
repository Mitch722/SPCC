%% Circle to fit round data
% Put a circle of the smallest radius which touches at least two points of
% the sample.

n = 100;
points =  randn(2,n);

plot(points(1,:),points(2,:),'.')

total = sum(points,2);
COM = total ./ length(points);

dist_sq = ( COM(1,1) - points(1,:) ).^2 + ( COM(2,1) - points(2,:) ).^2;

[max_dist_sq,index] = max(dist_sq);

dist = sqrt(max_dist_sq);

%% draw a Circle

t = linspace(0,2*pi);
hold on

axis equal
plot(dist*cos(t) + COM(1,1), dist*sin(t) + COM(2,1),'r')

%% Move the Circle towards the intersection point

int_sect = points(:,index);

division = ( ( int_sect(2,1) - COM(2,1) ) / (int_sect(1,1) - COM(1,1)) ); 

ang = atan( division );

factor = 1-0.1;
new_dist = dist*factor;

x_val = new_dist*cos(ang);
y_val = new_dist*sin(ang); 

if int_sect(1,1) > COM(1,1)
    x_val = -x_val;
    y_val = -y_val;    
end

new_c = [int_sect(1,1) + x_val ; int_sect(2,1) + y_val ];

plot(new_dist*cos(t) + new_c(1,1), new_dist*sin(t) + new_c(2,1), 'b' )

%% Find the next furthest point from the same center

points_1 = points;
points_1(:,index) = [];

dist_sq_1 = ( COM(1,1) - points_1(1,:) ).^2 + ( COM(2,1) - points_1(2,:) ).^2;

[max_dist_sq_1,index] = max(dist_sq_1);

% max distance to the next largest point from the COM
dist_1 = sqrt(max_dist_sq_1);





