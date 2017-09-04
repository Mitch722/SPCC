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



%% Find the next furthest point from the same center

points_1 = points;
points_1(:,index) = [];

dist_sq_1 = ( COM(1,1) - points_1(1,:) ).^2 + ( COM(2,1) - points_1(2,:) ).^2;

[max_dist_sq_1,index_1] = max(dist_sq_1);

% max distance to the next largest point from the COM
dist_1 = sqrt(max_dist_sq_1);

Q_point = points_1(:,index_1); 

%% Move the circle along the line towards P until the Distance from the centre are equal
int_sect = points(:,index);

division = ( ( int_sect(2,1) - COM(2,1) ) / (int_sect(1,1) - COM(1,1)) ); 

ang = atan( division );

new_dist = dist;

for i = 1:100000
    diff = 0.0001;
    new_dist = new_dist - diff;

    x_val = new_dist*cos(ang);
    y_val = new_dist*sin(ang); 

    if int_sect(1,1) > COM(1,1)
        x_val = -x_val;
        y_val = -y_val;    
    end

    new_c = [int_sect(1,1) + x_val ; int_sect(2,1) + y_val ];
    
    % Perform the test to see if the point is near the other one
    Q_to_new_c = ( new_c(1,1) - Q_point(1,1) ).^2 + ( new_c(2,1) - Q_point(2,1) ).^2;
    Q_to_new_c = sqrt(Q_to_new_c);
    P_to_new_c = ( new_c(1,1) - int_sect(1,1) ).^2 + ( new_c(2,1) - int_sect(2,1) ).^2;
    P_to_new_c = sqrt(P_to_new_c);
    
    if P_to_new_c <= Q_to_new_c*1.001
        break
    end
end

plot(P_to_new_c*cos(t) + new_c(1,1), P_to_new_c*sin(t) + new_c(2,1),'m')
plot(Q_to_new_c*cos(t) + new_c(1,1), Q_to_new_c*sin(t) + new_c(2,1),'c')


