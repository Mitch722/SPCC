%% Circle to fit round data
% Put a circle of the smallest radius which touches at least two points of
% the sample.
tic

n = 100;
points =  randn(2,n);

plot(points(1,:),points(2,:),'.')

total = sum(points,2);

% COM is the centre of mass of the points 
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

%% Move the circle along the line towards P until the Distance from the centre are equal
% int_sect is the first point of intersection with the circle. 

% P point %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int_sect = points(:,index);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    % Find the next furthest point away from the new centre
    dist_sq_1 = ( new_c(1,1) - points_1(1,:) ).^2 + ( new_c(2,1) - points_1(2,:) ).^2;

    [max_dist_sq_1,index_1] = max(dist_sq_1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Q_point = points_1(:,index_1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

plot(new_c(1,1),new_c(2,1),'+')

%% Move Circle towards Bisector of Two points until it hits the third

b_len = find_distance( int_sect, Q_point );
division = ( ( int_sect(2,1) - Q_point(2,1) ) / (int_sect(1,1) - Q_point(1,1)) ); 

ang2 = atan( division );

% the X value for mid point 
eX = 0.5*b_len*cos(ang2);
eY = 0.5*b_len*sin(ang2);

if int_sect < Q_point
    eX = -eX;
    eY = -eY;
    
    fprintf('Hello, the point gets added. \n')
end

% Bisector for the two points
E_point = [ Q_point(1,1) + eX; Q_point(2,1) + eY ];

hold on
plot(E_point(1,1), E_point(2,1),'O')

%% Check the code for the bisector

bisector_line = [int_sect(1,1), Q_point(1,1) ; int_sect(2,1), Q_point(2,1) ];
hold on 
plot(bisector_line(1,:), bisector_line(2,:), 'm')

%% Move the Radius Closer

c_len = find_distance(new_c, E_point);
c_len_2 = c_len;

division = ( ( new_c(2,1) - E_point(2,1) ) / (new_c(1,1) - E_point(1,1)) );

ang3 = atan(division);

new_c_2 = new_c;

%% New Radius is closer

diff = 0.00001;
c_len_2 = c_len_2 - diff;

x_val = c_len_2*cos(ang3);
y_val = c_len_2*sin(ang3); 

% if int_sect(1,1) > COM(1,1)
%     x_val = -x_val;
%     y_val = -y_val;    
% end

new_c_2 = new_c_2 + [ x_val ; y_val ];

toc

plot(new_c_2(1,1),new_c_2(2,1),'+')
