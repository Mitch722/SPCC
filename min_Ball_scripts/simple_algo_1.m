

try load('simple_alg_1.mat')
    
catch  
    %% Generate the Multisample
    % x.M, no. of global samples in the Multi-sample
    x.M = 600;
    % x.no_sub_samp, the no of sets of sub-samples
    x.no_sub_samp = 1;
    % x.n, no. of samples in the sub-sample
    x.n = x.M / x.no_sub_samp;

    assert(x.n == norm(x.n), 'The number of Samples is not divisable by the number of sub-samples')

    % x.dim, the dimension of the samples
    x.dim = 2;
    % x.sample, the global mutlisample
    x.sample = randn(x.M, x.dim);
    
    save('simple_alg_1.mat', '-v7.3')
    
end

% r = x.M - 100;
r = 463;

sub_s = x.sample((1 : r),:);

[R, c, ~] = min_R_SOCP(sub_s);

data = x.sample;
data(1:r,:) = [];

radii = zeros(1, length(data));

hold off

for i = 1 : length(data)
    fprintf('%d\n',i)
    
    radii(1,i) = R;
    
    poin = data(i,:);
    
    if norm(poin - c) < R
        
        arr = (sub_s - c);
        val = arr.^2;
        val = sum(val, 2);
        
        val = val - R^2;
        val = abs(val);
        
        [~, index] = min(val);
        
        sub_s(index,:) = poin;
      
        [R, c, ~] = min_R_SOCP(sub_s);
        
    else 
        fprintf('Outside Hypersphere \n')
    end
    
    t = linspace(0,2*pi);
    
    
    plot(R*cos(t) + c(1,1), R*sin(t) + c(1,2))
    hold on
    plot(sub_s(:,1),sub_s(:,2),'.')

    plot(c(1,1),c(1,2),'+')
    
    plot(poin(1), poin(2), 'x')
    grid on
    axis equal
    
    if i == x.M - r
        break 
    end
    
end

difs = diff(radii);

plot(R*cos(t) + c(1,1), R*sin(t) + c(1,2),'b')
hold on
plot(sub_s(:,1),sub_s(:,2),'.')

plot(c(1,1),c(1,2),'X')

plot(poin(1), poin(2), 'x')
grid on
axis equal

%% new figure
plot(x.sample(:,1), x.sample(:,2),'+')
hold on
plot(R*cos(t) + c(1,1), R*sin(t) + c(1,2), 'k')

plot(sub_s(:,1),sub_s(:,2),'.')

plot(c(1,1),c(1,2),'X')

plot(poin(1), poin(2), 'x')


grid on
axis equal

figure 

plot(x.sample(:,1), x.sample(:,2),'.')
hold on
plot(R*cos(t) + c(1,1), R*sin(t) + c(1,2), 'b')

plot(sub_s(:,1),sub_s(:,2),'.')

plot(c(1,1),c(1,2),'X')

plot(poin(1), poin(2), 'x')


grid on
axis equal

hold off
%%
arr = (sub_s - c);
val = arr.^2;
val = sum(val, 2);

rad_samp = max(val);

rad_samp = sqrt(rad_samp);

%
% [R, c, no_viol] = min_R_SOCP(sub_s);
% 
figure 

plot(x.sample(:,1), x.sample(:,2),'.')
hold on
plot(rad_samp*cos(t) + c(1,1), rad_samp*sin(t) + c(1,2),'b')

plot(sub_s(:,1),sub_s(:,2),'.')

plot(c(1,1),c(1,2),'X')

plot(poin(1), poin(2), 'x')


grid on
axis equal

hold off
