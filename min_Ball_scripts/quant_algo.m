
M = 600;
r = 463;

runs = 10;
no_outside = zeros(1, runs);

for i = 1:runs

    [R, c, sub_s, x] = Algo(M,r);

    [ no_violate, output_points, viol_fact ] = violation_function( R, c, x.sample, x);

%     plot(x.sample(:,1), x.sample(:,2), '.')
%     hold on
%     plot(sub_s(:,1), sub_s(:,2),'.')

    no_outside(i) = no_violate;

end
%%
save('quant_algo_1.mat', '-v7.3')

ave = sum(no_outside) / length(no_outside);

