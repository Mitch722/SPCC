
function [R, c, sub_s, x] = Algo(M,r)


% Generate the Multisample
% x.M, no. of global samples in the Multi-sample
x.M = M;
% x.no_sub_samp, the no of sets of sub-samples
x.no_sub_samp = 1;
% x.n, no. of samples in the sub-sample
x.n = x.M / x.no_sub_samp;

assert(x.n == norm(x.n), 'The number of Samples is not divisable by the number of sub-samples')

% x.dim, the dimension of the samples
x.dim = 2;
% x.sample, the global mutlisample
x.sample = randn(x.M, x.dim);
%%

sub_s = x.sample((1 : r),:);

[R, c, ~] = min_R_SOCP(sub_s);

data = x.sample;
data(1:r,:) = [];

radii = zeros(1, length(data));

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
    
    if i == x.M - r
        break 
    end
    
end

%%
arr = (sub_s - c);
val = arr.^2;
val = sum(val, 2);

[~, index] = max(val);

sub_s(index, :) = [];

[R, c, ~] = min_R_SOCP(sub_s);

end