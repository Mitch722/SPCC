%% Blurb
% Build an FIR model coefficients using least squares.
% yk = D a + e;
% output know is the sum of the prior inputs plus a disturbance
%
%% The function
% n_c: number of coefficients
% m_c: length of D
% Yk: output sequence
% Uk: input sequence
% 
%
function F = buildFIR(Yk, Uk, n_c, m_c)

assert(m_c >= n_c, 'the number of rows should be greater than columns, necessary that m > n')
assert(length(Yk)>= n_c + m_c, 'y vector is too small: make n and m smaller'  )

[no_outs, ~] = size(Yk);    % shape recorded output
% [no_ins, lenU] = size(Uk);     % shape recorded input

Yk = Yk(:, end-m_c+1:end);        % find the end point @ k-m to k
Yk = fliplr(Yk);                % flip Yk data: the newest on the left of 
Y_data = reshape(Yk, [], 1);     % column vector n_c 

Uk = fliplr(Uk);
Uk_stacked = repmat(Uk, no_outs, 1);

Data = zeros(no_outs*m_c, n_c);
Data(1: no_outs, :) = Uk_stacked(:, 1: n_c); 


for i = 1: m_c -1
   
    Data(i*no_outs+1 : (i+1)*no_outs, :) = Uk_stacked(:, i+1: i+n_c);
    
end

psudo = (Data'*Data)\Data';

F = psudo*Y_data;


