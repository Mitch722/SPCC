%% test:

e = 0: 0.01 : 1;

test = zeros(size( e ));
% Work using betaln
for i = 1 : length( e )
   
    b = 1000;
    q_fun = e(i)*b;
    
    a = q_fun - 2;
    % c is equivalent to 1 -  epsilon
    c = q_fun / b;
    
    test(i) = binopdf( a,b,c );
    
    
end

plot(e, test)