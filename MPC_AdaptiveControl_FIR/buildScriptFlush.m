
review_k = 10;
flush_time = 1000;

k_init = k;
flush_state = 0;

if mean( abs( u(:, k-review_k: k-1) )) > 0.75*100 || flush_state == 1
    % allows if statement to be completed 
    flush_state = 1;
    
    maxF = 100;
    % MPC prediction window
    p = 20;
    [H, f, Ac, Ax, b1, lb, ub, opt] = MPC_vars(A-B*K_opt, B, C, K_opt, R, p, main_bounds, maxF);
    % cholesky for mpcqpsolver
    [L2, ~] = chol(H,'lower');
    Linv = inv(L2);

    % options for mpcqpsolver:
    options = mpcqpsolverOptions;
    
    k0 = k/ ratioTs;
    
    X = xhat(:, k0);

    b = b1 + Ax*X;
    
    ck = mpcqpsolver(Linv, f, Ac, b, [], zeros(0,1), false(size(b)), options);
    % ck = quadprog(H, f, -Ac, -b, [], [], lb, ub, [], options);

    if isempty(ck) || abs(ck(1)) > 10
        c = 0;
    else
        c = ck(1);
    end
    
    if k > k_init + flush_time
        
        flush_state = 0;
    end
    
    
else
    flush_state = 0;
    k_init = k;
    
    if k/25 == round(k/25)
            % selects new constraints
            [~, ~, RstarModel, entry, ~, H1] = ACSA_FIR(M_samps, y(:, 1:k-1), Ck(:, 1:k-1), p, params, bnds, A-B*K_opt, B, C, K_opt);
                       
    end   

        b2 = cell2mat(RstarModel(entry.b0, :)');
        Dc2 = cell2mat(RstarModel(entry.Dc, :)');
        Dx2 = cell2mat(RstarModel(entry.Dx, :)');
        
%         [H2, f, Dc2, Dx2, b2, lb, ub, op] = MPC_varsFIR(A-B*K_opt, B, C, K_opt, R, p, bnds, maxF, aveFir);
        [L21, ~] = chol(H1,'lower');
        Linv1 = inv(L21);
        no_coefs = params.n;
        % X = xhat(:, k0);
        chat = Ck(k-no_coefs:k-1)';
        b = b2 + Dx2*chat;
                
        ck = mpcqpsolver(Linv1, zeros(length(Linv1), 1), Dc2, b, [], zeros(0,1), false(size(b)), options);
        if isempty(ck) || abs(ck(1)) > 10
            c = 0;
        else
            c = ck(1);
        end
end




