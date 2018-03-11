try 
    
    a  = y;
    
catch
    
    load('output_input.mat')
    
end

[ck, Mmodels] = ACSA_1(M, y, u, p, params, Q_bar, bnds, Ts);