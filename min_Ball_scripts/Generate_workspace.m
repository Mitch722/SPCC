% Try to load previous workspace
% if no previously stored workspace is avaliable then make and save a new
% one

try load('out_put_data.mat')

catch 
    disp('No output from previous run')
    
    multisample
    save('out_put_data.mat')
end