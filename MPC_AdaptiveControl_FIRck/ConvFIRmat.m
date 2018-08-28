%% ConFirmat builds a convolutional matrix for prediction output of an FIR model
%
% Inputs:  F: FIR coefficients 
%          p: length of output series for constraints 
%
% Output:  F_big: convolution matrix with constraints
%
function F_big = ConvFIRmat(F, ncc)

p = ncc+1-length(F(1, :)); % put in the new length takes 
% no_outputs: number of values in output vector y
% Flen length of the F coefficients
[no_outputs, Flen] = size(F);   
% the convolutional matrix of F
bigF = zeros(no_outputs*p, Flen + p - 1);
bigF(1:no_outputs, 1:Flen) = F;
% intermediate Fs
interbigF = bigF;
% start counter
i = 1;
% while loop
while i < p
    % shift intermediate no_outputs down and 1 across
    interbigF = circshift(interbigF, [no_outputs, 1]);
    bigF = bigF + interbigF;
    % update counter
    i = i + 1;
    
end
% the convolutional matrix 
F_big = bigF;
