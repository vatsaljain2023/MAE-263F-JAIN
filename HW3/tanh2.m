%% Vatsal Jain
% 605343009

function A = tanh2(Z)
% Tanh applies the tanh activation function on the input to get a nonlinear output.
%    Inputs:
%         X: A M x N matrix representing the output of the neurons, which serves as the input of the activation function. M is the number of neurons, and N is the number of examples
%   Outputs:
%         Z: a M x N matrix representing the output after the tanh activation function

    A = 2./(1+exp(-2*Z))-1;
end