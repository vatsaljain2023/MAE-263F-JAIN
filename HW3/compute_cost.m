%% Vatsal Jain
% 605343009

function cost = compute_cost(AL, Y)
    cost = -sum(Y .* log(AL));
end