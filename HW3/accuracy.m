%% Vatsal Jain
% 605343009
function acc = accuracy(Y_pred, Y)
% Returns the accuracy, which is the percentage of correct predictions.
% Inputs:
%         Y_pred: a 10 x N array containing the predicted labels of each input image.
%         Y: a 10 x N array containing the actual labels of each input image.
% Output: 
%         acc: The percentage of accurately predicted samples.
    acc = mean(all(Y_pred == Y, 1));
end