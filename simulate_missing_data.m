function [Y_hidden] = simulate_missing_data(Y)
C = load('../../data_upmc/percent_missing.mat');
if(length(C.percent_missing_top_data) ~= length(Y))
    fprintf('The data on percent missing is not the same length as Y!\n');
    throw(MException('simulate_missing_data:unequalLength','Y is not 4 times the length of missing data'));
end
% missing_prob = ones(4,length(C.percent_missing_top_data));
% missing_prob(1,:) = C.percent_missing_top_data;
% for i = 1:4
%     missing_prob(:,i) = missing_prob(1,i);
% end
% missing_prob = reshape(missing_prob,1,length(Y));
missing_prob = C.percent_missing_top_data;
r = rand(size(missing_prob));
keep = r > missing_prob;
Y_hidden = Y;
Y_hidden(:,~keep) = NaN;