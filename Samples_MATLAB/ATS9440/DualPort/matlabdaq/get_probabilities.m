function [np,p] = get_probabilities(idx,histo)

for i = 1:length(idx)-1
    p(i)  = sum(histo(idx(i):idx(i+1)))/sum(histo);
    np(i) = i - 1;
end
p=p./sum(p);

