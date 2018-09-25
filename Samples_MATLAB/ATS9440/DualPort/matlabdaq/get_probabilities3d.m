function [np,p] = get_probabilities3d(idx,histo)

for i = 1:length(idx(1,:))-1
    for k = 1:length(idx(2,:))-1
    p(i,k)  = sum(sum(histo(idx(1,i):idx(1,i+1),idx(2,k):idx(2,k+1))))/sum(sum(histo));
    np(i,k) = [i - 1 + k-1];    
    end
end
p=p./sum(sum(p));