function [mu,F] = fcn_fastflx(S)
% computes pairwise flexibility on node x layer x repetition matrix, S
%
%   inputs:
%       S,      node x layer x repetition matrix of community labels
%
%   outputs:
%       mu,     node x 1 vector of average pairwise flexibility 
%       F,      node x layer x layer x repetition distance matrix
%
%   Richard Betzel, University of Pennsylvania, 2018
%
[N,T,NReps] = size(S);
F = zeros(N,T,T,NReps);
for iRep = 1:NReps
    X = repmat(S(:,:,iRep),[1,1,T]);
    Y = permute(X,[1,3,2]);
    F(:,:,:,iRep) = double(X ~= Y);
end
for i = 1:T
    F(:,i,i,:) = nan;
end
mu = nanmean(nanmean(nanmean(F,2),3),4);