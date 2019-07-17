clear all
close all
clc
%% load example multisubject FC matrices
load ../data/example_mats.mat
A = A(:,:,1:5);
[N,~,T] = size(A);  % number of nodes/layers
C = cell(T,1);      % create cell array
for i = 1:T         % populate array
    C{i} = A(:,:,i);
end
%% set some parameters
gammarange = [min(A(:)),max(A(:))]; % range of gamma values to consider
omegarange = [-0.5,1.5];            % range of omega values to consider
samplesPerStage = 250;              % number of samples per stage (ideally should be large)
nstage = 5;                         % number of stages (ideally should be large)
keep = 20;                          % number of nearest neighbors
couplingtype = 'categorical';       % nature of interlayer coupling (could also be ordered)
%% get initial estimate of boundaries
[gbound,obound,outputdata] = ...    % run boundary estimator
    fcn_get_bounds(C,gammarange,omegarange,samplesPerStage,nstage,keep,couplingtype);
%% sample from within those boundaries
samplesTotal = 1000;                % total number of samples
[S,G,O] = fcn_staged_multilayer(C,gbound,obound,samplesTotal,couplingtype);
%% calculate flexibility
switch couplingtype
    case 'categorical'
        [~,flx] = fcn_fastflx(S);
        flx = squeeze(nanmean(nanmean(flx,2),3));
    case 'ordered'
        flx = squeeze(nanmean(S(:,1:end - 1,:) ~= S(:,2:end,:),2));
end
%% calculate node entropy and nonsingleton communities
ent = zeros(N,samplesTotal);        % preallocate array for entropy
nonsingle = zeros(samplesTotal,1);  % preallocate array for nonsingletons
for isample = 1:samplesTotal        
    mx = max(max(S(:,:,isample)));  % max communities per sample
    bins = 1:mx;                    % bins
    h = hist(S(:,:,isample)',bins); % node histogram
    p = bsxfun(@rdivide,h,sum(h));  % probabilities
    e = -nansum(p.*log2(p));        % entropy
    ent(:,isample) = e;             % store
    h = hist(S(:,:,isample),bins);  % layer histogram
    nonsingle(isample,:) = ...      % number of nonsingletons
        nanmean(sum(h > 1));
end
%% make some figures
sz = 10;                            % size of points
f = fcn_rickplot([2,2,9,3]);        % create figure
vals = ...                          % concatenate variables
    [nanmean(ent)',squeeze(mean(max(S))),nonsingle];
for i = 1:3
    ax = axes('outerposition',[(i - 1)/3,0,1/3,1]);
    scatter(G,O,sz,vals(:,i),'filled');
end
%% if you use this code, please cite:
%
%  1.      Jutla et al (2011). A generalized Louvain method for community 
%          detection implemented in MATLAB.
%
%  2.      Betzel et al (2019). The community structure of functional brain
%          networks exhibits scale-specific patterns of variability across 
%          individuals and time. Neuroimage.