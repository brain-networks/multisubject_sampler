function [S,G,O] = fcn_staged_multilayer(A,gbound,obound,nsamples,couplingtype)
N = length(A{1});
T = length(A);
%%
gammarange = [min(gbound),max(gbound)];
omegarange = [min(obound),max(obound)];
g = unifrnd(gammarange(1),gammarange(2),nsamples*100,1);
o = unifrnd(omegarange(1),omegarange(2),nsamples*100,1);
in = inpolygon(g,o,gbound,obound);
gammalistbig = g(in);
omegalistbig = o(in);
%%
X = spalloc(N*T,N*T,N*N*T+2*N*T);
Y = X;
switch couplingtype
    case 'categorical'
        all2all = N*[(-T+1):-1,1:(T-1)];
        Z = spdiags(ones(N*T,2*T-2),all2all,N*T,N*T);
    case 'ordered'
        Z = spdiags(ones(N*T,2),[-N,N],N*T,N*T);
end
%%
for s = 1:T
    idx = (1:N) + (s - 1)*N;
    X(idx,idx) = A{s};
    Y(idx,idx) = ~eye(N);
end
X = (X + X')/2;
Y = (Y + Y')/2;
%%
sample = 0;
S = zeros(N,T,nsamples);
G = zeros(nsamples,1);
O = G;
tries = 0;
tic;
count = 0;
while sample < nsamples
    count = count + 1;
    gamma = gammalistbig(count);
    omega = omegalistbig(count);
    B = (X - (gamma*Y)) + (Z*omega);
    s = genlouvain(B,[],false,true);
    s = reshape(s,[N,T]);
    flx = nanmean(nanmean(s(:,1:end - 1) ~= s(:,2:end)));
    num = mean(max(fcn_relabel_partitions(s)));
    isgood = flx > 0 & flx < 1 & num > 1 & num < N;
    if isgood
        sample = sample + 1;
        S(:,:,sample) = s;
        G(sample) = gamma;
        O(sample) = omega;
        if mod(sample,1) == 0
            fprintf('sample %i/%i, %i try, elapsed time of %.2f s\n',sample,nsamples,tries,toc);
            plot(G(1:sample),O(1:sample),'ko',gbound,obound,'r');
            drawnow;
        end
        tries = 0;
    end
    
end