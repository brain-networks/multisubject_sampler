function [gbound,obound,outputdata] = fcn_get_bounds(A,gammarange,omegarange,samplesPerStage,nstage,keep,couplingtype)
% %%
% gammarange = [-3,3];
% omegarange = [-10,50];
% samplesPerStage = 2500;
% nstage = 5;
% keep = 100;
%%
% load cellA
% A = A(1:4);
%%
outputdata.samplesPerStage = samplesPerStage;
outputdata.nstage = nstage;
outputdata.gammarange0 = gammarange;
outputdata.omegarange0 = omegarange;
%%
N = length(A{1});
T = length(A);
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
[mg,mo] = meshgrid(gammarange,omegarange);
mg = mg(:);
mo = mo(:);
k = boundary(mg(:),mo(:));
gbound = mg(k);
obound = mo(k);
%% for storing information
G = zeros(nstage*samplesPerStage,1);
O = G;
I = G;
S = zeros(N,T,nstage*samplesPerStage);
%% initial stages
tot = 0;
for istage = 1:nstage
    %% generate points within polygon
    
    gammarange = [min(gbound),max(gbound)];
    omegarange = [min(obound),max(obound)];
    g = unifrnd(gammarange(1),gammarange(2),samplesPerStage*100,1);
    o = unifrnd(omegarange(1),omegarange(2),samplesPerStage*100,1);
    in = inpolygon(g,o,gbound,obound);
    gammalistbig = g(in);
    omegalistbig = o(in);
    gammavals = gammalistbig(1:samplesPerStage);
    omegavals = omegalistbig(1:samplesPerStage);
    
%     count = 0;
%     gammavals = zeros(samplesPerStage,1);
%     omegavals = gammavals;
%     while count < samplesPerStage
%         gtest = unifrnd(gammarange(1),gammarange(2));
%         otest = unifrnd(omegarange(1),omegarange(2));
%         if inpolygon(gtest,otest,gbound,obound)
%             count = count + 1;
%             gammavals(count) = gtest;
%             omegavals(count) = otest;
%         end
%     end
    %% run modualrity maximization
    for isample = 1:samplesPerStage
        %% update counter
        tot = tot + 1;
        %%
        gamma = gammavals(isample);
        omega = omegavals(isample);
        B = (X - (gamma*Y)) + (Z*omega);
        s = genlouvain(B,[],false,true);
        s = reshape(s,[N,T]);
        S(:,:,tot) = s;
        %%
        flx = nanmean(nanmean(s(:,1:end - 1) ~= s(:,2:end)));
        num = mean(max(fcn_relabel_partitions(s)));
        isgood = flx > 0 & flx < 1 & num > 1 & num < N;
        I(tot) = isgood;
        G(tot) = gamma;
        O(tot) = omega;
        if mod(isample,25) == 0
            fprintf('stage = %i/%i, sample = %i/%i\n',istage,nstage,isample,samplesPerStage);
        end
    end
    %% find bad points near good ones
    idx = I(1:tot) == 1;
    g = G(1:tot);
    o = O(1:tot);
    gi = G(idx);
    oi = O(idx);
    gb = zeros(sum(idx)*keep,1);
    ob = gb;
    count = 0;
    for i = 1:sum(idx)
        dg = tiedrank(abs(g - gi(i)));
        do = tiedrank(abs(o - oi(i)));
        d = dg + do;
        d(idx) = inf;
        [~,jdx] = sort(d,'ascend');
        kdx = (count + 1):(count + keep);
        gb(kdx) = g(jdx(1:keep));
        ob(kdx) = o(jdx(1:keep));
        count = count + keep;
    end
    k = boundary(gb,ob);
    gbound = gb(k);
    obound = ob(k);
    scatter(g,o,100,I(1:tot)); hold on; scatter(gb,ob,50','ro','filled'); hold off; drawnow;
end
outputdata.gamma = G;
outputdata.omega = O;
outputdata.S = S;
outputdata.isgood = I;
outputdata.params.keep = keep;