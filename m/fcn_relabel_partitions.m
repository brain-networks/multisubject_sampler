function cinew = fcn_relabel_partitions(ci)
[n,m] = size(ci);

cinew = zeros(n,m);
for i = 1:m
    c = ci(:,i);
    d = zeros(size(c));
    count = 0;
    while sum(d ~= 0) < n
        count = count + 1;
        ind = find(c,1,'first');
        tgt = c(ind);
        rep = c == tgt;
        d(rep) = count;
        c(rep) = 0;
    end
    cinew(:,i) = d;
end