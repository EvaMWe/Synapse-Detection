function selected = discrimination(parameter)

sig0 = parameter.sigm0;
sig2 = parameter.sigm2;
moments = parameter.moments;
totalN = parameter.regionNb;

factor = 1/(2*pi*sig0*sig2*totalN);


m0q = zeros(totalN-1,1);
m2q = zeros(totalN-1,1);
S = zeros(totalN,1);

for p = 1:totalN
    
    m0p = moments(p,1);
    m2p = moments(p,2);
    
    m0q(1:p-1,1) = moments(1:p-1,1);
    m0q(p:end,1) = moments(p+1:end,1);
    
    
    m2q(1:p-1,1) = moments(1:p-1,2);
    m2q(p:end,1) = moments(p+1:end,2);
    
    S(p,1) = factor*(sum(exp(-(m0q - m0p).^2./(2*sig0)-(m2q-m2p).^2./(2*sig2))));
end

pth = 0.95;
[cnts,edge,~] = histcounts(S,400);
l = length(cnts);
k = 1;
while sum(cnts(l-k:l))/sum(cnts) < pth
    k = k + 1;
end
thresh =edge(l-k+1);

selected = find(S >= thresh);
