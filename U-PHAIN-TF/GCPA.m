function [x] = GCPA(param, paramsolver, oracle)
%% initialization
x = paramsolver.x0;
y = paramsolver.y0;
z = paramsolver.z0;

tau = paramsolver.tau;
sigma = paramsolver.sigma;
alpha = paramsolver.alpha;
eta = paramsolver.eta;

%% iteration

for i = 1:paramsolver.I

    tmp_r = y + eta*param.G(x - tau*(param.L_adj(z)+param.G_adj(y)));
    r = tmp_r - eta*param.proj(tmp_r/eta);
    
    p = x - tau*(param.L_adj(z) + param.G_adj(r));

    tmp_q = z + sigma*param.L(2*p - x);
    q = tmp_q - sigma*param.prox(tmp_q/sigma);
    
    x = x + alpha*(p - x);
    y = y + alpha*(r - y);
    z = z + alpha*(q - z);
    
end
