%% Code for smooth min entropy

function [H,lambda0,lambda1] = hae_kcbs(betaQ,delta,m,N)


% Effective Quantum violation

betaEff = betaQ - delta;


% KCBS N-cycle randomness SDP

d  = 3;
nA = 2;

[w,t] = gauss_radau(m);
t(t<1e-16) = 1e-16;

Pi = kcbs_projectors(N);



% Generation measurement

xstar = 1;

M{1} = Pi(:,:,xstar);
M{2} = eye(d) - Pi(:,:,xstar);



% YALMIP variables

yalmip('clear')

sigma = sdpvar(d,d,'hermitian','real');

for a = 1:nA
    for i = 1:m

        zeta{a,i}  = sdpvar(d,d,'full','real');
        eta{a,i}   = sdpvar(d,d,'hermitian','real');
        theta{a,i} = sdpvar(d,d,'hermitian','real');

    end
end



% Objective

objective = 0;

for i = 1:m

    for a = 1:nA

        z  = zeta{a,i};
        et = eta{a,i};
        th = theta{a,i};

        sym_part = z + z' + (1-t(i))*et;

        objective = objective ...
            + (w(i)/(t(i)*log(2))) ...
            * real(trace(M{a}*sym_part)) ...
            + (w(i)/log(2)) ...
            * real(trace(th));

    end
end



% Constraints

Constraints = [];

ConNorm = (trace(sigma) == 1);

Constraints = [Constraints,...
    ConNorm,...
    sigma >= 0];



% KCBS inequality

kcbs_expr = 0;

for k = 1:N
    kcbs_expr = kcbs_expr + real(trace(Pi(:,:,k)*sigma));
end

kcbs_expr = kcbs_expr/N;

ConKCBS = (kcbs_expr == betaEff);

Constraints = [Constraints,...
    ConKCBS];



% PSD Blocks

for i = 1:m

    for a = 1:nA

        z  = zeta{a,i};
        et = eta{a,i};
        th = theta{a,i};

        Gamma1 = [sigma z;
                  z'    et];

        Gamma2 = [sigma z';
                  z     th];

        Constraints = [Constraints,...
            Gamma1 >= 0,...
            Gamma2 >= 0];

    end
end



% SDP Solver

ops = sdpsettings('solver','mosek','verbose',1);

sol = optimize(Constraints,objective,ops);

if sol.problem ~= 0

    disp(sol.info)

    H       = NaN;
    lambda0 = NaN;
    lambda1 = NaN;

    return

end



% Dual variables

lambda0 = value(dual(ConNorm));
lambda1 = value(dual(ConKCBS));



% Entropy bound

obj_val = value(objective);

cm = sum(w./t)/log(2);

H = cm + obj_val;



% Results

fprintf('\n');
fprintf('N              = %d\n',N);
fprintf('KCBS value     = %.12f\n',value(kcbs_expr));
fprintf('Objective      = %.12f\n',obj_val);
fprintf('c_m            = %.12f\n',cm);
fprintf('H              = %.12f\n',H);

fprintf('lambda0        = %.12f\n',lambda0);
fprintf('lambda1        = %.12f\n',lambda1);

f_beta = cm - lambda0 - lambda1*betaEff;

fprintf('f(betaEff)     = %.12f\n',f_beta);

end





%% Code for Key rate

function r = RE_rate(betaQ,delta,q,N,n,lambda,m,epsEA,epss,gamma)


% Effective Quantum violation

betaEff = betaQ - delta;


% Quadrature constant c_m

[w,t] = gauss_radau(m);

cm = sum(w./t)/log(2);



% Min-tradeoff function

omegaQ = [1 betaEff];

f = (1-q)*(cm - lambda*omegaQ');



% Bounds on g

lambda_min = min(lambda);
lambda_max = max(lambda);

gmax = (1-q)*(cm - lambda_min);
gmin = (1-q)*(cm - lambda_max);



% Variance bound

varf = (gmax-gmin)^2/q;



% EAT constants

V = (log(2)/2) ...
    *(log2(9) + sqrt(varf + 2))^2;

K = (1/(6*(1-gamma)^3*log(2))) ...
    * 2^(gamma*(1+gmax-gmin)) ...
    * log(2^(1+gmax-gmin) + exp(2));



% Smooth min-entropy

Hmin = n*f ...
      -(1/gamma)*(1 - 2*log2(epsEA*epss)) ...
      -n*(gamma*V + gamma^2*K);



% Input randomness

lin = n*(binary_entropy(q) + q*log2(N)) + 6;



% Net randomness rate

r = (Hmin - lin)/n;



% Results

fprintf('\n');
fprintf('f(betaEff)     = %.10f\n',f);
fprintf('var(f)         = %.10f\n',varf);
fprintf('V              = %.10f\n',V);
fprintf('K              = %.10f\n',K);
fprintf('Hmin/n         = %.10f\n',Hmin/n);
fprintf('Input cost/n   = %.10f\n',lin/n);

end





%% KCBS orthogonal projectors

function Pi = kcbs_projectors(N)

% Odd N-cycle KCBS realization

if mod(N,2)==0
    error('KCBS realization requires odd N');
end

c = sqrt(cos(pi/N)/(1+cos(pi/N)));
s = sqrt(1-c^2);

Pi = zeros(3,3,N);

for j = 1:N

    phi = j*pi*(N-1)/N;

    u = [...
        c;
        s*sin(phi);
        s*cos(phi)];

    Pi(:,:,j) = u*u';

end

end




%% Gauss Radau quadrature

function [w,t] = gauss_radau(m)

J = zeros(m,m);

for n = 1:m-1

    J(n,n) = 1/2;

    J(n,n+1) = n/(2*sqrt(4*n^2-1));

    J(n+1,n) = J(n,n+1);

end

J(m,m) = (3*m-1)/(4*m-2);

[V,D] = eig(J);

w = (V(1,:).^2)';
t = diag(D);

end




%% Binary Entropy

function h = binary_entropy(x)

if x == 0 || x == 1
    h = 0;
else
    h = -x*log2(x) - (1-x)*log2(1-x);
end

end
