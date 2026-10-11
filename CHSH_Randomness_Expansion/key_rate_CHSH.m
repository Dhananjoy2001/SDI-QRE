%% Main Code

function [H,lambda0,lambda1] = hae_chsh(betaEff,m,N)


% Parameters

d  = 4;
nA = 2;
[w,t] = gauss_radau(m);
t(t<1e-16) = 1e-16;


% Cycle measurements

Proj = cycle_projectors_even(N);


% Generation Measurement

xstar = 1;
Gen{1} = Proj(:,:,xstar);
Gen{2} = eye(d) - Proj(:,:,xstar);

% Precompute real and imaginary parts of Gen{1} for real arithmetic
Gen_re = real(Gen{1});
Gen_im = imag(Gen{1});



% YALMIP variables (Purely Real Formulation to preserve duals)

yalmip('clear')


% Density matrix sigma = X + i*Y (X symmetric, Y skew-symmetric)

X = sdpvar(d, d, 'symmetric', 'real');
Y = sdpvar(d, d, 'skew', 'real');

for a = 1:nA
    for i = 1:m
        % Decompose complex zeta, eta, theta into real blocks:
        % z = z_re + i*z_im (full matrices)

        z_re{a,i}   = sdpvar(d, d, 'full', 'real');
        z_im{a,i}   = sdpvar(d, d, 'full', 'real');

        % et = et_re + i*et_im (Hermitian -> et_re symmetric, et_im skew)

        et_re{a,i}  = sdpvar(d, d, 'symmetric', 'real');
        et_im{a,i}  = sdpvar(d, d, 'skew', 'real');

        % th = th_re + i*th_im (Hermitian -> th_re symmetric, th_im skew)

        th_re{a,i}  = sdpvar(d, d, 'symmetric', 'real');
        th_im{a,i}  = sdpvar(d, d, 'skew', 'real');
    end
end



% Objective (Using expanded real trace components)

objective = 0;
for i = 1:m
    for a = 1:nA
        zR = z_re{a,i}; zI = z_im{a,i};
        eR = et_re{a,i}; eI = et_im{a,i};
        tR = th_re{a,i}; tI = th_im{a,i};


        % Real part of sym_part: (z + z') + (1-t(i))*et

        sym_part_re = zR + zR' + (1 - t(i)) * eR;
        sym_part_im = zI - zI' + (1 - t(i)) * eI;  % Note: z' transpose flips skew property for imag part

        if a == 1
            G_re = Gen_re; G_im = Gen_im;
        else
            G_re = eye(d) - Gen_re; G_im = -Gen_im;
        end

        % Trace of product expansion: Re(Tr(A * B)) = Tr(A_re*B_re - A_im*B_im)

        tr_gen_sym = trace(G_re * sym_part_re - G_im * sym_part_im);
        tr_th      = trace(tR); % Imaginary part of symmetric/hermitian trace is zero

        objective = objective + (w(i)/(t(i)*log(2))) * tr_gen_sym + (w(i)/log(2)) * tr_th;
    end
end



% Constraints

Constraints = [];
ConNorm = (trace(X) == 1);
Constraints = [Constraints, ConNorm];

SigmaPSD = [X -Y;
            Y  X] >= 0;

Constraints = [Constraints, ConNorm, SigmaPSD];


% Complex PSD constraint [sigma z; z' et] >= 0 mapped to real symmetric block matrix:
% [ X  -Y  |  zR -zI ]
% [ Y   X  |  zI  zR ]
% [-----------------]
% [ zR' zI |  eR -eI ]
% [ -zI zR |  eI  eR ]

for i = 1:m
    for a = 1:nA
        zR = z_re{a,i}; zI = z_im{a,i};
        eR = et_re{a,i}; eI = et_im{a,i};
        tR = th_re{a,i}; tI = th_im{a,i};

        Gamma1 = [  X,    -Y,    zR,   -zI;
                    Y,     X,    zI,    zR;
                   zR',   zI',   eR,   -eI;
                  -zI',   zR',   eI,    eR ];

        Gamma2 = [  X,    -Y,    zR',   zI';
                    Y,     X,   -zI',   zR';
                   zR,   -zI,    tR,   -tI;
                   zI,    zR,    tI,    tR ];

        Constraints = [Constraints, Gamma1 >= 0, Gamma2 >= 0];
    end
end



% Cycle inequality & Bell constraint

beta_cycle = 0;
for x = 1:N-1
    P00_tr = trace((eye(d) - Proj(:,:,x)) * (eye(d) - Proj(:,:,x+1)) * (X + 1i*Y));
    P11_tr = trace(Proj(:,:,x) * Proj(:,:,x+1) * (X + 1i*Y));
    beta_cycle = beta_cycle + real(P00_tr) + real(P11_tr);
end
P01_tr = trace((eye(d) - Proj(:,:,N)) * Proj(:,:,1) * (X + 1i*Y));
P10_tr = trace(Proj(:,:,N) * (eye(d) - Proj(:,:,1)) * (X + 1i*Y));
beta_cycle = (beta_cycle + real(P01_tr) + real(P10_tr)) / N;


ConCycle = (beta_cycle == betaEff);
Constraints = [Constraints, ConCycle];



% Solve SDP

ops = sdpsettings('solver', 'mosek', 'verbose', 1);
sol = optimize(Constraints, objective, ops);
if sol.problem ~= 0
    disp(sol.info)
    H       = NaN;
    lambda0 = NaN;
    lambda1 = NaN;
    return
end



% Dual variables (Now fully functional without NaNs)

lambda0 = value(dual(ConNorm));
lambda1 = value(dual(ConCycle));



% Entropy bound

obj_val = real(value(objective));
cm = sum(w./t)/log(2);
H = cm + obj_val;


% Results

fprintf('\n');
fprintf('N              = %d\n',N);
fprintf('Cycle value    = %.12f\n',value(beta_cycle));
fprintf('Objective      = %.12f\n',obj_val);
fprintf('c_m            = %.12f\n',cm);
fprintf('H              = %.12f\n',H);
fprintf('lambda0        = %.12f\n',lambda0);
fprintf('lambda1        = %.12f\n',lambda1);
f_beta = cm - lambda0 - lambda1*betaEff;
fprintf('f(betaEff)     = %.12f\n',f_beta);
end




%% Code for Key rate

function r = RE_rate(betaEff,q,N,n,lambda,m,epsEA,epss,gamma)


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



%% Binary entropy

function h = binary_entropy(x)

if x == 0 || x == 1
    h = 0;
else
    h = -x*log2(x) - (1-x)*log2(1-x);
end

end
