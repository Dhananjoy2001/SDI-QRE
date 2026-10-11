

function M = cycle_projectors_even(N)

sx = [0 1;
    1 0];

sy = [0 -1i;
    1i 0];

sz = [1 0;
    0 -1];

I2 = eye(2);

Pi_single = zeros(2,2,N);

for x = 0:N-1

    S = [cos(pi*x/N);
        sin(pi*x/N);
        0];

    Pi_single(:,:,x+1) = 0.5*( ...
        I2 ...
        + S(1)*sx ...
        + S(2)*sy ...
        + S(3)*sz );

end

M = zeros(4,4,N);

for x = 0:N-1

    if mod(x,2)==0
      M(:,:,x+1) = kron(conj(Pi_single(:,:,x+1)),I2);
    else
      M(:,:,x+1) = kron(I2,Pi_single(:,:,x+1));
    end

end

end