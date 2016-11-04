function f = wig(r, k, rho, M, sq_f_tab, sq_p_tab)
% helper function for wigner()

if(k < 0)
    f = conj(wig(r, -1 * k, rho, M, sq_f_tab, sq_p_tab));
else
    w = zeros(1, M - k + 2);
%    w(2) = (1/pi) .* (factorial(k) .^ (-1/2)) .* ((r .* sqrt(2)) .^ k) .* exp(-1 .* (r .^ 2));
    w(2) = (1/pi) .* sq_f_tab(k+1) .* ((r .* sqrt(2)) .^ k) .* exp(-1 .* (r .^ 2));
    for n = 1:(M-k)
%        w(n+2) = ((n * (n+k)) ^ (-1/2)) * ((2 * r^2 + 1 - k - 2*n) * w(n+1) - sqrt((n-1)*(n-1+k)) * w(n));
        w(n+2) = (1 / sq_p_tab(n+1,k+1)) * ((2 * r^2 + 1 - k - 2*n) * w(n+1) - sq_p_tab(n,k+1) * w(n));
    end
    
    f = sum( w(2:(M-k+2)) .* (diag( rho((k+1):(M+1),1:(M-k+1)) )).' );
end