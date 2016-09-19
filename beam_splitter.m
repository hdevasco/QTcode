function bsb = beam_splitter( eta, S )
%beam_splitter matrix for beam splitter in fock basis
%   bsb = beam_splitter(eta, S)
%
%   bsb = beam_splitter(eta, S) will calculate the beam splitter unitary
%   evolution matrix in the Fock basis.  The beam splitter's transmissivity
%   is eta.  The matrix acts on the direct product space of the two input
%   modes.  To apply beam splitter to density matrices rho1 and rho2, first
%   use [rho12, S] = double_kron(rho1, rho2).  Then rho12 = bsb*rho12*bsb'

b=[sqrt(eta),sqrt(1-eta);-sqrt(1-eta),sqrt(eta)];

doublePhotons = 2*S.photons;
if doublePhotons >= 171
    err = MException('TooMany:Photons', ...
        'Hilbert space contains too many photons.  Cannot construct beam splitter');
    throw(err)
end
factorialList = cumprod(1:doublePhotons);
factorialList = [1,factorialList];

% beam splitter operator in four dimensional space labeled by photon
% numbers in both input and output modes.
bs = zeros(S.dimHilbertSpace,S.dimHilbertSpace,S.dimHilbertSpace,S.dimHilbertSpace);
for a1=0:S.photons
    for a2=0:S.photons
        for n1=0:S.photons
            n2 = a1+a2-n1;
            if n2>=0 && n2<= S.photons
                s = 0;
                for k1 = 0:n1
                    if a1-k1>=0 && a1-k1<=n2
                        s = s + S.binom(n1+1,k1+1)*S.binom(n2+1,a1-k1+1)*b(1,1).^k1*b(2,1).^(n1-k1)*b(1,2).^(a1-k1)*b(2,2).^(n2-a1+k1);
                    end
                end
                bs(a1+1,a2+1,n1+1,n2+1) = realsqrt(factorialList(a1+1)*factorialList(n1+n2-a1+1)/(factorialList(n1+1)*factorialList(n2+1)))*s;
            end
        end
    end
end

% transform beam splitter into 2 dimensional kronecker product space
bsb = zeros(S.dimHilbertSpace^2);
for k = 1:S.dimHilbertSpace^2
    for j = 1:S.dimHilbertSpace^2
        a1 = floor((j-1)/S.dimHilbertSpace);
        a2 = mod(j-1,S.dimHilbertSpace);
        n1 = floor((k-1)/S.dimHilbertSpace);
        n2 = mod(k-1,S.dimHilbertSpace);
        bsb(j,k) = bs(a1+1,a2+1,n1+1,n2+1);
    end
end

% TODO: the above could probably be done with some kind of reshaping of the
% four dimensional bs

end

