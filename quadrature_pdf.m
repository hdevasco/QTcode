function xPDPairs = quadrature_pdf( state, phase, xVector, S )
%calculates quadrature probability density function
%   xPDPairs = quadrature_pdf(state, phase, xVector, S) will calculate the
%   quadrature probability density for state along direction given by
%   phase. The quadrature measurement has efficiency 1.  It returns a
%   2xlength(xVector) array; the nth row contains [xVector(n), the pdf at
%   xVector(n)]. state may be a density matrix or pure state vector in Fock
%   space.  S is the structre made by init_tables.

nX = length(xVector);
p=zeros(nX,1);
for n=1:nX
    h=homodyne_loss_measurement([phase,xVector(n)], 1, S, 'return vector');
    if isvector(state);
        p(n)=abs(h'*state)^2;
    else
        p(n)=h'*state*h;
    end
end
p = real(p);
xVector=xVector';

xPDPairs = [xVector,p];

end