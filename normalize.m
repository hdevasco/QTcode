function state = normalize(unState, check)
% normalizes quantum state and warns if unnormalized
%   state = normalize(unState) returns a normalized state. unState
%   may be any unnormalized quantum state represented as a column
%   vector or density matrix in the Fock basis.
% 
%   state = normalize(unState, check) will check to see if unState
%   was already normalized.  If not, it will give a warning and
%   return a normalized state.  If check is 'check' (or any
%   character array), then it will
%   use default accuracy of 10^(-4).  If check is numeric, then
%   that number will be used for the accuracy.


if isvector(unState)
  normalization = real(unState'*unState);
  state = unState/realsqrt(normalization);
else
  normalization = trace(unState);
  state = unState/normalization;
end


if nargin==2
  if ischar(check)
    check = 10^(-4);
  end
  if abs(normalization-1) > check
     warning('fockSpace:tooSmall','state needs normalization.')
  end
end