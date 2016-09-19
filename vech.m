function h = vech(m)
% h = vech(m)
% h is the column vector of elements on or below the main diagonal of m.
% m must be square.

t = tril(ones(size(m)));
h = m(logical(t));
