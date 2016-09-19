function x = von_neumann_entropy(rho)
% von_neumann_entropy von-Neumann entropy
%    x = von_neumann_entropy(rho)
%
%    x = von_neumann_entropy(rho) calculates the von-Neumann entropy of
%    state rho, which may be a pure state vector or density matrix.

if isvector(rho)
    x = 0;
else
    e = eig(rho);
    x = -e'*log2(e+(e==0)); %0*log2(0)=0
end

x = real(x);

% author: Toby Cubitt, modified by Scott Glancy
% license: GPL2
%% Copyright (C) 2004-2009 Toby Cubitt
%%
%% This program is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public License
%% as published by the Free Software Foundation; either version 2
%% of the License, or (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
%% MA 02110-1301, USA.


