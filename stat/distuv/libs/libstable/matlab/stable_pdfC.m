function pdf = stable_pdfC(x, pars, parametrization)
% pdf = stable_pdfC(x, pars, parametrization)
% Code for computing the PDF of an alpha-estable distribution.
% Expresions presented in [1] are employed.
%  Inputs:
%    pars: parameters of the alpha-estable distributions.
%                 params=[alpha,beta,sigma,mu];
%
%    x:    vector of evaluataion points.
%
%    parametrization:  parameterization employed (view [1] for details)
%             0:   mu=mu_0
%             1:   mu=mu_1
%
% [1] Nolan, J. P. Numerical Calculation of Stable Densities and
%     Distribution Functions Stochastic Models, 1997, 13, 759-774
%
% Copyright (C) 2013. Javier Royuela del Val
%                     Federico Simmross Wattenberg

% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; version 3 of the License.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program; If not, see <http://www.gnu.org/licenses/>.
%
%
%  Javier Royuela del Val.
%  E.T.S.I. Telecomunicación
%  Universidad de Valladolid
%  Paseo de Belén 15, 47002 Valladolid, Spain.
%  jroyval@lpi.tel.uva.es    
%

if nargin < 3
    parametrization = 0;
end
if nargin < 2
    pars = [2, 0, 1, 0];
end
if length(pars) ~= 4
    error('pars must be a four elements vector');
end
if (parametrization<0 || parametrization >1)
    error('parametrization must be 0 or 1');
end

dist=calllib('libstable','stable_create',pars(1),pars(2),pars(3),pars(4),parametrization);

calllib('libstable','stable_set_THREADS',0);
calllib('libstable','stable_set_relTOL',1e-12);
calllib('libstable','stable_set_absTOL',1e-12);

n=length(x);
pdf=zeros(1,n);
[~,~,pdf,~]=calllib('libstable','stable_pdf',dist,x,n,pdf,[]);

calllib('libstable','stable_free',dist);

end
