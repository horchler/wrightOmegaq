function Z=wrightOmegaqInv(W)
%WRIGHTOMEGAQINV  Inverse of the Wright omega function.
%   Z = WRIGHTOMEGAQINV(W) evaluates the inverse of the Wright omega function
%   numerically. W is an array and may be complex. If W is an array of symbolic
%   values, it is converted to double-precision for computation and then recast
%   as symbolic.
%
% 	The inverse Wright omega function is defined as (Corless & Jeffrey, 2002):
%       Z = W+LOG(W)-2*pi*1i,               -Inf < REAL(W) < -1, IMAG(W) = 0
%       Z = -1+/-pi*1i (-1+pi*1i returned),	W = -1
%       Z = W+LOG(W),                       otherwise
%
%   Note that Z = -1+pi*1i is returned here for the singularity at REAL(W) = -1
%   as this is a solution to W+LOG(W) = Z. The inverse of the Wright omega
%   function is not single-valued for all W, even though the Wright Omega
%   function itself is single-valued for all Z. Additionally, note that Z = -Inf
%   for W = 0.
%
%   Class support for Z:
%       float: double, single
%       symbolic
%
%   See also: WRIGHTOMEGAQ, WRIGHTOMEGA, LAMBERTW

%	Robert M. Corless and David J. Jeffrey, "The Wright omega Function," In:
%   Artificial Intelligence, Automated Reasoning, and Symbolic Computation,
%   Joint International Conferences, AISC 2002 and Calculemus 2002, Marseille,
%   France, July 2002, (Jacques Calmet, Belaid Benhamou, Olga Caprotti, Laurent
%   Henocque, and Volker Sorge, Eds.), Berlin: Springer-Verlag, pp. 76-89, 2002.
%   http://orcca.on.ca/TechReports/2000/TR-00-12.html

%   Andrew D. Horchler, horchler @ gmail . com, Created 3-11-13
%   Revision: 1.0, 3-12-12


% Convert symbolic input, converted back at end
isSym = isa(W,'sym') || ischar(W);
if isSym
    try
        W = double(sym(W));
    catch ME
        if strcmp(ME.identifier,'symbolic:sym:double:cantconvert')
            error('SHCTools:wrightOmegaqInv:InvalidSymbolicW',...
                 ['Symbolic input W must contain only numeric values and no '...
                  'expressions containing variables.'])
        else
            rethrow(ME);
        end
    end
elseif ~isfloat(W)
    error('SHCTools:wrightOmegaqInv:InvalidW',...
          'Input W must be an array of floating point or symbolic values.');
end

if isempty(W) || all(isnan(W(:)))
    Z = W;
else
    % Support for single precision: single(pi) ~= pi
	dataType = class(W);
    
    Y = imag(W);    
    if isscalar(W)
        % Special values, Z(1) used in order retain datatype
        if W == 0.5671432904097838                	% Omega constant
            Z(1) = 0;
        elseif W > 2^59                             % W self-saturates: X > 2^59
            Z = W;
        elseif W < -1 && Y == 0
            Z = W+log(W)-2*cast(pi,dataType)*1i;
        else
            Z = W+log(W);                           % Handles singularities too
        end
    else
        if (isreal(W) || all(Y(:) == 0)) && all(W(:) >= 0)
            W = real(W);
            Z = NaN(size(W),dataType);
            
            iw = (W ~= 0.5671432904097838 & W <= 2^59);
            Z(iw) = W(iw)+log(W(iw));               % Handles singularities too
        else
            Z = complex(NaN(size(W),dataType));
            
            iw = (W < -1 & Y == 0);
            Z(iw) = W(iw)+log(W(iw))-2*cast(pi,dataType)*1i;
            
            iw = (~iw & W ~= 0.5671432904097838 & W <= 2^59);
            Z(iw) = W(iw)+log(W(iw));               % Handles singularities too
        end
        
        % Special values
        iw = (W == 0.5671432904097838);
      	Z(iw) = 0;                                  % Omega constant
        
        iw = (W > 2^59);                            % W self-saturates: X > 2^59
        Z(iw) = W(iw);
    end
end

% Reconvert symbolic input
if isSym
    Z = sym(Z,'d');
end