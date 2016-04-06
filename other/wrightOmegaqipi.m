function W=wrightOmegaqipi(Z)
%WRIGHTOMEGAQIPI  Wright omega function, solution of the equation W+LOG(W) = Z.
%   W = WRIGHTOMEGAQIPI(Z) performs floating point evaluation of the Wright
%   omega function for complex inputs that return a real-valued solution, i.e.,
%   along the "lines of discontinuity" Z = X+i*pi and Z = X-i*pi for X <= -1.
%   The input, Z, must have this form and may be an array of complex values. The
%   output, W, will always be real and have the same dimensions as Z.
%
%   Example:
%       % Plot absolute value of the function
%     	x = -20:0.01:-1; w = wrightOmegaqipi([x+pi*1i;x-pi*1i]);
%       figure; semilogy(x,abs(w)); xlabel('X'); ylabel('|\omega(X\pmi\pi)|');
%       title('Wright \omega Function Along Lines of Discontinuity');
%
%   Note:
%       Due to numerical precision and the nature of this equation, the inverse
%       of WRIGHTOMEGAQIPI, may not always return a value close (in an absolute
%       sense) to Z, e.g., Z <= -714.84989998498+pi*1i.
%
%   Class support for Z:
%       float: double, single
%
%   See also: WRIGHTOMEGAQ, WRIGHTOMEGA, LAMBERTW

%   Based on:
%
%   [1] Piers W. Lawrence, Robert M. Corless, and David J. Jeffrey, "Algorithm
%   917: Complex Double-Precision Evaluation of the Wright omega Function," ACM
%   Transactions on Mathematical Software, Vol. 38, No. 3, Article 20, pp. 1-17,
%   Apr. 2012. http://dx.doi.org/10.1145/2168773.2168779
%
%	[2] Robert M. Corless and David J. Jeffrey, "The Wright omega Function," In:
%   Artificial Intelligence, Automated Reasoning, and Symbolic Computation,
%   Joint International Conferences, AISC 2002 and Calculemus 2002, Marseille,
%   France, July 2002, (Jacques Calmet, Belaid Benhamou, Olga Caprotti, Laurent
%   Henocque, and Volker Sorge, Eds.), Berlin: Springer-Verlag, pp. 76-89, 2002.
%   http://orcca.on.ca/TechReports/2000/TR-00-12.html
%
%	Numbers in parentheses below refer to equations in Lawrence, et al. 2012.

% 	The inverse Wright omega function is defined as (Corless & Jeffrey, 2002):
%       W+LOG(W)+2*pi*1i,   -Inf < W < -1
%       -1 +/- pi*1i,       W = 1
%       W+LOG(W),           otherwise

%   WRIGHTOMEGAQIPI is up to three times faster than WRIGHTOMEGAQ, which itself
%   is three to four orders of magnitude faster than the builtin WRIGHTOMEGA for
%   double-precision arrays. Additionally, it has less numeric error, properly
%   evaluates values greater than 2^28, supports single-precision evaluation,
%   and handles NaN inputs.

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-12-12
%   Revision: 1.0, 3-11-13


% Check input
if ~isfloat(Z)
    error('SHCTools:wrightOmegaqipi:InvalidZ',...
          'Input Z must be a floating point array.');
end
isnanZ = isnan(Z(:));
if isreal(Z) && ~(isempty(Z) || all(isnanZ))
    error('SHCTools:wrightOmegaqipi:NonComplexZ',...
         ['Input Z must be complex. NaN is permitted and the empty matrices '...
          'are permitted too.']);
end

x = real(Z);
y = imag(Z);
if ~all(x <= -1)
    error('SHCTools:wrightOmegaqipi:InvalidRealPart',...
          'Real part of input Z must be less than or equal to -1.');
end

% Support for single precision: single(pi) ~= pi
dataType = class(Z);
PI = cast(pi,dataType);
if ~all(abs(y) == PI)
    error('SHCTools:wrightOmegaqipi:InvalidImaginaryPart',...
          'Imaginary part of input Z must equal to Pi or -Pi.');
end

tol = eps(dataType);

if isempty(Z) || all(isnanZ)
    W = Z;
elseif isscalar(Z)
    % Special values
    if x < log(eps(realmin(dataType)))-log(2) && y == PI
        W(1) = 0;                                                   % Z -> -Inf
    elseif x == -1
        W(1) = -1;
    elseif Z == log(1/3)-1/3+PI*1i
        W(1) = -1/3;
    elseif Z == log(2)-2-PI*1i
        W(1) = -2;
    else
        % Regularization, all solutions are real
      	Z = x;                                                          % (26)
        
        % Upper and lower lines of discontinuity
        if x > -2
            % Regions 1 and 2: near z = -1+pi*1i and z = -1-pi*1i
            x = sign(y)*sqrt(-2*(Z+1));                             % (20, 22)
            W(1) = -1+x*(1-x*(1440-x*(120+x*(16+x)))/4320);         % (21, 23)
        elseif y == PI
            % Region 3: series about -Inf
            x = exp(Z);
            W(1) = ((((-125*x-64)*x-36)*x/24-1)*x-1)*x;                 % (24)

            % Series is exact, x < -exp(2)
            if Z < -7.38905609893065
                return;
            end
        else
            % Region 6: negative log series about z = z+pi*1i
            x = log(-Z);
            lzi = x/Z;
            W(1) = Z-x+lzi*(1+lzi*(x/2-1+lzi*((x/3-3/2)*x+1)));         % (28)
        end
        
        % Residual (can be zero)
        r = Z-(W+log(-W));                                              % (14)
        
        if abs(r) > tol
            % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
            W = W*(1+(r/(1+W))*((1+W)*(1+W+(2/3)*r)...
                -r/2)/((1+W)*(1+W+(2/3)*r)-r));                         % (15)
            
            % Test residual
            r = Z-(W+log(-W));                                          % (14)
            
            % Second iterative improvement via FSC method, if needed
            if abs(r) > tol
                W = W*(1+(r/(1+W))*((1+W)*(1+W+(2/3)*r)...
                    -r/2)/((1+W)*(1+W+(2/3)*r)-r));                     % (15)
            end
        end
    end
else
    W = NaN(size(Z),dataType);
    
    % Special values
    W(x < log(eps(realmin(dataType)))-log(2) & y == PI) = 0;        % Z -> -Inf
    W(x == -1) = -1;
    W(Z == log(1/3)-1/3+PI*1i) = -1/3;
    W(Z == log(2)-2-PI*1i) = -2;
    
    % Upper and lower lines of discontinuity
    iz = (isnan(W(:)) & ~isnanZ);
    if any(iz)
        w = W(iz);
        y = y(iz);
        
        % Regularization, all solutions are real
    	Z = x(iz);                                                      % (26)
        
        % Regions 1 and 2: near z = -1+pi*1i and z = -1-pi*1i
        c1 = (Z > -2);
        if any(c1)
            x = sign(y(c1)).*sqrt(-2*(Z(c1)+1));                    % (20, 22)
            w(c1) = -1+x.*(1-x.*(1440-x.*(120+x.*(16+x)))/4320);   	% (21, 23)
        end
        
        % Region 3: series about -Inf
        c = (~c1 & y == PI);
        if any(c)
            x = exp(Z(c));
            w(c) = ((((-125*x-64).*x-36).*x/24-1).*x-1).*x;             % (24)
            
            % Series is exact, x < -exp(2)
            if all(c) && all(Z < -7.38905609893065)
                W(iz) = w;
                return;
            end
        end
        
        % Region 6: negative log series about z = z+pi*1i
        c = (~c1 & y == -PI);
        if any(c)
            zc = Z(c);                                                  % (26)
            x = log(-zc);
            lzi = x./zc;
            w(c) = zc-x+lzi.*(1+lzi.*(x/2-1+lzi.*((x/3-3/2).*x+1)));	% (28)
        end
        
        % Residual (can be zero)
        r = Z-(w+log(-w));                                              % (14)
        
        % FSC-type iteration, N = 3, (Fritsch, Shafer, & Crowley, 1973)
        rc = (abs(r) > tol);
        if any(rc)
            wr = w(rc);
            r = r(rc);
            
            w(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                    -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));                 % (15)
            
            % Test residual
            r = Z-(w+log(-w));                                          % (14)
            
            rc = (abs(r) > tol);
            if any(rc)
                wr = w(rc);
                r = r(rc);
                
                % Second iterative improvement via FSC method, if needed
                w(rc) = wr.*(1+(r./(1+wr)).*((1+wr).*(1+wr+(2/3)*r)...
                        -r/2)./((1+wr).*(1+wr+(2/3)*r)-r));             % (15)
            end
        end
        
        W(iz) = w;
    end
end