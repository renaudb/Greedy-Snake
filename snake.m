function [x, y, it] = snake(Im, x, y, alpha, beta, gamma, thresh, wsize, tol, sigma)
% Given an image and a list of points, computes the locally optimal snake.
%
% Input: 
%   Im        The input image.
%   [x, y]    The snake starting points.
%
% Output:
%   [x, y]    The snake ending points.
%
% Parameters:
%   alpha     The continuity energy weighting term.
%   beta      The curvature energy weighting term.
%   gamma     The image energy weighting term.
%   thresh    A tuple containing the angle and image energy thresholds for
%             relaxation.
%   wsize     The window size.
%   tol       The maximum number of points that can change at any
%             iteration.
%   sigma     The gaussian filter standard deviation.
%

    % Validate number of arguments.
    if nargin < 3
        error('Not enough input arguments.');
    end
    
    % Extract sizes.
    [n, m] = size(x);

    % Fill default arguments.
	if nargin < 4, alpha = 1.0; end
	if nargin < 5, beta = 1.0; end
	if nargin < 6, gamma = 1.2; end
	if nargin < 7, thresh = [60.0, 0.15]; end
	if nargin < 8, wsize = 5; end
	if nargin < 9, tol = ceil(0.5 * n); end
    if nargin < 10, sigma = sqrt(2); end
    
    % Edge filter image.
    Imf = edge(Im, sigma);
        
    % Convert threshold from degree.
    thresh(1) = (2 * sin(deg2rad(thresh(1))/2))^2;
        
    % Initialize coefficients vectors.
    alphas = alpha * ones(n, 1);
    betas  = beta  * ones(n, 1);
    gammas = gamma * ones(n, 1);

    % Count number of iterations.
    it = 0;
    
    % Track number of points moved.
    moved = Inf;
    
    while moved > tol
        % Increment number of iteration.
        it = it + 1;
        
        % Compute average distance between points.
        d = davg(x, y);

        % Reset move count.
        moved = 0;
                
        % Update each point.
        for i = 1:n
            % Compute next and previous points.
            xp = circshift(x, 1); yp = circshift(y, 1);
            xn = circshift(x, -1); yn = circshift(y, -1);

            % Compute energy terms for point neighborhood.
            Cont = econtinuity(x(i), y(i), xp(i), yp(i), wsize, d);
            Curv = ecurvature(x(i), y(i), xp(i), yp(i), xn(i), yn(i), wsize);
            Img  = eimage(x(i), y(i), wsize, Imf);
        
            % Compute total energy for each points in point neighborhood.
            En = alphas(i) * Cont + betas(i) * Curv + gammas(i) * Img;
                    
            % Find point in neighborhood of minal energy.
            [r, c] = find(En == min(min(En)));
            xnew = x(i) + (c(1) - (wsize+1)/2);
            ynew = y(i) + (r(1) - (wsize+1)/2);
            
            % Check if the point has moved.
            if (xnew ~= x(i) || ynew ~= y(i))
                moved = moved + 1;
            end
            
            % Update the values in the point vectors.
            x(i) = xnew;
            y(i) = ynew;
        end
                
        % Relax corners.
        betas(corners(x, y, Imf, thresh)) = 0;
    end
end

function [davg] = davg(x, y)
    % Compute pairwise distance in x and y.
	dx = x - circshift(x, 1);
	dy = y - circshift(y, 1);
        
    % Average the pairwise total distance.
	davg = mean(sqrt(dx.^2 + dy.^2));
end

function [C] = econtinuity(x, y, xp, yp, w, d)
    % Compute half window size.	
    h = (w - 1) / 2;
        
    % Compute X and Y neighborhood matrices.
	Mx = ones(1, w)' * [-h:1:h] + x;
    My = [-h:1:h]' * ones(1, w) + y;
        
    % Compute position change in x and y.
    Dx = Mx - xp;
    Dy = My - yp;
        
    % Compute squared error from average distance.
	C = (d - sqrt(Dx.^2 + Dy.^2)).^2;

    % Normalize the result.
    M = max(max(C));
    if M ~= 0
        C = C / M;
    end
end

function [C] = ecurvature(x, y, xp, yp, xn, yn, w)
    % Compute half window size.
	h = (w - 1) / 2;
        
    % Compute X and Y neighborhood matrices.
    Nx = ones(1, w)' * [-h:1:h] + x;
    Ny = [-h:1:h]' * ones(1, w) + y;

    % Compute curvature change in x and y.
    Cx = xp - 2 * Nx + xn;
    Cy = yp - 2 * Ny + yn;
        
    % Compute curvature norm.
    C = Cx.^2 + Cy.^2;
    
    % Normalize the result.
    M = max(max(C));
    if M ~= 0
        C = C / M;
    end
end

function [E] = eimage(x, y, w, F)
    % Compute half window size.
    h = (w - 1) / 2;
        
    % Select neighborhood of point in filtered image.
    Fpad = padarray(F, [h, h], 0);                        
    E = Fpad(y:y+2*h, x:x+2*h);

    % Normalize the result.
    m = min(min(E));
    M = max(max(E));    
    if M - m ~= 0
        E = (E - m) / (M - m);
    end
    
    % Return negative energies.
    E = -E;
end

function [t] = corners(x, y, Imf, thresh)
	% Align previous and next points.
    xp = circshift(x, 1); yp = circshift(y, 1);
    xn = circshift(x, -1); yn = circshift(y, -1);

    % Compute x and y values for previous and next vectors.
    dx1 = x - xp; dy1 = y - yp;
    dx2 = xn - x; dy2 = yn - y;
        
    % Compute previous and next vectors norm.
    n1 = sqrt(dx1.^2 + dy1.^2);
    n2 = sqrt(dx2.^2 + dy2.^2);
        
    % Compute normalized curvature at point.
    c = (dx1 ./ n1 - dx2 ./ n2).^2 + (dy1 ./ n1 - dy2 ./ n2).^2;
    
    % Find local maxima.
    t = (c > circshift(c, 1)) .* (c > circshift(c, -1));
        
    % Discard points with curvature smaller than threshold.
    t = t .* (c > thresh(1));
    
    % Discard points with edge strenght smaller than threshold.
    t = t .* (Imf(sub2ind(size(Imf), y, x)) > thresh(2));
    
    % Return a logical array.
    t = logical(t);
end

function [Imf] = edge(Im, sigma)
    % Create separate difference of gaussian functions.
    [f, df] = dofgauss(sigma);
    
    % Filter the image.
	Imfx = conv2(df, f, Im, 'same');
    Imfy = conv2(f, df, Im, 'same');
    
    % Return the modulus.
    Imf = sqrt(Imfx.^2 + Imfy.^2);        
end

function [g, dg] = dofgauss(sigma)
    % Create a window of size equal to three standard deviations.
    wsize = ceil(sigma * 3) * 2 + 1;

    % Define functions range.
    r = -(wsize-1)/2:1:(wsize-1)/2;
        
    % Compute functions.
    c = 1 / (sqrt(2 * pi) * sigma);
    e = exp(-(r.^2) / (2 * sigma^2));
    g  = c .* e;
    dg = c .* (-r / sigma^2) .* e;
end