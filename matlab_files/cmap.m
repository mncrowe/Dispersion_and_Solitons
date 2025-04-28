function out = cmap(n,f,c,scale,color,c_lim)
% Creates custom colormap of length n, with scaling function f, centred on c
% - n: map length, default: 256
% - f: scaling function, f:[0, 1] -> [0, 1], used either side of center, default: f(x) = @(x) x
% - c: center value, default: 1/2*(colorbar_max+colorbar_min), cannot be outside range
% - scale: 1; rescales colormap each side of c, 0; don't rescale colorbar, default: 1
% - color: 2N+1 x 3 matrix with color values for bottom, bottommiddle, middle, topmiddle, top
% - c_lim: vector of colormap limits, chooses from active figure if unspecified
%
% ----------------------------------------------------------------------------
% Note: Run after creating figure if c_lim not specified. This script gets
%       the colorbar limits from the active figure.
%
% Note: Default color is [blue white red].
%
% Note: By default c = (c_lim(1)+c_lim(2))/2. The regions either side of c
%       are then scaled according to the function f. If scale = 1 then each
%       side is scaled seperately, if scale = 0 then the scaling of the
%       larger side is applied to both.
% ----------------------------------------------------------------------------

% default arguments
if nargin < 1 || n == 0; n = 256; end
if nargin < 2 || isfloat(f) == 1; f = @(x) x; end
if nargin < 4; scale = 1; end
if nargin < 5 || color == 0; color = [0 0 0.5; 0 0.5 1; 1 1 1; 1 0 0; 0.5 0 0]; end
if nargin < 6; c_lim = get(gca,'CLim'); end
if nargin < 3; c = (c_lim(1)+c_lim(2))/2; end

map_size = size(color); N = (map_size(1)-1)/2;
if N ~= floor(N); error('incorrect color matrix'); end

if c_lim(1) < c && c_lim(2) > c
    % both sides
    
    ratio = abs(c_lim(2)-c)/(c_lim(2)-c_lim(1));
    n_top = ceil(ratio*n);
    n_bot = n - n_top;
    
    % set top_max, bot_max and indexes
    if scale == 1
        top_max = 1;
        bot_max = 1;
    else
        if abs(c_lim(1)-c) > abs(c_lim(2)-c)
            bot_max = 1;
            top_max = abs(c_lim(2)-c)/abs(c_lim(1)-c);
        else
            bot_max = abs(c_lim(1)-c)/abs(c_lim(2)-c);
            top_max = 1;
        end
    end
    top_index = linspace(0,top_max,n_top); top_index = f(top_index);
    bot_index = linspace(0,bot_max,n_bot+1); bot_index = f(bot_index(2:end));
    
    % interpolation sections
    top_N = 1; bot_N = 1;
    for iN = 1:N-1
        top_N = top_N+(top_index>iN/N);
        bot_N = bot_N+(bot_index>iN/N);
    end
    map_top = zeros(n_top,3);
    for i_top = 1:n_top
        for i_col = 1:3
            map_top(i_top,i_col) = interp1([(top_N(i_top)-1)/N top_N(i_top)/N],[color(N+top_N(i_top),i_col) color(N+top_N(i_top)+1,i_col)],top_index(i_top));
        end
    end
    map_bot = zeros(n_bot,3);
    for i_bot = 1:n_bot
        for i_col = 1:3
            map_bot(i_bot,i_col) = interp1([(bot_N(i_bot)-1)/N bot_N(i_bot)/N],[color(N+2-bot_N(i_bot),i_col) color(N+1-bot_N(i_bot),i_col)],bot_index(i_bot));
        end
    end
    out = [map_bot(end:-1:1,:); map_top];
end

if c_lim(1) >= c
    % top only
    index = f(linspace(0,1,n));
    top_N = 1;
    for iN = 1:N-1
        top_N = top_N+(index>iN/N);
    end
    out = zeros(n,3);
    for i_top = 1:n
        for i_col = 1:3
            out(i_top,i_col) = interp1([(top_N(i_top)-1)/N top_N(i_top)/N],[color(N+top_N(i_top),i_col) color(N+top_N(i_top)+1,i_col)],index(i_top));
        end
    end
end

if c_lim(2) <= c
    % bottom only
    index = f(linspace(0,1,n));
    bot_N = 1;
    for iN = 1:N-1
        bot_N = bot_N+(index>iN/N);
    end
    out = zeros(n,3);
    for i_bot = 1:n
        for i_col = 1:3
            out(n+1-i_bot,i_col) = interp1([(bot_N(i_bot)-1)/N bot_N(i_bot)/N],[color(N+2-bot_N(i_bot),i_col) color(N+1-bot_N(i_bot),i_col)],index(i_bot));
        end
    end  
end


end