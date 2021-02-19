function [boo] = inPoly(x,y,S)
    % Function to use create mask for region inside a polygon S
    n = size(S,1);
    boo = false(size(x));
    for i = 1:n
        boo = max(boo, inpolygon(x,y,S(i).X,S(i).Y));
    end
end