function yy = smooth(y, span)
    yy = y;
    l = length(y);

    for i = 1 : l
        if i < span
            d = i;
        else
            d = span;
        end

        w = d - 1;
        p2 = floor(w / 2);

        if i > (l - p2)
           p2 = l - i; 
        end

        p1 = w - p2;

        yy(i) = sum(y(i - p1 : i + p2)) / d;
    end
end