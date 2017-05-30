
function [Xmaxes, Ymaxes, Xmins, Ymins] = calculateExtremes(X, Y, TStamp, EndPts)

    [Xi, Yi] = fillInterpolation(X, Y, TStamp);

%     DrawPlot(Xi, 'X', -Yi, 'Y', EndPts);

    [Xmaxes, Ymaxes, Xmins, Ymins] = getExtremes(Xi, Yi);

end

function [Xi, Yi] = fillInterpolation(X, Y, TStamp)

    Xi = [X(1)];
    Yi = [Y(1)];
    for i = 1 : (length(X) - 1)
        a = (Y(i) - Y(i+1)) / (X(i) - X(i+1));
        b = (Y(i)*X(i+1) - X(i)*Y(i+1)) / (X(i+1) - X(i));

        if TStamp(i+1) - TStamp(i) > TStamp(2) - TStamp(1)
            continue
        end

        if X(i) > X(i+1)
            for xn = X(i) : -1 : X(i+1)
                Xi = [Xi; xn];
                Yi = [Yi; a*xn+b];
            end
        else
            for xn = X(i) : X(i+1)
                Xi = [Xi; xn];
                Yi = [Yi; a*xn+b];
            end
        end                    

    end

end

function [Xmaxes, Ymaxes, Xmins, Ymins] = getExtremes(Xi, Yi)

    Xmaxes = [];
    Ymaxes = [];
    Xmins = [];
    Ymins = [];
    for i = min(Xi) : max(Xi)

        indexes = Xi==i;
        colY = Yi(indexes);

        maxx = max(colY);
        minn = min(colY);
        if isempty(maxx)
            continue
        end

        Xmaxes = [Xmaxes; i];
        Ymaxes = [Ymaxes; maxx];
        Xmins = [Xmins; i];
        Ymins = [Ymins; minn];

    end

end