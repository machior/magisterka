function [Xmeans, Ymeans, fullArea] = calculateMeans(Xmins, Ymins, Xmaxes, Ymaxes)

    Xmeans = [];
    Ymeans = [];
    fullArea = 0;
    for i = 1 : length(Xmins)
        if isnan(Xmins(i)) || isnan(Xmaxes(i))
            continue
        end

        if ( Ymaxes(i) - Ymins(i) ) < inf && Ymaxes(i) ~= Ymins(i)
           fullArea = fullArea + ( Ymaxes(i) - Ymins(i) );

           Xmeans = [ Xmeans; Xmins(i) ];
           Ymeans = [ Ymeans; (Ymaxes(i) + Ymins(i))/2 ];
        end

    end

end