function [locExtr, locExtrInd] = getExtremes(X)

    [locMaxes, locMaxesInd] = findpeaks(X);
    [locMins, locMinsInd] = findpeaks(-X);
    locExtr = [locMaxes; -locMins];
    locExtrInd = [locMaxesInd; locMinsInd];
    
    n = length(locExtr);
    while (n > 0)
        % Iterate through x
        nnew = 0;
        for i = 2:n
            % Swap elements in wrong order
            if (locExtr(i) < locExtr(i - 1))
                locExtr = swap(locExtr,i,i - 1);
                locExtrInd = swap(locExtrInd,i,i - 1);
                nnew = i;
            end
        end
        n = nnew;
    end

end

function x = swap(x,i,j)
    % Swap x(i) and x(j)
    % Note: In practice, x xhould be passed by reference

    val = x(i);
    x(i) = x(j);
    x(j) = val;

end