function ReadSignature(FileName,ShowSig)

    [X Y TStamp Pressure EndPts] = GetParameters(FileName);
    
    [Xmaxes, Ymaxes, Xmins, Ymins] = calculateExtremes(X, Y, TStamp, EndPts);
    [Xmeans, Ymeans, fullArea] = calculateMeans(Xmins, Ymins, Xmaxes, Ymaxes);

    disp( strcat('Center X: ', num2str(mean(Xmeans)), ' Y: ', num2str(mean(Ymeans))) )

%     DrawPlot(Xmaxes, 'X', -Ymaxes, 'Y', EndPts);
%     DrawPlot(Xmins, 'X', -Ymins, 'Y', EndPts);

    %derivatives
    TSteps = diff(TStamp);
    V_X = diff(X)./TSteps;

    %additional
    Slant = ( diff(Y) ./ diff(X) );
    Slant(find(isnan(Slant))) = [];
    Slant = radtodeg( atan(Slant) );
    
    PathTraveled = cumsum( sqrt( diff(X).^2 + diff(Y).^2 ) );
    PathTraveled = [0; PathTraveled];
    PenVelocity = diff(PathTraveled)./TSteps;
    PenAcceleration = diff(PenVelocity)./TSteps(2:end);
    
    if ShowSig

%         DrawPlot(X, 'X', -Y, 'Y', EndPts);
%                DrawPlot(TStamp, 't', X, 'X', EndPts);
%                DrawPlot(TStamp, 't', Y, 'Y', EndPts);
%                DrawPlot(TStamp(2:end), 't', V_X, 'V_X', EndPts);
%                DrawPlot(TStamp(2:end-1), 't', A_X, 'A_X', EndPts);
%                DrawPlot(TStamp, 't', PathTraveled, 'PathTraveled', EndPts);
%                DrawPlot(TStamp(2:end), 't', PenVelocity, 'PenVelocity', EndPts);
%                DrawPlot(TStamp(2:end-1), 't', PenAcceleration, 'PenAcceleration', EndPts);
%                DrawPlot(TStamp, 't', -Y, 'Y', TStampEndPts);

%                DrawPlot(TStamp, 't', Distance, 'Dist', EndPts);
%                DrawPlot(TStamp(2:end), 't', Velocity, 'Vel', EndPts);
%                
%                DrawPlot(TStamp(2:end), 't', V_X, 'V_X', EndPts);
%                DrawPlot(TStamp(2:end), 't', V_Y, 'V_Y', EndPts);
%                DrawPlot(TStamp(2:end), 't', V_P, 'V_P', EndPts);
%                figure, plot(TStamp);

%                 signature duration
        disp( strcat('time: ', num2str(TStamp(end)), ' ms') )

%                 pen-down ratio
        penUpTime = sum( TSteps( TSteps > TSteps(1) ) );
        penDownTime = TStamp(end)-penUpTime;
        disp(strcat('pen-down ratio: ', num2str(penDownTime/TStamp(end))))

%                 horizontal length
        hLength = max(X)-min(X);
        disp( strcat('horizontal length: ', num2str(hLength)) )

%                 aspect ratio
        aRatio = hLength / (max(Y)-min(Y));
        disp( strcat('aspect ratio: ', num2str(aRatio)) )

%                 pen-ups
        penUps = sum(TSteps > TSteps(1)) + 1;
        disp(strcat('number of pen-ups: ', num2str( sum( penUps ) )))

%                 cursiviness
        disp(strcat('cursiviness: ', num2str( hLength / penUps )))

%                 top heaviness
        topHeav = mean(Y) / median(Y);
        disp( strcat('top heaviness: ', num2str(topHeav)) )

%                 Horizontal Dispersion
        horDisp = mean(X) / median(X);
        disp( strcat('top heaviness: ', num2str(horDisp)) )

%                 Curvature
        curvature = sum( pitZ(X, Y) ) / hLength;
        disp( strcat('curvature: ', num2str(curvature)) )
        smoothedX = smoothenPlot(X, 5);
        
        [locExtr, locExtrInd] = getExtremes(smoothedX);
        
        % plot it.
        figure;
        plot(X);
        hold on;
        plot(smoothedX);
        plot(locExtrInd, locExtr, 'o');
        xlswrite('pps.xls',smoothedX);
        
%                 Maximum velocity
        disp( strcat('Maximum velocity: ', num2str( max(PenVelocity) )) )

%                 Average velocity
        disp( strcat('Mean velocity: ', num2str( mean(PenVelocity) )) )

%                 Standard Deviation of the Velocity
        disp( strcat('Standard Deviation of the Velocity: ', num2str( std(PenVelocity) )) )

%                 Average Absolute Acceleration
        disp( strcat('Average Absolute Acceleration: ', num2str( mean(PenAcceleration) )) )

%                 Standard Deviation of the Absolute Acceleration
        disp( strcat('Standard Deviation of the Absolute Acceleration: ', num2str( std(PenAcceleration) )) )

%                 Maximum acceleration
        disp( strcat('Maximum Acceleration: ', num2str( max(PenAcceleration) )) )

%                 Maximum deceleration
        disp( strcat('Maximum Deceleration: ', num2str( min(PenAcceleration) )) )

%                 Handwriting Slant Using All Points
        disp( strcat('Mean Slant: ', num2str( mean(Slant) )) )
        
%                 Horizontal Velocity
        disp( strcat('Horizontal Velocity: ', num2str( mean(V_X) )) )

%                 Mean Pen-Tip Pressure
        disp( strcat('Mean Pen-Tip Pressure: ', num2str( mean(Pressure) )) )

%                 Standard Deviation of Pen-Tip Pressure
        disp( strcat('Standard Deviation of Pen-Tip Pressure: ', num2str( std(Pressure) )) )

%                 Maximum Pen-Tip Pressure
        disp( strcat('Maximum Pen-Tip Pressure: ', num2str( max(Pressure) )) )

%                 Minimum Pen-Tip Pressure
        disp( strcat('Minimum Pen-Tip Pressure: ', num2str( min(Pressure) )) )

%                 Circularity
        disp( strcat('Circularity: ', num2str(fullArea/hLength)) )
        
%                 Area
        disp( strcat('Area: ', num2str(fullArea)) )
    
%                 Middle-Heaviness
        disp( strcat('Middle-Heaviness: ', num2str(100*fullArea/(hLength * (max(Y)-min(Y)))), '%') )
        
%                 Component Physical Spacing
        avSpace = averageSpace(X, Y, TSteps);
        disp( strcat('Component Physical Spacing: ', num2str( avSpace )) )
        
%                 Component Time Spacing
        disp( strcat('Component Time Spacing: ', num2str( penUpTime )) )

    end          

end

function smoothedVector = smoothenPlot(X, windowSize)

       % Construct blurring window.
    windowWidth = int16(windowSize);
    halfWidth = windowWidth / 2;
    gaussFilter = gausswin(5);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.

    % Do the blur.
    smoothedVector = conv(X, gaussFilter);
    smoothedVector = smoothedVector(halfWidth:end-halfWidth+1);
    a = ( smoothedVector(halfWidth) - X(1) ) / ( halfWidth - 1 );
    b = X(1) - a;
    for i=1:halfWidth
        smoothedVector(i)=a*i+b;
    end
    a = ( X(end) - smoothedVector(end-halfWidth) ) / ( halfWidth );
    b = X(end) - a*length(X);
    for i=(length(smoothedVector)-halfWidth):length(smoothedVector)
        smoothedVector(i)=a*i+b;
    end
        
end

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

function avSpace = averageSpace(X, Y, TSteps)

    spacings = (TSteps > TSteps(1));
    shiftedBif = find( spacings>0 );
    spaces = [];
    for i = 1:length(shiftedBif)
        space = sqrt( ( X(shiftedBif(i))-X(shiftedBif(i)+1) )^2 + ( Y(shiftedBif(i))-Y(shiftedBif(i)+1) )^2 );
        spaces = [spaces; space];
    end
    if length(spaces)
        avSpace = sum(spaces) / length(spaces);
    else
        avSpace = 0;
    end

end

function Z = pitZ(X, Y)

    Z = sqrt( diff(X).^2 + diff(Y).^2 );

end


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