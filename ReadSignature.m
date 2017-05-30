
%   Parameters:
%   FileName: file name of the signature.
%   ShowSig : if set to non 0 value displays signature.
function ReadSignature(FileName, getFeatures, drawPlots, writeData)

    [X Y TStamp Pressure EndPts] = GetParameters(FileName);

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
    
    if getFeatures
        
        [Xmaxes, Ymaxes, Xmins, Ymins] = calculateExtremes(X, Y, TStamp, EndPts);
        [Xmeans, Ymeans, fullArea] = calculateMeans(Xmins, Ymins, Xmaxes, Ymaxes);
        
%                 center of signature
        xCenter = mean(Xmeans);
        yCenter = mean(Ymeans);

%                 signature duration
        signatureDuration = TStamp(end);
        
%                 Component Time Spacing
        penUpTime = sum( TSteps( TSteps > TSteps(1) ) );
        
%                 pen-down ratio
        penDownTime = TStamp(end)-penUpTime;
        penDownRatio = penDownTime/TStamp(end);

%                 horizontal length
        hLength = max(X)-min(X);

%                 aspect ratio
        aRatio = hLength / (max(Y)-min(Y));

%                 pen-ups
        penUps = sum(TSteps > TSteps(1)) + 1;

%                 cursiviness
        cursiviness = hLength / penUps;

%                 top heaviness
        topHeav = mean(Y) / median(Y);

%                 Horizontal Dispersion
        horDisp = mean(X) / median(X);

%                 Curvature
        curvature = sum( pitZ(X, Y) ) / hLength;
        
        % Strokes
        smoothedX = smoothenPlot(X, 5);
        smoothedY = smoothenPlot(Y, 5);
        [XlocExtr, XlocExtrInd] = getExtremes(smoothedX);
        [YlocExtr, YlocExtrInd] = getExtremes(smoothedY);
        hStrokes = length(XlocExtrInd) + 1;
        vStrokes = length(YlocExtrInd) + 1;
%         xlswrite('pps.xls',smoothedX);
        
%                 Maximum velocity
        maxVelocity = max(PenVelocity);

%                 Average velocity
        meanVelocity = mean(PenVelocity);

%                 Standard Deviation of the Velocity
        stdPenVelocity = std(PenVelocity);

%                 Average Absolute Acceleration
        meanPenAcceleration = mean(PenAcceleration);

%                 Standard Deviation of the Absolute Acceleration
        stdPenAcceleration = std(PenAcceleration);

%                 Maximum acceleration
        maxPenAcceleration = max(PenAcceleration);

%                 Maximum deceleration
        minPenAcceleration = min(PenAcceleration);

%                 Handwriting Slant Using All Points
        meanSlant = mean(Slant);
        
%                 Horizontal Velocity
        meanV_X = mean(V_X);

%                 Mean Pen-Tip Pressure
        meanPressure = mean(Pressure);

%                 Standard Deviation of Pen-Tip Pressure
        stdPressure = std(Pressure);

%                 Maximum Pen-Tip Pressure
        maxPressure = max(Pressure);

%                 Minimum Pen-Tip Pressure
        minPressure = min(Pressure);

%                 Circularity
        circularity = fullArea/hLength;
        
%                 Area
        fullArea;
    
%                 Middle-Heaviness
        midHeaviness = fullArea/(hLength * (max(Y)-min(Y)));
        
%                 Component Physical Spacing
        avSpace = averageSpace(X, Y, TSteps);
        
    end
        
        
    if drawPlots
        
%         DrawPlot(Xmaxes, 'X', -Ymaxes, 'Y', EndPts);
%         DrawPlot(Xmins, 'X', -Ymins, 'Y', EndPts);
        DrawPlot(X, 'X', -Y, 'Y', EndPts);
        DrawPlot(TStamp, 't', X, 'X', EndPts);
        DrawPlot(TStamp, 't', Y, 'Y', EndPts);
        DrawPlot(TStamp(2:end), 't', V_X, 'V_X', EndPts);
        DrawPlot(TStamp, 't', PathTraveled, 'PathTraveled', EndPts);
        DrawPlot(TStamp(2:end), 't', PenVelocity, 'PenVelocity', EndPts);
        DrawPlot(TStamp(2:end-1), 't', PenAcceleration, 'PenAcceleration', EndPts);
        DrawPlot(TStamp(2:end), 't', V_X, 'V_X', EndPts);
        extremesPlot(smoothedX, XlocExtr, XlocExtrInd)
        extremesPlot(smoothedY, YlocExtr, YlocExtrInd)
        
    end
        
    if writeData
        
        disp( strcat('Center X: ', num2str(xCenter), ' Y: ', num2str(yCenter)) )
        disp( strcat('time: ', num2str(signatureDuration), ' ms') )
        disp( strcat('Component Time Spacing: ', num2str( penUpTime )) )
        disp( strcat('pen-down ratio: ', num2str(penDownRatio)) )
        disp( strcat('horizontal length: ', num2str(hLength)) )
        disp( strcat('aspect ratio: ', num2str(aRatio)) )
        disp( strcat('number of pen-ups: ', num2str( penUps )) )
        disp( strcat('cursiviness: ', num2str( cursiviness )) )
        disp( strcat('top heaviness: ', num2str(topHeav)) )
        disp( strcat('top heaviness: ', num2str(horDisp)) )
        disp( strcat('curvature: ', num2str(curvature)) )
        disp( strcat('Number of Horizontal Strokes: ', num2str( max(hStrokes) )) )
        disp( strcat('Number of Vertical Strokes: ', num2str( max(vStrokes) )) )
        disp( strcat('Maximum velocity: ', num2str( maxVelocity )) )
        disp( strcat('Mean velocity: ', num2str( meanVelocity )) )
        disp( strcat('Standard Deviation of the Velocity: ', num2str( stdPenVelocity )) )
        disp( strcat('Average Absolute Acceleration: ', num2str( meanPenAcceleration )) )
        disp( strcat('Standard Deviation of the Absolute Acceleration: ', num2str( stdPenAcceleration )) )
        disp( strcat('Maximum Acceleration: ', num2str( maxPenAcceleration )) )
        disp( strcat('Maximum Deceleration: ', num2str( minPenAcceleration )) )
        disp( strcat('Mean Slant: ', num2str( meanSlant )) )
        disp( strcat('Horizontal Velocity: ', num2str( meanV_X )) )
        disp( strcat('Mean Pen-Tip Pressure: ', num2str( meanPressure )) )
        disp( strcat('Standard Deviation of Pen-Tip Pressure: ', num2str( stdPressure )) )
        disp( strcat('Maximum Pen-Tip Pressure: ', num2str( maxPressure )) )
        disp( strcat('Minimum Pen-Tip Pressure: ', num2str( minPressure )) )
        disp( strcat('Circularity: ', num2str(circularity)) )
        disp( strcat('Area: ', num2str(fullArea)) )
        disp( strcat('Middle-Heaviness: ', num2str(100*midHeaviness), '%') )
        disp( strcat('Component Physical Spacing: ', num2str( avSpace )) )

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

function extremesPlot(mainPlot, locExtr, locExtrInd)

    % plot it.
    figure;
    plot(mainPlot);
    hold on;
    plot(locExtrInd, locExtr, 'o');

end