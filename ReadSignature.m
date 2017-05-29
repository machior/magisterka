%This function is provided as it is, without any warranty.
%   Parameters:
%   FileName: file name of the signature.
%   ShowSig : if set to non 0 value displays signature.

%   Output:
%   X           : x coordinates
%   Y           : y coordinates
%   TStamps     : point acquisition time in msec (relative to 1'st point) 
%   Pressure    : self explanatory
%   EndPts      : 1 indicates end of a signature segment, 0 other wise

function ReadSignature(FileName,ShowSig)

    [X Y TStamp Pressure EndPts]=GetParameters(FileName);
            
    [Xmaxes, Ymaxes, Xmins, Ymins] = calculateExtremes(X, Y, TStamp, EndPts);
    [Xmeans, Ymeans, fullArea] = calculateMeans(Xmins, Ymins, Xmaxes, Ymaxes);

    disp( strcat('Center X: ', num2str(mean(Xmeans)), ' Y: ', num2str(mean(Ymeans))) )

    DrawPlot(Xmaxes, 'X', -Ymaxes, 'Y', EndPts);
    DrawPlot(Xmins, 'X', -Ymins, 'Y', EndPts);

    %derivatives
    TSteps = diff(TStamp);
    V_X = diff(X)./TSteps;
    V_Y = diff(Y)./TSteps;
    V_P = diff(Pressure)./TSteps;

    A_X = diff(V_X)./TSteps(2:end);
    A_Y = diff(V_Y)./TSteps(2:end);
    A_P = diff(V_P)./TSteps(2:end);

    %additional
    Distance = sqrt( X.^2 + Y.^2 );
    Slant = ( diff(Y) ./ diff(X) );
    Slant(find(isnan(Slant))) = [];
    Slant = radtodeg( atan(Slant) );
%             diff( Y(isnan(Slant)))
%             diff( X(isnan(Slant)))
%             diff( Y(isnan(Slant))) ./ diff(X(isnan(Slant)))
%             Slant(isnan(Slant))
    PathTraveled = cumsum( sqrt( diff(X).^2 + diff(Y).^2 ) );
    PathTraveled = [0; PathTraveled];
    PenVelocity = diff(PathTraveled)./TSteps;
    PenAcceleration = diff(PenVelocity)./TSteps(2:end);
    Velocity = sqrt( V_X.^2 + V_Y.^2 );
%             Acceleration = sqrt( A_X.^2 + A_Y.^2 );
%             A_X(1:10)
%             A_Y(1:10)

    if ShowSig

        DrawPlot(X, 'X', -Y, 'Y', EndPts);
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

        TSteps = diff(TStamp);

%                 signature duration
%                 disp(TStamp(end))
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
%                 disp(sum(TSteps > TSteps(1)))
        disp(strcat('number of pen-ups: ', num2str( sum( penUps ) )))

%                 cursivity

%                 cursiviness
        disp(strcat('cursiviness: ', num2str( hLength / penUps )))

%                 top heaviness
        topHeav = mean(Y) / median(Y);
        disp( strcat('top heaviness: ', num2str(topHeav)) )

%                 Horizontal Dispersion
        horDisp = mean(X) / median(X);
        disp( strcat('top heaviness: ', num2str(horDisp)) )

%                 Curvature
%                 curvature = sum( sqrt(diff(X).^2 + diff(Y).^2) ) / hLength;
        curvature = sum( pitZ(X, Y) ) / hLength;
        disp( strcat('curvature: ', num2str(curvature)) )

%                 Average Curvature per Stroke
%         ups = [strfind( TSteps' > TSteps(1), 1 ); 1];
%         curvSum = sum( pitZ( X(1:ups(1)), Y(1:ups(1)) ) ) / (max(Y(1:ups(1))) - min(Y(1:ups(1))));
% 
%         for i = 1 : ( length(ups) )
% 
%             if i ~= length(ups)
% %                         disp( TStamp(ups(i)+1 : ups(i+1)) )
%                 curvSum = curvSum + sum( pitZ( X(ups(i)+1 : ups(i+1)), Y(ups(i)+1 : ups(i+1)) ) ) / (max(Y(ups(i)+1 : ups(i+1))) - min(Y(ups(i)+1 : ups(i+1))));
%             else
% %                         disp( TStamp(ups(i)+1 : end) )
%                 curvSum = curvSum + sum( sqrt(diff( X(ups(i)+1 : end) ).^2 + diff( Y(ups(i)+1 : end) ).^2) ) / (max(Y(ups(i)+1 : end)) - min(Y(ups(i)+1 : end)));
%                 curvSum = curvSum + sum( pitZ( X(ups(i)+1 : end), Y(ups(i)+1 : end ) ) ) / (max(Y(ups(i)+1 : end)) - min(Y(ups(i)+1 : end)));
%             end
% 
%         end
% 
%         disp( strcat('average curvature per stroke: ', num2str( curvSum / (length(ups)+1) )) )

%                 Number of strokes
                % Construct blurring window.
        windowWidth = int16(5);
        halfWidth = windowWidth / 2
        gaussFilter = gausswin(5)
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

        [locMaxes, locMaxesInd] = findpeaks(smoothedVector);
        [locMins, locMinsInd] = findpeaks(-smoothedVector);
        
        [locExtr, locExtrInd] = getExtremes(X);
        
        % plot it.
        figure;
        plot(X);
        hold on;
        plot(smoothedVector);
        plot(locExtrInd, locExtr, 'o');
        
        
%                 Mean Ascender Height:
%                 Mean Descender Depth:
%                 Maximum height

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

function [locExtr, locExtrInd] = getExtremes(X)

    % Construct blurring window.
    windowWidth = int16(5);
    halfWidth = windowWidth / 2
    gaussFilter = gausswin(5)
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

    [locMaxes, locMaxesInd] = findpeaks(smoothedVector);
    [locMins, locMinsInd] = findpeaks(-smoothedVector);
    locExtr = [locMaxes; -locMins];
    locExtrInd = [locMaxesInd; locMinsInd];

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