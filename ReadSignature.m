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

function [X Y TStamp Pressure EndPts]=ReadSignature(FileName,ShowSig)

        FID=fopen(FileName,'r'); 
        
      if FID==-1
            error(strcat('could not open signature file: ', FileName));
      else
            fgets(FID); %discard info line
            fgets(FID,10); %discard
            NumOfPoints=fscanf(FID,'%d',1);
                
            %read sig. data
            SigData=fscanf(FID,'%d%d%d%d%d',NumOfPoints*5);
            
            %parse data
            i=5*[0:NumOfPoints-1];
            X=SigData(1+i); 
            Y=SigData(2+i);
            TStamp=SigData(3+i);
            Pressure=SigData(4+i);
            EndPts=SigData(5+i);
            
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
            
           DrawPlot(Xi, 'X', -Yi, 'Y', EndPts);
           
           Xmaxes = [];
           Ymaxes = [];
           Xmins = [];
           Ymins = [];
           for i = min(Xi) : max(Xi)
               
%                indexes = Xi==Xi(i);
               indexes =  Xi==i;
               colX = Xi(indexes);
               colY = Yi(indexes);
               
               maxx = max(colY);
               minn = min(colY);
               
               Xmaxes = [Xmaxes; i];
               Ymaxes = [Ymaxes; maxx];
               Xmins = [Xmins; i];
               Ymins = [Ymins; minn];
               
           end
           
           area = 0;
           for i = min(Xmins) : max(Xmins)
               if Xmins(i)==nan
                   continue
               end
               index = find(Xmaxes==Xmins(i))
               if index
                   area = area + (Ymaxes(index) - Ymins(i));
                   
               end
               
           end
           
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
                ups = [strfind( TSteps' > TSteps(1), 1 ); 1]
                curvSum = sum( pitZ( X(1:ups(1)), Y(1:ups(1)) ) ) / (max(Y(1:ups(1))) - min(Y(1:ups(1))));
                
                for i = 1 : ( length(ups) )
                    
                    if i ~= length(ups)
%                         disp( TStamp(ups(i)+1 : ups(i+1)) )
                        curvSum = curvSum + sum( pitZ( X(ups(i)+1 : ups(i+1)), Y(ups(i)+1 : ups(i+1)) ) ) / (max(Y(ups(i)+1 : ups(i+1))) - min(Y(ups(i)+1 : ups(i+1))));
                    else
%                         disp( TStamp(ups(i)+1 : end) )
                        curvSum = curvSum + sum( sqrt(diff( X(ups(i)+1 : end) ).^2 + diff( Y(ups(i)+1 : end) ).^2) ) / (max(Y(ups(i)+1 : end)) - min(Y(ups(i)+1 : end)));
                        curvSum = curvSum + sum( pitZ( X(ups(i)+1 : end), Y(ups(i)+1 : end ) ) ) / (max(Y(ups(i)+1 : end)) - min(Y(ups(i)+1 : end)));
                    end
                    
                end
                
                disp( strcat('average curvature per stroke: ', num2str( curvSum / (length(ups)+1) )) )
                
%                 Number of strokes
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
                
%                 Handwriting Slant Using “Long Stroke” End-points:
%                 Handwriting Slant Using All Points of “Long Strokes”
%                 Handwriting Slant Through Regression of “Long Strokes”:
%                 Handwriting Slant Using Cai and Liu Technique
%                 Handwriting Slant Based on Vertical Overlap
%                 Stroke Concavity
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
                
%                 Degree of Parallelism
%                 Baseline Consistency
%                 Ascender-line Consistency
%                 Circularity
%                 Area
%                 Middle-Heaviness
%                 Component Physical Spacing
%                 Component Time Spacing
                disp( strcat('Component Time Spacing: ', num2str( penUpTime )) )

           end
        end
           fclose(FID);            

end

function DrawPlot(x, xBadge, y, yBadge, EndPts)

    eInd=find(EndPts==1); 
    Ytmp=y; Xtmp=x;
    Ytmp(eInd)=nan; Xtmp(eInd)=nan;
    figure,plot(Xtmp,Ytmp,'.'); hold on; 
    plot(Xtmp,Ytmp);
    xlabel(xBadge); ylabel(yBadge);
    warning off;
    title(strcat(yBadge, '(', xBadge, ')'));

end

function Z = pitZ(X, Y)

    Z = sqrt( diff(X).^2 + diff(Y).^2 );

end