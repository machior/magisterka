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

function [X Y TStamp Pressure EndPts]=ReadSignature(FileName)

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

    end
    fclose(FID);
end