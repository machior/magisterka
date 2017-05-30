function ExtractFeatures(getFeatures, drawPlots, writeData)

    if strcmp(FileName, 'all')
        
        D = dir(['BlindSubCorpus/GENUINE', '/*.sig']);
        num = length(D(not([D.isdir])));
        userNr=1;
        userNrStr = num2str(userNr, '%03i');
        signatureNr=0;
        i=0;
        while i<num
            signatureNr = signatureNr + 1;
            fileName = strcat('BlindSubCorpus/GENUINE/', userNrStr, '_g_', num2str(signatureNr), '.sig');
            if exist(fileName, 'file') ~= 2
                signatureNr = 1;
                userNr = userNr + 1;
                userNrStr = num2str(userNr, '%03i');
                continue;
            end
            ReadSignature(fileName, getFeatures, drawPlots, writeData)
            i = i+1;
        end
        WriteFileName = 'ExtractedFeatures';
        
    end

end