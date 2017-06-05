function TrainNetwork(FileName, getFeatures, drawPlots, writeData)

    if strcmp(FileName, 'all')
        
        filePath = 'ExtractedFeatures/VisualSubCorpus/FORGERY/';
        D = dir([filePath, '/*.sig']);
        num = length(D(not([D.isdir])));
        userNr=1;
        userNrStr = num2str(userNr, '%03i');
        signatureNr=0;
        i=0;
        csvFileContent = [];
        while i<num
            signatureNr = signatureNr + 1;
            fileName = strcat(filePath, userNrStr, '_f_', num2str(signatureNr), '.sig');
            disp(strcat(num2str(i+1), '/', num2str(num), ' ', fileName))
            if exist(fileName, 'file') ~= 2
                signatureNr = 0;
                writeFileName = strcat( 'ExtractedFeatures/', filePath, userNrStr );
                csvwrite(writeFileName,csvFileContent);
                csvFileContent = [];
                userNr = userNr + 1;
                userNrStr = num2str(userNr, '%03i');
                continue;
            end
%             writeFileName = strcat( 'ExtractedFeatures/BlindSubCorpus/GENUINE/', userNrStr, '_g_', num2str(signatureNr) );
            csvFileContent = [ csvFileContent, ReadSignature(fileName, getFeatures, drawPlots, writeData) ];
            i = i+1;
        end
        
        writeFileName = strcat( 'ExtractedFeatures/', filePath, userNrStr );
        csvwrite(writeFileName,csvFileContent);
%         WriteFileName = strcat( 'ExtractedFeatures/BlindSubCorpus/GENUINE', userNrStr, '_g_', num2str(signatureNr) );
        
    end

end