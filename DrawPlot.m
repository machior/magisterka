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