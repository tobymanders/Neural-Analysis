noiseGauge(1,length(smoothDiff)*10) = 0;
for y = 1:(length(borders)-1);
    start   = index(borders(y),1);
    stop    = index(borders(y+1)-1,1);
    sectionLength = (stop-start);
    subData = smoothDiff(start:stop);
    if median(subData)>(meanVal+s)
        for p = (10*start):(10*stop)
            noiseGauge(p)=1;
            hippoStruct(m).denoised(p)=0;
        end
    end
end

for z = 2:length(smoothDiff);
    noiseEdge(z)=(noiseGauge(z)~=noiseGauge(z-1));
end
edges = find(noiseEdge==1);

% Plot noise removal regions
fig=figure;
hax=axes;
obj = plot (smoothDiff, 'LinestopWidth', 2);
for num = 1:length(edges);
    SP = edges(num);
    line([SP SP],get(hax,'YLim'),'Color',[1 0 0], 'LineWidth', 0.25, 'LineStyle', '--');
end
uistack(obj);
hline1 = refline (0, meanVal+2*s);
hline1.Color = 'r';
hline1.LineStyle = '--';
hline2 = refline (0, meanVal-2*s);
hline2.Color = 'r';
hline2.LineStyle = '--';
