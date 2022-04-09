function barCtrs = plotFBAxSAData(dataAll, SEMs, doLabels)

dylims=[0 2.05];

ct1ToPlot=1:size(dataAll,1);
ct2ToPlot=1:size(dataAll,2);

p.groupByColor=true;
ncolors=length(ct2ToPlot);

hues = ones(1,ncolors)*0.78-0.75*p.groupByColor;
sats = zeros(1,ncolors); %now grascale
vals = linspace(.3,1,ncolors);

if ncolors==1
    vals=.5;
end

hsvs = [hues' sats' vals'];

myCols=hsv2rgb(hsvs);

myFontName='Arial';
axisFontSize=18;
axisLabFontSize=18;
legendFontSize=10;

bwid=.9;
bsep=.25;

bsepgroup=0.9;

bi=0; barCtrs=[];

xtickvals = [];
groupctrs = [];
bctr = bsep;

cueValTitles={'Spatial Invalid','Spatial neutral','Spatial valid'};
colCueValLegend={'Invalid','Neutral','Valid'};



for ct1 = ct1ToPlot
    thisgroupctrs = [];
    for ct2 = ct2ToPlot
        
        if p.groupByColor
            locv=ct2;
            colv=ct1;
        else
            locv=ct1;
            colv=ct2;
        end
        
        hold on;
        bi=bi+1;
        
        if ct1>ct1ToPlot(1) && ct2 == ct2ToPlot(1)
            bctr=bctr+bwid+bsepgroup;
        else
            bctr=bctr+bwid+bsep;
        end
        
        
        bx=bctr+[-bwid/2 bwid/2];
        by=[0 dataAll(locv, colv)];
        
        vertx=[bx; bx];
        verty=[by fliplr(by)];
        
        colri=ct2;
        
        
        cvHandles(colri)=fill(vertx(:), verty(:),myCols(colri,:));
        
        plot([bctr bctr],by(end)+[-1 1]*SEMs(locv,colv),'k-');
        
        barCtrs(bi)=bctr;
        thisgroupctrs=[thisgroupctrs bctr];
        
        
    end
    groupctrs=[groupctrs mean(thisgroupctrs)];
    
    xtickvals = [xtickvals mean(thisgroupctrs)];
    
end


set(gca,'XTick',xtickvals);
if ~p.groupByColor
    set(gca,'XTickLabel',locCueValLegend);
    xlabel('Location cue validity','FontSize',axisLabFontSize);
else
    set(gca,'XTickLabel',colCueValLegend);
    xlabel('Color cue validity','FontSize',axisLabFontSize);
end
xlim([0 max(bx)+bsep]);

if dylims(1)<0
    plot([0 max(bx)+bsep], [0 0],'k-');
end

ylim(dylims);

set(gca,'YTick',dylims(1):.4:dylims(2));

ylabtext='d''';
legendLoc='NorthWest';

if doLabels
    ylabel(ylabtext,'FontSize',axisLabFontSize);
    legend(cvHandles,cueValTitles,'Location',legendLoc,'FontSize',legendFontSize);
end

set(gca, 'FontSize', axisFontSize);
set(gca,'FontName', myFontName);

