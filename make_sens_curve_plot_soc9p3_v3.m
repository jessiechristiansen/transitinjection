
% This script takes the outputs of the SOC 9.3 pipeline run on injected targets
% Injections are on-target, FGK dwarfs, and have at least three transits falling on valid cadences
% The results are split into P<100 and P>100 day bins, since the behaviour is still a function of period
% The gamma CDF is fit separately to the P<100 and P>100 and overplotted on the data

load dr25alldata_lists

h = figure;
xbar = 0.75:0.5:20.75;
xbarp = 0.75:0.5:50.75;

% show theoretical erf
MES = 0:0.1:25;
y = (erf(MES-7.1))/2 + 0.5;
h1 = plot(MES,y,'r-.','LineWidth',3);
hold on;
ylim([0,1.01]);
ylim([0,1.005]);

fractionWindowedMatchedAlllt100 = zeros(101,1);
numMesesAll3lt100 = zeros(101,1);
detectedMesesAll3lt100 = zeros(101,1);
binomialErrorlt100 = zeros(101,1);
fitFlaglt100 = false(101,1);

for i = 1:40
    numMesesAll3lt50(i) = sum(fgkMes > ((i*0.5)) & fgkMes <=(((i+1)*0.5)) & fgkPeriods>400 & fgkThreeTransits);% & nooffsetIndx);
    detectedMesesAll3lt50(i) = sum(fgkMatch(fgkMes > ((i*0.5)) & fgkMes <=(((i+1)*0.5)) & fgkPeriods>400 & fgkThreeTransits));% & nooffsetIndx));
    thisFraction = detectedMesesAll3lt50(i)/numMesesAll3lt50(i);
    fractionWindowedMatchedAlllt100(i) = thisFraction;
    if(thisFraction<1)
        binomialErrorlt100(i) = sqrt((thisFraction*(1-thisFraction))/numMesesAll3lt50(i));
    else
        tmpFraction = (detectedMesesAll3lt50(i)-1)/numMesesAll3lt50(i);
        binomialErrorlt100(i) = sqrt((tmpFraction*(1-tmpFraction))/numMesesAll3lt50(i));
    end
    if(binomialErrorlt100(i)>0)
        fitFlaglt100(i) = true;
    end
    display(['For ' num2str(i*0.5) '-' num2str((i+1)*0.5) ': ' num2str(detectedMesesAll3lt50(i)) ' of ' num2str(numMesesAll3lt50(i)) ' found'])
    display(['For ' num2str(i*0.5) '-' num2str((i+1)*0.5) ': binomial error = ' num2str(binomialErrorlt100(i))])
end

%plot lt100 solution
uplim = 32; % only fitting through to MES=16, too noisy between 16 and 20
xdata = xbar(1:uplim);
xdata = xdata(fitFlaglt100(1:uplim));
ydata = fractionWindowedMatchedAlllt100(1:uplim);
ydata = ydata(fitFlaglt100(1:uplim));
ydata(isnan(ydata)) = 0;

%show gamma function fit
gammodel=@(x,alp,bet,c) c.*gamcdf(x,alp,bet);
sumresid=@(vec) sum((gammodel(xdata,vec(1),vec(2),vec(3))-ydata').^2);
[resgam,fval,exflg,outstruct]=fminsearch(sumresid,[35.0,0.2,0.95]);
h3 = plot(xbar(1:41),gammodel(xbar(1:41),resgam(1),resgam(2),resgam(3)),'LineWidth',3);

display(['resgam(1) = ' num2str(resgam(1)) ', resgam(2) = ' num2str(resgam(2)) ', resgam(3) = ' num2str(resgam(3))]);
xtest = 0:0.01:20;
gamtest = gammodel(xtest,resgam(1),resgam(2),resgam(3));
fiftypercent = min(xtest(gamtest>0.5));
fiftypercentString = ['Detection crosses 50% = ' num2str(fiftypercent)];
text(1,0.65, ['4000-7000K, logg>4.0, P<100 days'])
text(1,0.60,fiftypercentString)
text(1,0.55,['Plateau at 16 sigma: ' num2str(resgam(3))])
h6 = plot([0.5,0.9],[0.65,0.65],'b','LineWidth',3);

% long period
fractionWindowedMatchedAllgt100 = zeros(101,1);
numMesesAll3gt100 = zeros(101,1);
detectedMesesAll3gt100 = zeros(101,1);
binomialErrorgt100 = zeros(101,1);
fitFlaggt100 = false(101,1);

for i = 1:40
    numMesesAll3gt50(i) = sum(fgkMes > ((i*0.5)) & fgkMes <=(((i+1)*0.5)) & fgkPeriods>=100 & fgkThreeTransits);% & nooffsetIndx);
    detectedMesesAll3gt50(i) = sum(fgkMatch(fgkMes > ((i*0.5)) & fgkMes <=(((i+1)*0.5)) & fgkPeriods>=100 & fgkThreeTransits));% & nooffsetIndx));
    thisFraction = detectedMesesAll3gt50(i)/numMesesAll3gt50(i);
    fractionWindowedMatchedAllgt100(i) = thisFraction;
    if(thisFraction<1)
        binomialErrorgt100(i) = sqrt((thisFraction*(1-thisFraction))/numMesesAll3lt50(i));
    else
        tmpFraction = (detectedMesesAll3lt50(i)-1)/numMesesAll3lt50(i);
        binomialErrorgt100(i) = sqrt((tmpFraction*(1-tmpFraction))/numMesesAll3lt50(i));
    end
    if(binomialErrorgt100(i)>0)
        fitFlaggt100(i) = true;
    end
    display(['For ' num2str(i*0.5) '-' num2str((i+1)*0.5) ': ' num2str(detectedMesesAll3gt50(i)) ' of ' num2str(numMesesAll3gt50(i)) ' found'])
    display(['For ' num2str(i*0.5) '-' num2str((i+1)*0.5) ': binomial error = ' num2str(binomialErrorgt100(i))])
end

%plot gt100 solution
uplim = 32; % only fitting through to MES=16, too noisy between 16 and 20
xdata = xbar(1:uplim);
xdata = xdata(fitFlaggt100(1:uplim));
ydata = fractionWindowedMatchedAllgt100(1:uplim);
ydata = ydata(fitFlaggt100(1:uplim));
ydata(isnan(ydata)) = 0;

%show gamma function fit
gammodel=@(x,alp,bet,c) c.*gamcdf(x,alp,bet);
sumresid=@(vec) sum((gammodel(xdata,vec(1),vec(2),vec(3))-ydata').^2);
[resgam,fval,exflg,outstruct]=fminsearch(sumresid,[35.0,0.2,0.95]);
h4 = plot(xbar(1:41),gammodel(xbar(1:41),resgam(1),resgam(2),resgam(3)),'r','LineWidth',3);

display(['resgam(1) = ' num2str(resgam(1)) ', resgam(2) = ' num2str(resgam(2)) ', resgam(3) = ' num2str(resgam(3))]);
xtest = 0:0.01:20;
gamtest = gammodel(xtest,resgam(1),resgam(2),resgam(3));
fiftypercent = min(xtest(gamtest>0.5));
fiftypercentString = ['Detection crosses 50% = ' num2str(fiftypercent)];
text(1,0.50, ['4000-7000K, logg>4.0, P>100 days'])
text(1,0.45,fiftypercentString)
text(1,0.40,['Plateau at 16 sigma: ' num2str(resgam(3))])
h5 = plot([0.5,0.9],[0.50,0.50],'r','LineWidth',3);


bh = bar(xbarp(1:41),fractionWindowedMatchedAlllt100(1:41));
ch = get(bh,'child');
set(ch,'facea',.3)
hold on;
bh = bar(xbarp(1:41),fractionWindowedMatchedAllgt100(1:41),'r');
ch = get(bh,'child');
set(ch,'facea',.3)
hold on;
fractionWindowedMatchedAlllt100 = fractionWindowedMatchedAlllt100';
fractionWindowedMatchedAllgt100 = fractionWindowedMatchedAllgt100';
binomialErrorlt100 = binomialErrorlt100';
binomialErrorgt100 = binomialErrorgt100';
h2 = plot([xbar(1:41); xbar(1:41)],[fractionWindowedMatchedAlllt100(1:41)-binomialErrorlt100(1:41); ...
    fractionWindowedMatchedAlllt100(1:41)+binomialErrorlt100(1:41)], 'b', 'LineWidth', 2);
h7 = plot([xbar(1:41); xbar(1:41)],[fractionWindowedMatchedAllgt100(1:41)-binomialErrorgt100(1:41); ...
    fractionWindowedMatchedAllgt100(1:41)+binomialErrorgt100(1:41)], 'r', 'LineWidth', 2);
xlim([0,20]);
%xlim([0,25]);
fontname = 'Helvetica';
set(0,'DefaultAxesFontName','CMU Serif Roman');
set(0,'defaulttextfontname',fontname);
set(gca,'FontSize',16);
xlabel(['Expected MES ($\sigma$)'],'interpreter','latex');
ylabel('Fraction detected','FontSize',16,'interpreter','latex');

%save sensitivityhistogram.mat xbar fractionWindowedMatchedAll fgkMes nooffsetIndx fgkMatches detectedMesesAll3 numMesesAll3

set(h,'PaperUnits','inches','PaperPosition',[0 0 10 7]);

xex = [7.1,7.1];
yex = [0,1];
plot(xex, yex, 'k--','LineWidth',2);
text(6,0.77,'7.1$\sigma$','FontSize',16,'interpreter','latex');

legend([h1 h3 h4], {'7.1$\sigma$ error function','$\Gamma$ CDF $<$100 days','$\Gamma$ CDF $>$100 days'},'Location','NorthWest', 'FontSize', 18,'interpreter','latex')


print(h, '-dpng', 'highres_senscurvwithfits.png', '-r300');
print(h, '-dpsc', 'highres_senscurvwithfits.eps');

