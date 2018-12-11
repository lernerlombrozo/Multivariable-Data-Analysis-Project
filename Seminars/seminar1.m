%% seminar 1 question 3
rng('default') % For reproducibility
x = normrnd(3,5,100,1);
[muHat,sigmaHat] = normfit(x)
bar(x)

rng('default') % For reproducibility
x = betarnd(3,5,100,1);
[p,ci] = betafit(x)
figure
bar(x)

%%Seminar 1 question 6
x=[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1];
y=[-1.9029 -0.2984 0.4047 0.5572 0.9662 2.0312 3.2286 5.722 10.0952];

p = polyfit(x,y,1)
%yfit =  p(1) * x + p(2);
yfit = polyval(p,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (length(y)-1) * var(y);
rsq = 1 - SSresid/SStotal
rsq_adj = 1 - SSresid/SStotal * (length(y)-1)/(length(y)-length(p))
%p2 = polyfit(X,Y,1)

p2 = polyfit(x,y,2)
%yfit2 =  p2(1) * x^2 + p2(2) * x + p2(3);
yfit2 = polyval(p2,x);
yresid2 = y - yfit2;
SSresid2 = sum(yresid2.^2);
SStotal2 = (length(y)-1) * var(y);
rsq2 = 1 - SSresid2/SStotal2
rsq_adj2 = 1 - SSresid2/SStotal2 * (length(y)-1)/(length(y)-length(p2))

p3 = polyfit(x,y,3)
%yfit3 =  p3(1) * x^3 + p3(2) * x^2 + p3(3) * x + p3(4);
yfit3 = polyval(p3,x);
yresid3 = y - yfit3;
SSresid3 = sum(yresid3.^2);
SStotal3 = (length(y)-1) * var(y);
rsq3 = 1 - SSresid3/SStotal3
rsq_adj3 = 1 - SSresid3/SStotal3 * (length(y)-1)/(length(y)-length(p3))

p4 = polyfit(x,y,4)
%yfit4 =  p4(1) * x^4 + p4(2) * x^3 + p4(3) * x^2 + p4(4) * x + p4(5);
yfit4 = polyval(p4,x);
yresid4 = y - yfit4;
SSresid4 = sum(yresid4.^2);
SStotal4 = (length(y)-1) * var(y);
rsq4 = 1 - SSresid4/SStotal4
rsq_adj4 = 1 - SSresid4/SStotal4 * (length(y)-1)/(length(y)-length(p4))