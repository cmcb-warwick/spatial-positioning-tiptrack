function [A,GG,RR] = DCK_Curve_Analysis(varargin)

% PURPOSE:
%   This summarises and plots all the data of velocity and relative EB peak
%   positions from dual-colour experiments. You need to run DCK_select.m and
%   DCK_Kymo_Run.m first.

% INPUT:  
%   provide 'path/DCA_Kymo_Run.mat' as input argument. 

% OUTPUT:
%   saves data files with superaveraged intersity profiles for red and
%   green channel and with peak differences and microtubule growth velocity
%   for each phase analysed
%   as well as pdf files with averaged curves and a stairs histogram of peak differences
%   for all combinations of EB proteins analysed. Additional figures with speed and peak difference
%   histograms are generated and can be saved individually.


% Copyright: Straube lab, University of Warwick, July 2016
%   shared under New BSD license as specified below
%
% Code was tested under Matlab release 2012b


%--------------------------------------------------
% Load and Save Directory - Asking
%--------------------------------------------------

if isempty(varargin)
    load_loc = uigetdir(pwd,'Pick Load Directory');
else
    load_loc = varargin{1};
end

MT = load(load_loc);
MT = MT.MT;

j = 0;
m = 0;
Pre = 'pol';
% Pre1 = '1';
L = length(MT);

[savefile,~,~] = fileparts(load_loc);
savefile = [savefile '/DCA_diff_velo.txt'];
savefile = fopen(savefile,'w+');

CG = zeros(L,310);
CR = zeros(L,310);
for i = 1:L                                             % go through all data
    
    
    if ~strcmp(Pre,MT(i).EB_Proteins)                                   % check whether we are looking at data from a new EB protein combination
        m
        if j > 0                                                        % if yes, and also not the first MT, produce averaged curves, save and plot
            t = -60:249;
            t = t*8;
            CG(isnan(CG(:,1)),:) = [];
            sG = std(CG,0,1)./sqrt(m);
            mG = mean(CG,1);
            sR = std(CR,0,1)./sqrt(m);
            mR = mean(CR,1);
            
            i % just for testing, can delete later
            MT(i).EB_Proteins           % just for testing, can delete later
            
            SC(j).curveg = mG;
            SC(j).errorg = sG;
            SC(j).curver = mR;
            SC(j).errorr = sR;

            figure(31);
            subplot(2,3,j); cla;
            fe = errorbar(t,mG,sG);
            set(fe,'Color',[0 1 0]);
            hold on
            fe = errorbar(t,mR,sR);
            set(fe,'Color',[1 0 0]);
            title(Pre); 
            xlabel('Distance (nm)')
            ylabel('Intensity (au)')
            xlim([-500 1000])
            
            
        end
        j = j+1;                                                        % increase counter for different experimental conditions and reset m
        m = 1;
        
        clear CG
        clear CR
        Desc{j} = MT(i).EB_Proteins;
    else
        m = m+1;
    end
    
        Pre = MT(i).EB_Proteins;
    
    
    R = MT(i).curve_C1;
    G = MT(i).curve_C2;
    
    figure(1); clf;
    plot(R./max(R),'r')
    hold on
    plot(G./max(G),'g')
    
    if length(R) < 20; continue; end
%     Rmax = find(R >= 0.3*max(R),1,'first');
%     Gmax = find(G >= 0.5*max(G),1,'first');
    Rmax = find(R == max(R),1,'first');
    if isempty(Rmax)
        continue
    end
    
    if Rmax < 61 || length(R) - Rmax < 21 
        continue
    end
    if isnan(G)
        continue
    end
    Gmax = find(G(1:Rmax+20) == max(G(1:Rmax+20)),1,'first');
    
   diff = Gmax-Rmax;         % calculate difference between peak positions between green and red curve      
    
    if isnan(diff) 
        continue
    end
    
    if isempty(diff)
        continue
    end
%     
%     if abs(diff) > 200
%         pause()
%     end
    
    
    
    xlim([Rmax-30 Rmax+100])
    
    
%     if Rmax < 51
%         Rmax = 51;
%     end
    
    mm = min(length(R),Rmax+249);
    l = mm-Rmax+61;
    
    CG(m,1:l) = G(Rmax-60:mm)./max(G(Gmax:Gmax+5));                     %write green curve relative to red peak position
    CR(m,1:l) = R(Rmax-60:mm)./max(R(Rmax:Rmax+5));                     %write red curve relative to red peak position
  
%     CG(m,1:l) = G(Rmax-50:mm);
%     CR(m,1:l) = R(Rmax-50:mm);


    A(j).diff(m) = diff;
    A(j).vHis(m) = MT(i).Velo;
    
    Number1 = find(MT(i).Green_Parent == '_',1,'last');
    Number1 = MT(i).Green_Parent(Number1+1:end);
    
    Number2 = find(MT(i).Green_File == '_',1,'last');
    Number2 = MT(i).Green_File(Number2+1:Number2+3);
    
    fprintf(savefile,'%s \t %s \t %s \t %s \t %02f \t %02f\n', ...
        MT(i).EB_Proteins,MT(i).Stack,Number1,Number2,diff*8, MT(i).Velo*16);
    
    
end
m
j

CG(isnan(CG(:,1)),:) = [];
t = -60:249;
t = t*8;
sG = std(CG,0,1)./sqrt(m);
mG = mean(CG,1);
sR = std(CR,0,1)./sqrt(m);
mR = mean(CR,1);

SC(j).curveg = mG;
SC(j).errorg = sG;
SC(j).curver = mR;
SC(j).errorr = sR;

figure(31);
subplot(2,3,j); cla;
fe = errorbar(t,mG,sG);
set(fe,'Color',[0 1 0]);
hold on
fe = errorbar(t,mR,sR);
set(fe,'Color',[1 0 0]);
title(Pre); 
xlabel('Distance (nm)')
ylabel('Intensity (au)')
xlim([-500 1000])

clear CG
clear CR


% GG = mG;
% RR = mR;

fclose(savefile);

figure(23);                 %plot histograms with peakdifferences for all conditions
clf;
for i = 1:j
    subplot(2,ceil(j/2),i);
    D = A(i).diff;
    X = -40:40;
    Y = hist(D,X)./length(D);
    x = gaussfit(X,Y);
    bar(X,Y,'k');
    hold on
    F = myfun(x,X);
    plot(X,F,'r', 'LineWidth', 0.75)
    set(gca,'FontSize',14)
    title(Desc{i}); 
    xlabel('Distance (nm)');
    xlim([-20 20])
    ylim([0 0.2])
    set(gca,'TickDir','Out','XTick',-18:6:18,'XTickLabel',(-18:6:18)*8); 
    ylabel('Relative Frequency');
    str = ['Mean = ' num2str(x(2)*8,'%5.2f') ];
    text(-19,0.18,str)
    str = ['Std = ' num2str(x(3)*8,'%5.2f')];
    text(-19,0.165,str)
    str = ['Median = ' num2str(median(D)*8,'%5.2f')];
    text(-19,0.15,str)
    str = ['n = ' num2str(length(D),'%5.2f')];
    text(-19,0.135,str)
    set(gca,'TickDir','out')
    
    
    
    % write file with average curves for all conditions
    
    [savefile,~,~] = fileparts(load_loc);
savefile = [savefile '/DCA_curve' Desc{i} '.txt'];
savefile = fopen(savefile,'w+');

fprintf(savefile,'%s \t %s \t %s \t %s \n','Mean_Green','StEr_Green','Mean_Red','StEr_Red');

for k = 1:310
    fprintf(savefile,'%08.7f \t %08.7f \t %08.7f \t %08.7f \n', ...
        SC(i).curveg(k),SC(i).errorg(k),SC(i).curver(k),SC(i).errorr(k));
end
fclose(savefile);
    
    
    
end

% Plot stairs with differences

H=figure(24);
clf; hold on
X = -320:8:320;
cc = hsv(j);

set(H,'PaperUnit','centimeters','PaperSize',[6 6])
a = get(H,'PaperSize');
set(H,'PaperPositionMode','Manual','PaperPosition',[0 0 a(1) a(2)])

for i = 1:j
    
    D = A(i).diff.*8;
    
    Y = hist(D,X)./length(D);
    stairs(X,Y,'Color',cc(i,:))
    fe = plot([(median(D)) (median(D))],[0 0.01]);
    set(fe,'Color',cc(i,:),'LineWidth',0.75);
    str = ['Median = ' num2str(median(D),'%5.1f') ];
    text(50,(0.19-0.02*i),str,'Color',cc(i,:))
end

FtSz = 10;
Font = 'Arial';
set(gca,'FontName',Font,'FontSize',FtSz)

ylabel('relative frequency')
xlabel('distance (nm)')
xlim([-96 144])
ylim([0 0.2])
set(gca,'TickDir','Out','XTick',-96:48:144,'XTickLabel',(-96:48:144));
set(gca,'TickDir','Out','YTick',0:0.05:0.2,'YTickLabel',(0:0.05:0.2));

[savefile,~,~] = fileparts(load_loc);
savefile = [savefile '/Diff_Stairs.pdf'];
saveas(H,savefile,'pdf')


% Plot curves with microtubule growth velocity distributions

figure(25); clf;
for i = 1:j
    subplot(2,ceil(j/2),i)
    X = -2:0.1:8;
    Y = hist(A(i).vHis,X)./length(A(i).vHis);
    x = gaussfit(X,Y);
    bar(X*16,Y,'k');
    hold on
    F = myfun(x,X);
    plot(X*16,F,'r', 'LineWidth', 0.75)
    set(gca,'FontSize',14)
    title(Desc{i}); 
    xlim([20 120])
    xlabel('Velocity (nm/s)');
    ylim([0 0.2])
    ylabel('Relative Frequency');
    set(gca,'TickDir','Out','XTick',20:10:120,'XTickLabel',20:10:120); 
    str = ['Mean = ' num2str(x(2)*16,'%5.2f') ];
    text(22,0.18,str)
    str = ['Std = ' num2str(x(3)*16,'%5.2f')];
    text(22,0.165,str)
    str = ['Median = ' num2str(median(A(i).vHis)*16,'%5.2f')];
    text(22,0.15,str)
    str = ['n = ' num2str(length(A(i).vHis),'%5.2f')];
    text(22,0.135,str)
    set(gca,'TickDir','out')
   
    
end



   
 % Plot scatters of peak difference versus microtubule growth velocity

figure(26); clf;
for i = 1:j
    subplot(2,ceil(j/2),i)
    X = A(i).vHis .* 16;
    Y = A(i).diff .* 8;
    hold on
    
    R=corrcoef(X,Y)
    R_squared=R(2)^2;
    scatter(X,Y,10,cc(i,:),'filled')
    
    % plot(X,Y,'.','Color',cc(i,:),'MarkerSize',5)
    set(gca,'FontSize',14)
    title(Desc{i}); 
    xlim([20 120])
    xlabel('Velocity (nm/s)');
    ylim([-48 96])
    ylabel('Peak difference (nm)');
    
    text(22, 90, ['R^2 = ' num2str(R_squared)])
    
    set(gca,'TickDir','Out','XTick',20:10:120,'XTickLabel',20:10:120); 
    set(gca,'TickDir','Out','YTick',-96:48:144,'YTickLabel',(-96:48:144));
   
    
end  
   
   
 CCC = [1,0,0;...
       0,1,0;...
       0,0,1;...
       0,1,1;...
       1,1,0];
  
   
   
% size(SC)
   
figure(40); clf; hold on
t = -60:249;
t = t*8;
for i = 1:j
    fe = errorbar(t,SC(i).curveg,SC(i).errorg);
    set(fe,'Color',CCC(i,:));
end

FtSz = 10;
Font = 'Arial';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting up to 5 curves pretty and one summary curve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

H = figure(41); clf;               

set(H,'PaperUnit','centimeters','PaperSize',[18 12])
a = get(H,'PaperSize');
set(H,'PaperPositionMode','Manual','PaperPosition',[0 0 a(1) a(2)])

for i = 1:j
    
    subplot(2,3,i);
    hold on

% set(H,'PaperUnit','centimeters','Color',[1,1,1])
% axes('Units','centimeters','OuterPosition',[0 5 6 5],'Position',[1.6 6.2 4 3.7])

mG = SC(i).curveg./max(SC(i).curveg);
sG = SC(i).errorg./max(SC(i).curveg);
mR = SC(i).curver./max(SC(i).curver);
sR = SC(i).errorr./max(SC(i).curver);

t = -60:249;
t = t*8;

% size(t)      % for testing
% size(mG)
% size(mR)

fe = errorbar(t,mR,sR,'r');
set(fe,'Color',[1 0 0],'LineWidth',0.75);
fe.CapSize = 2;

peakR = find(mR == max(mR),1,'last');
plot([t(peakR) t(peakR)],[0 max(mR)],'r','LineWidth',0.75);

fe = errorbar(t,mG,sG,'g');
set(fe,'Color',[0 1 0],'LineWidth',0.75);
fe.CapSize = 2;

peakG = find(mG == max(mG),1,'last');
plot([t(peakG) t(peakG)],[0 max(mG)],'g','LineWidth',0.75);

xlabel('relative position (nm)','FontSize',FtSz,'FontName',Font)
ylabel('relative intensity','FontSize',FtSz,'FontName',Font)
xlim([-400 800])
ylim([0 1.1])
set(gca,'box','off','LineWidth',0.75)
set(gca,'TickDir','Out','YTick',[(0:0.1:1.1)'],'YTickLabel',{'0';'';'0.2';'';'0.4';'';'0.6';'';'0.8';'';'1.0';''},'TickLength',[0.0175 0.01])
set(gca,'TickDir','Out','XTick',[(-400:200:800)'],'XTickLabel',{'-400';'';'0';'';'400';'';'800';},'TickLength',[0.0175 0.01])

title(Desc{i})
title([Desc{i} ' Dist: ' num2str(t(peakG),'%05.1f') ' N: ' num2str(size(A(i).diff,2),'%05d')],'FontSize',8)




end

subplot(2,3,6); hold on
for i = 1:j
    fe = errorbar(t,(SC(i).curveg./max(SC(i).curveg)),(SC(i).errorg./max(SC(i).curveg)));
    set(fe,'Color',CCC(i,:),'LineWidth',0.75);
    fe.CapSize = 2;
    
    peakG = find(SC(i).curveg == max(SC(i).curveg),1,'last');
    fe = plot([t(peakG) t(peakG)],[0 1.0]);
    set(fe,'Color',CCC(i,:),'LineWidth',0.75);
end
xlabel('relative position (nm)','FontSize',FtSz,'FontName',Font)
ylabel('relative intensity (au)','FontSize',FtSz,'FontName',Font)
xlim([-400 800])
ylim([0 1.1])
set(gca,'box','off','LineWidth',0.75)
set(gca,'TickDir','Out','YTick',[(0:0.1:1.1)'],'YTickLabel',{'0';'';'0.2';'';'0.4';'';'0.6';'';'0.8';'';'1.0';''},'TickLength',[0.0175 0.01])
set(gca,'TickDir','Out','XTick',[(-400:200:800)'],'XTickLabel',{'-400';'';'0';'';'400';'';'800';},'TickLength',[0.0175 0.01])



[savefile,~,~] = fileparts(load_loc);
savefile = [savefile '/AveragedCurves.pdf'];
saveas(H,savefile,'pdf')

end



function [x] = gaussfit(X,Y)

x0(1) = max(Y);    
x0(2) = X(find(Y == max(Y),1,'first'));
x0(3) = 4;

lb = [0.5*x0(1),x0(2)-4,0];
ub = [1.5*x0(1),x0(2)+4,10];

options = optimset('Display','off','TolFun',1e-8);

x = lsqcurvefit(@myfun,x0,X,Y,lb,ub,options);
end

function F = myfun(x,X)
    F = x(1)*exp( -((X-x(2)).^2)./(2*x(3)^2));
end
    
function [estimates r rr] = DecayFit(xdata, ydata)
    
x1 = find(ydata > 0.5,1,'first');
x2 = find(ydata > 0.5*(1-ydata(end)),1,'last');

x0 = [1   x1    10  1   x2    10];
lb = [0.9 x1-1  1   0.5 x1+10 1];
ub = [2.0 x2-10 50  1.0 x2+1  50];

options   = optimset('Display','off','TolFun',1e-6);
[estimates r rr] = lsqcurvefit(@DecayFtn, x0,xdata,ydata, lb, ub,options);
end


function [Y] = DecayFtn(x,xdata)
Y = x(1).*0.5.*(1+erf((xdata-x(2))./(x(3).*sqrt(2)))) - x(4).*0.5.*(1+erf((xdata-x(5))./(x(6).*sqrt(2))));
end

   % Copyright (c) 2016, Straube lab, University of Warwick
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in the
%       documentation and/or other materials provided with the distribution.
%     * Neither the name of the organization nor the
%       names of its contributors may be used to endorse or promote products
%       derived from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
   % LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    