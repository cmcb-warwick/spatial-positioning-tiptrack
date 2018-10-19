function [] = PairedNucleotides()

% PURPOSE:
%   Computes distribution of nucleotides as a function of the distance from the microtubule tip
%   assuming first order kinetics for GTP hydrolysis (k1) and phosphate release (k2).
%
% INPUT:  
%   no input required to run the function
%   However, desired rates can be changed in lines 43-50.
%
% OUTPUT:
%   A grid of 9 subplots with calculated distributions for GTP, GDP/Pi and GDP as function from 
%   microtubule tip for three different combinations of reaction rates (left column). Distribution of pairs of
%   nucleotides as number of dimers (middle column) and normalised and colvolved distributions (right column).
%   This will be saved in current folder unless other destination specified in line 68.
%
% Copyright: Straube lab, University of Warwick, July 2016
% shared under New BSD license as specified below
%
% Code was tested under Matlab release 2012b




FH = 1;
figure(FH); clf; hold on
set(FH,'Units','Centimeters')
set(FH,'Position',[10,10,16,15])

for i = 1:9
    subplot(3,3,i)
    l = (mod(i-1,3)*5)+2;
    b = (3-ceil((i)/3))*4.7;
    set(gca,'Units','Centimeters')
%     set(gca,'OuterPosition',[l,b,6,4.5]);
    if l == 10
        set(gca,'Position',[l+1.2,b+1,3.2,3.3]);
    else
        set(gca,'Position',[l+1  ,b+1,3.4,3.3]);
    end
end
global k1 k2                % rates are given per tubulin layer
k1 = 0.13;                  % this equals 0.325 s-1 at 20nm/s assembly
k2 = 0.052;                 % this equals 0.130 s-1 at 20nm/s assembly
EndPosition_v2(1,[3,3,1])
k1 = 0.26;
k2 = 0.026;
EndPosition_v2(1,[3,3,4])
k1 = 0.065;
k2 = 0.104;
EndPosition_v2(1,[3,3,7])

for i = 9:-1:1
    subplot(3,3,i)
    l = (mod(i-1,3)*5);
    b = (3-ceil((i)/3))*4.7;
    set(gca,'Units','Centimeters','LineWidth',0.75)
    if l == 10
        set(gca,'OuterPosition',[l b 5 4.7],'Position',[l+1.6,b+1.4,3.4,3.1]);
    else
        set(gca,'OuterPosition',[l b 5 4.7],'Position',[l+1.4,b+1.4,3.6,3.1]);
    end
end

set(FH,'Units','Centimeters')
set(FH,'Position',[10,10,19,15])

savename = [pwd '/PairedNucleotides.eps'];
saveas(FH, savename,'eps2c')


% text('string','k$_1$ = 0.025 k$_2$ = 0.1','FontWeight','bold','FontSize',9,'FontName','Helvetica', ...
%     'VerticalAlignment','middle','HorizontalAlignment','center', ...
%     'Units','centimeters','Rotation',90,'Position',[-1,1.65],)

end

function [] = EndPosition_v2(FH,SP)

FS = 9;

options = odeset('RelTol',1e-4,'AbsTol',[1e-4 1e-4 1e-5]);
[T,Y] = ode45(@rigid,[0 200],[13 0 0],options);
L = length(T);
Z = zeros(L,6);
for i = 1:L
p1 = Y(i,1);
p2 = Y(i,2);
p3 = Y(i,3);
C  = StateDist_v2(13,p1,p2,p3);
Z(i,:) = C;
end

figure(FH); subplot(SP(1),SP(2),SP(3)); hold on
t = 0:0.1:200;
X = spline(T,Y(:,1),t);
plot(t*8,X,'Color',[1,0,0],'LineWidth',1)
X = spline(T,Y(:,2),t);
plot(t*8,X,'Color',[0,128/255,0],'LineWidth',1)
X = spline(T,Y(:,3),t);
plot(t*8,X,'Color',[0,0,1],'LineWidth',1)
set(gca,'TickDir','Out')
ylim([0 14])
xlim([0 800])

% legend('GTP','GDP+Pi','GDP','Location','NorthEast')
set(gcf,'color','white')
set(gca,'FontSize',FS,'FontName','Helvetica')
ylabel('number of dimers','VerticalAlignment','baseline','Margin',3)
xlabel('distance from MT tip (nm)','VerticalAlignment','middle','Margin',3)

set(gca,'TickDir','Out','YTick',(0:2:14)','YTickLabel',(0:2:14)','TickLength',[0.0175 0.01])
set(gca,'TickDir','Out','XTick',(0:200:800)','XTickLabel',(0:200:800)','TickLength',[0.0175 0.01])


cc = [255/255,  0/255,  0/255; ...  
	  0/255,  0/255,  0/255; ...
      255/255,  0/255,255/255; ...  
	  0/255,128/255,  0/255; ...
	  0/255,  1/255,255/255];

figure(FH); subplot(SP(1),SP(2),SP(3)+1); hold on
hold on;
for i = 1:5
X = spline(T,Z(:,i),t);
plot(t*8,X,'Color',cc(i,:),'LineWidth',1)
end
% plot(T,Z(:,1),'r-',T,Z(:,2),'r-.',T,Z(:,3),'r.')
% plot(T,Z(:,4),'b-',T,Z(:,5),'b-.',T,Z(:,6),'g.')
X = spline(T,Z(:,6),t);
h = line(t*8,X);
set(h,'Color','blue','LineWidth',1)
X = spline(T,sum(Z(:,1:3),2),t);
h = line(t*8,X);
set(h,'Color',cc(1,:),'LineWidth',1)
X = spline(T,sum(Z(:,[2,4,5]),2),t);
h = line(t*8,X);
set(h,'Color',cc(4,:),'LineWidth',1)
% plot(T,sum(Z(:,1:3),2),'k-')
% plot(T,sum(Z(:,[2,4,5]),2),'m-')
% legend('GTP','GTP & GDP+Pi','GTP & GDP','GDP+Pi','GDP+Pi & GDP','GDP','Location','NorthEast')
xlim([0 800])
ylim([0 14])
set(gca,'TickDir','Out')
set(gcf,'color','white')
set(gca,'FontName','Helvetica')
set(gca,'FontSize',FS)
ylabel('number of dimers','Margin',3,'VerticalAlignment','baseline')
xlabel('distance from MT tip (nm)','Margin',3,'VerticalAlignment','middle')
set(gca,'TickDir','Out','YTick',(0:2:14)','YTickLabel',(0:2:14)','TickLength',[0.0175 0.01])
set(gca,'TickDir','Out','XTick',(0:200:800)','XTickLabel',(0:200:800)','TickLength',[0.0175 0.01])



t = (0:0.1:200)';
for i = 1:6
z(:,i) = spline(T,Z(:,i),t);
end


sig = 130;
T = [(-100:0.1:-0.1)'; t]*8;
Z = [zeros(1000,6); z];
L = length(T);

figure(FH); subplot(SP(1),SP(2),SP(3)+2); hold on
hold on;
cc = [255/255,  0/255,  0/255; ...  
	  0/255,  0/255,  0/255; ...
      255/255,  0/255,255/255; ...  
	  0/255,128/255,  0/255; ...
	  0/255,  0/255,255/255];
for j = 1:5
A = zeros(L,1);
B = sum(Z(:,j),2);
for i = 1:L
A = A + B(i)*normpdf(T,T(i),sig);
end
A = A./max(A);

plot(T,A,'Color',cc(j,:),'LineWidth',1);
end


% legend('GTP','GTP & GDP+Pi','GTP & GDP','GDP+Pi','GDP+Pi & GDP','Location','NorthEast')

xlim([-300 800])
ylim([0 1.05])
set(gcf,'color','white')
set(gca,'FontSize',FS,'FontName','Helvetica')
ylabel('relative intensity (au)','Margin',100,'VerticalAlignment','baseline')
xlabel('distance from MT tip (nm)','Margin',3,'VerticalAlignment','middle')
set(gca,'TickDir','Out','YTick',(0:0.2:1)','YTickLabel',{'0.0';'0.2';'0.4';'0.6';'0.8';'1.0'},'TickLength',[0.0175 0.01])
	set(gca,'TickDir','Out','XTick',(-200:200:800)','XTickLabel',(-200:200:800)','TickLength',[0.0175 0.01])
	
	
end	
	
	
	function dy = rigid(t,y)
	
	global k1 k2
	
	% k1 = 0.05;
	% k2 = 0.05;
	
	dy = zeros(3,1);    % a column vector
	dy(1) = -k1*y(1);
	dy(2) = k1*y(1) - k2*y(2);
	dy(3) = k2*y(2);
	
end


function [C] = StateDist_v2(n, p1, p2, p3)

t1 = p1/(p1+p2+p3);
t2 = p2/(p2+p3+p1);
t3 = p3/(p2+p3+p1);

C = zeros(1,6);

C(1) =   t1*t1*n;
C(2) = 2*t1*t2*n;
C(3) = 2*t1*t3*n;
C(4) =   t2*t2*n;
C(5) = 2*t2*t3*n;
C(6) =   t3*t3*n;

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