function [ ] = LineScans (PathName, scale)

% PURPOSE:
%   analyses and averages line scan data. the algorithm expects a folder (PathName) with
%   multiples of 3 txt files, with one for each set of EB1, EB2 and EB3.

% INPUT:  
%   PathName = give path to folder that contains linescan data
%   scale = give pixel resolution, default is 0.069061 for 100x objective on
%   spinning disk.

% OUTPUT:
%   saves plots and data of averaged linescans and S.E.M. once normalised to 
%   EB1 only, once normalised to EB3 only and once normalised to average of EB1 and EB3 positions.
%   All normalisation is to first half-maximal point.

% Copyright: Straube lab, University of Warwick, July 2016
%   shared under New BSD license as specified below
%
% Code was tested under Matlab release 2012b



if nargin<2
    scale = 0.069061;
end

datafiles = dir([PathName '*linescan*.txt'])    % Lists the linescan txt files in the Sample directory in a vector

lines = length(datafiles)/3;

dataset = cell(14,1);
uploc1 = zeros(lines,1);
uploc3 = zeros(lines,1);
ylength1 = zeros(lines,1);
ylength3 = zeros(lines,1);
uploc2 = zeros(lines,1);
ylength2 = zeros(lines,1);


for k = 1:lines;

filename1 = [PathName datafiles(3*k-2).name];
filename2 = [PathName datafiles(3*k-1).name];
filename3 = [PathName datafiles(3*k).name];

fid  = fopen(filename1);
xh=fgets(fid);
data1 = textscan(fid,'%f\t%f');
fclose(fid);

ylength = length(data1{1});  % Gets number of rows

data = zeros(ylength,4);
data(:,1)=[0:scale:(scale*(ylength-1))];
data(:,2) = data1{2};

fid  = fopen(filename2);
xh=fgets(fid);
data1 = textscan(fid,'%f\t%f');
fclose(fid);

data(:,3) = data1{2};

fid  = fopen(filename3);
xh=fgets(fid);
data1 = textscan(fid,'%f\t%f');
fclose(fid);

data(:,4)=data1{2};

[minvalue] = min(data);

data(:,2)=data(:,2)-minvalue(2);
data(:,3)=data(:,3)-minvalue(3);
data(:,4)=data(:,4)-minvalue(4);

uploc1(k) = find(data(:,2) > (max(data(:,2))/2),1,'first');   %finds first value half of the maximum EB1 peak to align
uploc3(k) = find(data(:,4) > (max(data(:,4))/2),1,'first');   %finds first value half of the maximum EB3 peak to align
uploc2 = round((uploc1 + uploc3)/2);
ylength1(k) = ylength-uploc1(k);
ylength3(k) = ylength-uploc3(k);
ylength2(k) = ylength-uploc2(k);

dataset{k}=data;
    
%     figure(1);  %this block makes and saves a figure of the 3 linescans
%     clf;
%     hold on;
%     plot(data(:,1),data(:,2),'g');
%     plot(data(:,1),data(:,3),'b');
%     plot(data(:,1),data(:,4),'r');
%     hold off;
%     legend('EB1','EB2', 'EB3')
%     title('linescan')
%     xlabel('length (\mum)')
%     ylabel('intensity (a.u.)')
%     str = int2str(k);
%     
%     newname = [PathName 'line' str];
%     saveas(gcf, [newname '_plot.eps'],'eps2c');
    
end

uploc = uploc1;
ylength = ylength1;
shortest = max(uploc);
longest = max(ylength);
data1 = zeros((longest+shortest),lines);
data2 = zeros((longest+shortest),lines);
data3 = zeros((longest+shortest),lines);

for k = 1:lines;    %make separate arrays for EB1, EB2 and EB3 data
    
    data = dataset{k};
    data1((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,2);
    data1(1:(shortest-uploc(k)),k) = NaN;
    data1((ylength(k)+shortest+1):end,k) = NaN;
    data2((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,3);
    data2(1:(shortest-uploc(k)),k) = NaN;
    data2((ylength(k)+shortest+1):end,k) = NaN;
    data3((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,4);
    data3(1:(shortest-uploc(k)),k) = NaN;
    data3((ylength(k)+shortest+1):end,k) = NaN;
end

data = zeros((longest+shortest),7);
data(:,1)=[-(scale*(shortest-1)):scale:(scale*longest)];
data(:,2) = nanmean(data1,2);
data(:,3) = nanstd(data1,0,2)./sqrt(sum(~isnan(data1),2));
data(:,4) = nanmean(data2,2);
data(:,5) = nanstd(data2,0,2)./sqrt(sum(~isnan(data2),2));
data(:,6) = nanmean(data3,2);
data(:,7) = nanstd(data3,0,2)./sqrt(sum(~isnan(data3),2));

    figure(1);  %this block makes and saves a figure of the means of linescans
    clf;
    hold on;
    plot(data(:,1),data(:,2),'g');
    plot(data(:,1),data(:,4),'b');
    plot(data(:,1),data(:,6),'r');
    hold off;
    legend('EB1','EB2', 'EB3')
    title('averaged data normalised to EB1')
    xlabel('length (\mum)')
    ylabel('intensity (a.u.)')
   
    
    saveas(gcf, [PathName 'EB1norm_plot.eps'],'eps2c');

savename = [PathName 'EB1norm_data.txt'];
dlmwrite (savename, data, 'delimiter','\t');     
  


uploc = uploc3;     %repeat all for EB3 normalisation
ylength = ylength3;
shortest = max(uploc);
longest = max(ylength);
data1 = zeros((longest+shortest),lines);
data2 = zeros((longest+shortest),lines);
data3 = zeros((longest+shortest),lines);

for k = 1:lines;    %make separate arrays for EB1, EB2 and EB3 data
    
    data = dataset{k};
    data1((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,2);
    data1(1:(shortest-uploc(k)),k) = NaN;
    data1((ylength(k)+shortest+1):end,k) = NaN;
    data2((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,3);
    data2(1:(shortest-uploc(k)),k) = NaN;
    data2((ylength(k)+shortest+1):end,k) = NaN;
    data3((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,4);
    data3(1:(shortest-uploc(k)),k) = NaN;
    data3((ylength(k)+shortest+1):end,k) = NaN;
end

data = zeros((longest+shortest),7);
data(:,1)=[-(scale*(shortest-1)):scale:(scale*longest)];
data(:,2) = nanmean(data1,2);
data(:,3) = nanstd(data1,0,2)./sqrt(sum(~isnan(data1),2));
data(:,4) = nanmean(data2,2);
data(:,5) = nanstd(data2,0,2)./sqrt(sum(~isnan(data2),2));
data(:,6) = nanmean(data3,2);
data(:,7) = nanstd(data3,0,2)./sqrt(sum(~isnan(data3),2));

    figure(2);  %this block makes and saves a figure of the means of linescans
    clf;
    hold on;
    plot(data(:,1),data(:,2),'g');
    plot(data(:,1),data(:,4),'b');
    plot(data(:,1),data(:,6),'r');
    hold off;
    legend('EB1','EB2', 'EB3')
    title('averaged data normalised to EB3')
    xlabel('length (\mum)')
    ylabel('intensity (a.u.)')
    
    saveas(gcf, [PathName 'EB3norm_plot.eps'],'eps2c');

savename = [PathName 'EB3norm_data.txt'];
dlmwrite (savename, data, 'delimiter','\t');  


uploc = uploc2;     %repeat all for EB13 normalisation
ylength = ylength2;
shortest = max(uploc);
longest = max(ylength);
data1 = zeros((longest+shortest),lines);
data2 = zeros((longest+shortest),lines);
data3 = zeros((longest+shortest),lines);

for k = 1:lines;    %make separate arrays for EB1, EB2 and EB3 data
    
    data = dataset{k};
    data1((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,2);
    data1(1:(shortest-uploc(k)),k) = NaN;
    data1((ylength(k)+shortest+1):end,k) = NaN;
    data2((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,3);
    data2(1:(shortest-uploc(k)),k) = NaN;
    data2((ylength(k)+shortest+1):end,k) = NaN;
    data3((shortest-uploc(k)+1):(ylength(k)+shortest),k) = data(:,4);
    data3(1:(shortest-uploc(k)),k) = NaN;
    data3((ylength(k)+shortest+1):end,k) = NaN;
end

data = zeros((longest+shortest),7);
data(:,1)=[-(scale*(shortest-1)):scale:(scale*longest)];
data(:,2) = nanmean(data1,2);
data(:,3) = nanstd(data1,0,2)./sqrt(sum(~isnan(data1),2));
data(:,4) = nanmean(data2,2);
data(:,5) = nanstd(data2,0,2)./sqrt(sum(~isnan(data2),2));
data(:,6) = nanmean(data3,2);
data(:,7) = nanstd(data3,0,2)./sqrt(sum(~isnan(data3),2));

    figure(3);  %this block makes and saves a figure of the means of linescans
    clf;
    hold on;
    errorbar(data(:,1),data(:,2),data(:,3),'g');
    errorbar(data(:,1),data(:,4),data(:,5),'b');
    errorbar(data(:,1),data(:,6),data(:,7),'r');
    hold off;
    legend('EB1','EB2', 'EB3')
    title('averaged data normalised to middle of EB1 and EB3')
    xlabel('length (\mum)')
    ylabel('intensity (a.u.)')
    xlim([-2 5]);
    
    saveas(gcf, [PathName 'EB13norm_plot.pdf']);

savename = [PathName 'EB13norm_data.txt'];
dlmwrite (savename, data, 'delimiter','\t');  



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
