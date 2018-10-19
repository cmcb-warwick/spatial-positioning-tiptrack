function [d] = DCK_select(varargin)

% PURPOSE:
%   This allows manual selection of constant growth regions from kymographs
%   for further analysis of relative EB peak positions in dual-colour
%   experiments.

% INPUT:  
%   none required. Calling the function will request first selecting a folder
%   with kymographs and second a folder in which to save the results. The
%   first image will be opened in ImageJ and user is required to use the
%   straight line tool in ImageJ to select growth phases and add to ROI
%   manager. Once all segments are selected for an image, enter anything
%   into the MATLAB prompt to save selected regions and open the next
%   image. If done, write Quit.

% OUTPUT:
%   saves cropped kymographs and MT.mat file for further analysis using the DCK_Kymo_Run function into the selected results folder.


% Copyright: Straube lab, University of Warwick, July 2016
%   shared under New BSD license as specified below
%
% Code was tested under Matlab release 2012b


%--------------------------------------------------
% Load and Save Directory - Asking
%--------------------------------------------------

if isempty(varargin)
    load_loc{1} = uigetdir(pwd,'Pick Load Directory');
else
    load_loc{1} = varargin{1};
end

if length(varargin) >= 2
    save_loc = varargin{2};
else
    save_loc = uigetdir(load_loc{1},'Pick Save Directory');
end

if ~exist(save_loc,'dir')
    mkdir(save_loc)
end

%--------------------------------------------------
% Finding Tiff Files
%--------------------------------------------------

r_files = loadData(load_loc,3,'/*red*.tif');
g_files = loadData(load_loc,3,'/*green*.tif');

r_files = sortrows(r_files)
g_files = sortrows(g_files);

L = length(r_files);

format compact

javaaddpath([matlabroot '/java/ij.jar'])
javaaddpath([matlabroot '/java/mij.jar'])
MIJ.start;
opener = ij.io.Opener();

mt = 0;

for i = 1:L
    
    fprintf(['Starting ' num2str(i,'%03d') ' out of ' num2str(L,'%03d') '\n']);
    
    % obtains file structure
    
    clear a b
    [a,b{1},~] = fileparts(r_files{i});
    [a,b{2},~] = fileparts(a);
    [~,b{3},~] = fileparts(a);
    clear a
    
    % loads and shows file
    
    fprintf('loading ... \n')
    name   = java.lang.String(r_files{i});
    opener.openImage(name).show();
    imp    = ij.IJ.getImage();
    xlim   = imp.getWidth();
    setZoom(imp,300)

    % starts ROI Manager
    
    ij.plugin.frame.RoiManager;
    ij.IJ.setTool('line');  
    rm = ij.plugin.frame.RoiManager.getInstance();
    rm.runCommand('Show All');
    
    % command line inputs
    
    fprintf('"Quit" will end process, any other input ask to load a new file \n')
    ins = input('Next Action: ','s');
    
    % obtains ROIs
    
    n    = rm.getCount();
    rois = rm.getROIs();
    list = rm.getList();
    
    % loops through ROIs saving cropped images
    clear Temp_MT
    
    % Red Channel
    for j = 0:n-1
        
        name        = list.getItem(j);
        roi         = rois.get(name);
        px          = [roi.x1; roi.x2];
        py          = [roi.y1; roi.y2];
        
        a1 = find(b{1} == '_')+1;
        a3 = b{1};
        a3 = a3(a1:end);
        
        s_loc   = [save_loc '/' b{3} '/' 'Cropped' '/'];
        if ~exist(s_loc,'dir'); mkdir(s_loc); end
        s_loc   = [s_loc b{2} '_' a3 '_' num2str(j,'%02d') '_Red.tif'];
        fprintf(['saving ... ' s_loc '\n'])
        s_loc   = java.lang.String(s_loc);
        
        [x,flip] = ROIcrop(px,py,xlim);
            
        imp.setRoi(x(1),x(2),x(3)-x(1),x(4)-x(2))
        ip = imp.duplicate();
        
        if flip == 1;
            ip.getProcessor().flipHorizontal();
        end
        
        saver = ij.io.FileSaver(ip);
        saver.saveAsTiff(s_loc);
        
        Temp_MT(j+1).Red_Parent  = r_files{i};
        Temp_MT(j+1).Red_File    = char(s_loc);
        Temp_MT(j+1).Crop_Box    = x;
        Temp_MT(j+1).Flipped     = flip;
        Temp_MT(j+1).EB_Proteins = b{3};
        Temp_MT(j+1).Stack       = b{2};
        
    end
    
    rm.close;
    ij.IJ.run('Close All');
    
    name   = java.lang.String(g_files{i});
    opener.openImage(name).show();
    imp    = ij.IJ.getImage();
    
    % Repeat for Green Channel
    for j = 0:n-1
        
        s_loc   = Temp_MT(j+1).Red_File;
        a       = find(s_loc == '_',1,'last');
        
        s_loc   = [s_loc(1:a) 'Green.tif'];
        fprintf(['saving ... ' s_loc '\n'])
        s_loc   = java.lang.String(s_loc);
        
        x       = Temp_MT(j+1).Crop_Box;
        flip    = Temp_MT(j+1).Flipped;
            
        imp.setRoi(x(1),x(2),x(3)-x(1),x(4)-x(2))
        ip = imp.duplicate();
        
        if flip == 1;
            ip.getProcessor().flipHorizontal();
        end
        
        saver = ij.io.FileSaver(ip);
        saver.saveAsTiff(s_loc);
        
        Temp_MT(j+1).Green_Parent = g_files{i};
        Temp_MT(j+1).Green_File   = char(s_loc);
        
        mt = mt+1;
        MT(mt) = Temp_MT(j+1);
        
    end
    
    rm.close;
    ij.IJ.run('Close All');
    
     
    switch ins
        case 'Quit'
            break
    end
    
end

save([save_loc '/MT.mat'],'MT');
MIJ.exit

% ** Based on Albert Cardona's ZoomExact plugin:
% http://albert.rierol.net/software.html */
function [] = setZoom(imp,mag) 

    win = imp.getWindow();
    mag = mag/100;
    ic = imp.getCanvas();

    win.getCanvas().setMagnification(mag);
    w = imp.getWidth()*mag;
    h = imp.getHeight()*mag;

    ic.setDrawingSize(w,h);
    win.pack();
    ic.repaint();

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