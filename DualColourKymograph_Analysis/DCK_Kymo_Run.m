function [d] = DCK_Kymo_Run(varargin)

% PURPOSE:
%   This corrects for temporal offset if using sequential imaging and collects 
%   all data from intensity profiles along the microtubule tip for further analysis 
%   of relative EB peak positions in dual-colour experiments.

% INPUT:  
%   none required. Calling the function will request first selecting the folder
%   containing the MT.mat file generated when cropping growth regions using
%   DCK_Kymo_select.m 
%   In the same folder shouls also be a file called offset.mat, which
%   contains information on all stacks in which the red channel was imaged 
%   before the green channel. This should be a structure with offset.P containing 
%   the name of the protein used and offset.F containing the stack numbers. 

% OUTPUT:
%   saves DCA_Kymo_Run.mat file for further analysis using the DCK_Curve_Analysis function.


% Copyright: Straube lab, University of Warwick, July 2016
%   shared under New BSD license as specified below
%
% Code was tested under Matlab release 2012b



%--------------------------------------------------
% Initial Parameters
%--------------------------------------------------

sres0 = 81;
tres0 = 1;
sbin  = 10;
tbin  = 2;
sres1 = sres0/sbin;
tres1 = tres0/tbin;

%--------------------------------------------------
% Load and Save Directory - Asking
%--------------------------------------------------

if isempty(varargin)
    load_loc = uigetdir(pwd,'Pick Load Directory');
else
    load_loc = varargin{1};
end

% if exist(load_loc,'dir')
%     if length(varargin) >= 2
%         save_loc = varargin{2};
%     else
%         save_loc = uigetdir(load_loc,'Pick Save Directory');
%     end
% 
%     if ~exist(save_loc,'dir')
%         mkdir(save_loc)
%     end
%     
%     load_loc = DCK_cropped(load_loc,save_loc);
% end


MT = load([load_loc '/MT.mat']);
MT = MT.MT;
[save_loc,~,~] = fileparts(load_loc);

%--------------------------------------------------
% Offsets in Imaging
%--------------------------------------------------

clear offset
% [parent_loc,~,~] = fileparts(save_loc);
offset = load([load_loc '/offset.mat']);
offset = offset.offset;



L = length(MT);
for i = 1:L
    
    fprintf(['Starting ' num2str(i,'%03d') ' out of ' num2str(L,'%03d') '\n']);
    
    O_r0 = imread(MT(i).Red_File);
    O_g0 = imread(MT(i).Green_File);
    
    if ~(size(O_r0) == size(O_g0))
        i
        continue
    end
        

    
% Comment this lot in if we have different stacks imaged in different order to correct for the time offset in imaging these    
    
%     sit = 0;
%     state = 0;
%     while sit == 0 && state < length(offset)
%         state = state + 1;
%         if strcmp(offset(state).P,MT(i).EB_Proteins)
%             sit = 1;
%         end
%     end
% 
%     if state == 7
%         continue
%     end
% 
%     off = 0;
%     
%     if sum(str2num(MT(i).Stack) == offset(state).F) == 1;
         off = 0;
%     end
    
    [X0,T0] = meshgrid(1:size(O_r0,2),1:size(O_r0,1));
    [X1,T1] = meshgrid(1:size(O_r0,2),tbin:(size(O_r0,1)*tbin));
    
    T1 = T1/tbin;

    O_r1 = interp2(X0,T0,O_r0,X1,T1);
    O_g1 = interp2(X0,T0,O_g0,X1,T1);
    
    [X1,T1] = meshgrid(1:size(O_r1,2),1:size(O_r1,1));
    [X2,T2] = meshgrid(sbin:(size(O_r1,2)*sbin),1:size(O_r1,1));

    X2 = X2/sbin;

    O_r2 = interp2(X1,T1,O_r1,X2,T2);
    O_g2 = interp2(X1,T1,O_g1,X2,T2);

    if off == 1
        O_g2(1,:)   = [];
        O_r2(end,:) = [];
    else
        O_g2(end,:) = [];
        O_r2(1,:)   = [];
    end
    
    MT(i).sres0 = sres0;
    MT(i).tres0 = tres0;
    MT(i).sres1 = sres1;
    MT(i).tres1 = tres1;
    
    MT(i).Reference = 'Green';
    MT(i).Save_Loc  = save_loc;
    
    Temp_MT = DCK_Straight_Kymo_Space_C1(MT(i),O_r2);
    Temp_MT = DCK_Straight_Kymo_Space_C2(Temp_MT,O_g2);
    
    MT(i).Velo              = Temp_MT.Velo;
    MT(i).Half_vector_C1    = Temp_MT.Half_vector_C1;
    MT(i).Half_Intensity_C1 = Temp_MT.Half_Intensity_C1; 
    MT(i).half_length_C1    = Temp_MT.half_length_C1;
    MT(i).au_curve_C1       = Temp_MT.au_curve_C1;
    MT(i).au_normalised_C1  = Temp_MT.au_normalised_C1;
    MT(i).peak_width_C1     = Temp_MT.peak_width_C1;
    MT(i).curve_C1          = Temp_MT.curve_C1;
    MT(i).error_C1          = Temp_MT.error_C1;
    MT(i).half_length_C2    = Temp_MT.half_length_C2;
    MT(i).au_curve_C2       = Temp_MT.au_curve_C2;
    MT(i).au_normalised_C2  = Temp_MT.au_normalised_C2;
    MT(i).peak_width_C2     = Temp_MT.peak_width_C2;
    MT(i).curve_C2          = Temp_MT.curve_C2;
    MT(i).error_C2          = Temp_MT.error_C2;
    
    if isfield(Temp_MT,'shift')
        MT(i).shift = Temp_MT.shift;
    end

    clear Temp_MT
    

end

 
save([save_loc '/DCA_Kymo_Run.mat'],'MT')

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