function [c_h] = ctmr_gauss_plot(cortex,elecmatrix,weights,hemi,do_lighting,cax,addl,CM,gsp)
% function [c_h]=ctmr_gauss_plot(cortex,elecmatrix,weights)
%
% projects electrode locationsm (elecmatrix) onto their cortical spots in 
% the left hemisphere and plots about them using a gaussian kernel
% for only cortex use:
% ctmr_gauss_plot(cortex,[0 0 0],0)
% rel_dir=which('loc_plot');
% rel_dir((length(rel_dir)-10):length(rel_dir))=[];
% addpath(rel_dir)

%     Copyright (C) 2009  K.J. Miller & D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

%   Version 1.1.0, released 26-11-2009
%
% Edited 2016 by Liberty Hamilton
% Edited 2017-21 by Jon Kleen


if isempty(elecmatrix); %nargin<2
    elecmatrix = [0 0 0];
end
if isempty(weights); %nargin<3 || isempty(weights)
    weights = zeros(size(elecmatrix,1),1);
end
if ~exist('hemi','var')
    hemi = 'diag';
end
if nargin<5
    do_lighting=1;
end

if ~exist('addl','var'); addl=0; end

%load in colormap
% load('loc_colormap')
% load('loc_colormap_thresh')
% load('BlWhRdYl_colormap')
% load('BlGyOrCp_colormap')
% cm = flipud(cbrewer('div','RdGy',64));

% cm = flipud(cbrewer('div','RdGy',64)); cm=cm(:,[1 3 1]); %purple brain, purple brain
% cm = flipud(cbrewer('div','RdGy',64)); cm=[cm(1:3:32,[1 3 1]); cm(27:32,[1 3 1]); [1 1 1]; flipud(cm(33:38,[3 1 1])); cm(33:3:64,[1 3 1])];
% cm = flipud(cbrewer('div','RdGy',64)); cm=[cm(1:32,[1 3 1]); ones(32,3); cm(33:64,[1 3 1])]; 
% cm = flipud(cbrewer('div','RdGy',64)); cm=cm(:,[3 1 3]); 

% cm=-cool+1; a=[0:1/32:1]'; a(end)=[]; cm(1:32,:)=cm(1:32,:).*[a a a]; cm=-[zeros(size(cm,1)+1,3); cm]+1;
if exist('CM','var') && ~isempty(CM); cm=CM; else cm=cmSz3D; end
% cm(65,:)=brainskin;

brain=cortex.vert;

if length(weights)~=length(elecmatrix(:,1))
    error('you sent a different number of weights than elecmatrix (perhaps a whole matrix instead of vector)')
end

%gaussian "cortical" spreading parameter - in mm, so if set at 10, its 1 cm
%- distance between adjacent elecmatrix
if ~exist('gsp','var') || isempty(gsp); gsp=10; end %default 10, edited from 50


% % original version
c=zeros(length(cortex(:,1)),1);
for i=1:length(elecmatrix(:,1))
    b_z=abs(brain(:,3)-elecmatrix(i,3));
    b_y=abs(brain(:,2)-elecmatrix(i,2));
    b_x=abs(brain(:,1)-elecmatrix(i,1));
    %d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2).^.5)/gsp^.5); %exponential fall off
    d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian
    c=c+d';
end

% new version but slow
% c=zeros(size(brain(:,1),1),1);
% for i=1:length(brain(:,1)); 
%     b_z=abs(brain(i,3)-elecmatrix(:,3));
%     b_y=abs(brain(i,2)-elecmatrix(:,2));
%     b_x=abs(brain(i,1)-elecmatrix(:,1));
%     w=weights.*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian
%    c(1,i)=max(w);
% end


% % new version but ?fast
% c=zeros(size(brain(:,1),1),size(elecmatrix,1));
% for i=1:length(elecmatrix(:,1))
%     b_z=abs(brain(:,3)-elecmatrix(i,3));
%     b_y=abs(brain(:,2)-elecmatrix(i,2));
%     b_x=abs(brain(:,1)-elecmatrix(i,1));
%     d=weights(i)*exp((-(b_x.^2+b_z.^2+b_y.^2))/gsp); %gaussian
%     c(:,i)=d';
% end
% c = max(c,[],2)';




% c=(c/max(c));
c_h=tripatch(cortex, 'nofigure', c');
if ~addl; shading interp; end
a=get(gca);

%%NOTE: MAY WANT TO MAKE AXIS THE SAME MAGNITUDE ACROSS ALL COMPONENTS TO REFLECT
%%RELEVANCE OF CHANNEL FOR COMPARISON's ACROSS CORTICES
d=a.CLim;
%set(gca,'CLim',[-max(abs(d)) max(abs(d))])
if exist('cax','var') && ~isempty(cax); set(gca,'CLim',[cax(1) cax(2)]); 
else set(gca,'CLim',[-max(abs(d)) max(abs(d))]); 
end
colormap(gca,cm)
lighting phong; %play with lighting...
%material shiny;
material dull;
% material([.3 .8 .1 10 1]);
% material([.2 .9 .2 50 1]); %  BF: editing mesh viewing attributes
axis off

if do_lighting; 
    %delete(findall(gcf,'Type','light')); %NOTE: this line is usually culprit for subplot/loop issues
    l=light;
    
    if strcmp(hemi,'lh')
        view(270, 0);
        set(l,'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi,'rh')
        view(90, 0);
        set(l,'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi, 'top')
        view(0, 90);
        set(l,'Position',[0 0 1],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi, 'bottom')
        view(0, -90);
        set(l,'Position',[0 0 -1],'Color',[0.8 0.8 0.8]);
    elseif strcmp(hemi, 'diag')
        view(135, 0);
        set(l,'Position',[1 0 1],'Color',[0.8 0.8 0.8]);
    end
end
