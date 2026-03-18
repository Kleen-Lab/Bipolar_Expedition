function [elecs,em,anatomy]=getregionelecs(pt,reg)
% index of electrodes that are in the anatomical area (region) specified
% for the patient (pt) specified. 
% Note: em's 3 columns are ML(R+), AP(A+), DV (D+) 
% __common region names below:__
% 
% 'stg' --- 'superiortemporal'
% 'mtg' --- 'middletemporal'
% 'itg' --- 'inferiortemporal'
%
% 'ent' --- 'entorhinal'
% 'fus' --- 'fusiform'
% 'ph' --- 'parahippocampal'
% 'tp' --- 'temporalpole'
%
% 'hp' --- 'Right-Hippocampus' or
%          'Left-Hippocampus'
% 'am' --- 'Right-Amygdala' or
%          'Left-Amygdala'
%
% 'sm' --- 'supramarginal'
% 'prec' --- 'precentral'
% 'postc' --- 'postcentral'
% 'mof' --- 'medialorbitofrontal'
% 'lof' --- 'lateralorbitofrontal'
% 'pt' --- 'parstriangularis'
% 'pop' --- 'parsopercularis'
% 'por' --- 'parsorbitalis'
% 'rmf' --- 'rostralmiddlefrontal'
% 'sf' --- 'superiorfrontal'
% 'fp' --- 'frontalpole'

elecs=[];
[em,~,anatomy]=getelecs(pt,2); %only using TDT here
% if any(strcmp(pt,{'EC129','EC137','EC163','EC166','EC173'})); hem='Right-'; else hem='Left-'; end

if ~exist('reg','var') || (exist('reg','var') && isempty(reg)); 
    [u,~,~]=unique(anatomy(:,4));
    for i=1:length(u); u(i,2)={length(find(strcmp(anatomy(:,4),u{i})))}; end
    [~,s]=sort(cell2mat(u(:,2))); u=u(flipud(s),:);
%     disp(strcat({'Available regions for '},pt,':'))
%     disp({'REGION','[NUMBER OF ELECTRODES]'})
%     disp(u)
    return
end
if nanmean(em(:,1))>0; hem='Right'; else hem='Left'; end; hem=[hem '-']; 
isc=0;
if ischar(reg); reg={reg}; isc=1; end
for i=1:length(reg);
z=strcmp(reg{i},'stg'); if any(z); reg(i)={'superiortemporal'}; end
z=strcmp(reg{i},'mtg'); if any(z); reg(i)={'middletemporal'}; end
z=strcmp(reg{i},'itg'); if any(z); reg(i)={'inferiortemporal'}; end
z=strcmp(reg{i},'fus'); if any(z); reg(i)={'fusiform'}; end
z=strcmp(reg{i},'tp'); if any(z); reg(i)={'temporalpole'}; end
z=strcmp(reg{i},'ph'); if any(z); reg(i)={'parahippocampal'}; end
z=strcmp(reg{i},'hp'); if any(z); reg(i)={[hem 'Hippocampus']}; end
z=strcmp(reg{i},'am'); if any(z); reg(i)={[hem 'Amygdala']}; end
z=strcmp(reg{i},'ent'); if any(z); reg(i)={'entorhinal'}; end
z=strcmp(reg{i},'pt'); if any(z); reg(i)={'parstriangularis'}; end
z=strcmp(reg{i},{'po','pop'}); if any(z); reg(i)={'parsopercularis'}; end
z=strcmp(reg{i},{'por','porb'}); if any(z); reg(i)={'parsorbitalis'}; end
z=strcmp(reg{i},'cmf'); if any(z); reg(i)={'caudalmiddlefrontal'}; end
z=strcmp(reg{i},'rmf'); if any(z); reg(i)={'rostralmiddlefrontal'}; end
z=strcmp(reg{i},'sf'); if any(z); reg(i)={'superiorfrontal'}; end
z=strcmp(reg{i},'sm'); if any(z); reg(i)={'supramarginal'}; end
z=strcmp(reg{i},'mof'); if any(z); reg(i)={'medialorbitofrontal'}; end
z=strcmp(reg{i},'lof'); if any(z); reg(i)={'lateralorbitofrontal'}; end
z=strcmp(reg{i},'fp'); if any(z); reg(i)={'frontalpole'}; end
z=strcmp(reg{i},'prec'); if any(z); reg(i)={'precentral'}; end
z=strcmp(reg{i},'postc'); if any(z); reg(i)={'postcentral'}; end
end
z=strcmp(reg,'mfg'); if any(z); reg={'caudalmiddlefrontal','rostralmiddlefrontal'}; isc=0; end
z=strcmp(reg,'cing'); if any(z); reg={'caudalanteriorcingulate','isthmuscingulate','posteriorcingulate','rostralanteriorcingulate'}; isc=0; end


if isc; reg=reg{1}; end
if ischar(reg); [~,~,anatomy]=getelecs(pt,2); elecs=find(strcmp(anatomy(:,4),reg))';
else
for i=1:length(reg); 
    [~,~,anatomy]=getelecs(pt,2); elecs=[elecs find(strcmp(anatomy(:,4),reg{i}))'];
end
end

