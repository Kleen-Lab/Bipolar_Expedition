function [elecs,elecmatrix,eleclabels,anatomy]=getelecs_region(pt,reg,clin1TDT2)
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
% 'ip' --- 'inferiorparietal'
% 'sp' --- 'superiorparietal'
% 'precun' --- 'precuneus'
% 'cun' --- 'cuneus'
% 'lin' --- 'lingual'

elecs=[];
[em,~,anatomy]=getelecs(pt,clin1TDT2); 
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
if mean(em(:,1))>0; hem='Right'; else hem='Left'; end; hem=[hem '-']; 
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
z=strcmp(reg{i},'rmf'); if any(z); reg(i)={'rostralmiddlefrontal'}; end
z=strcmp(reg{i},'sf'); if any(z); reg(i)={'superiorfrontal'}; end
z=strcmp(reg{i},'sm'); if any(z); reg(i)={'supramarginal'}; end
z=strcmp(reg{i},'mof'); if any(z); reg(i)={'medialorbitofrontal'}; end
z=strcmp(reg{i},'lof'); if any(z); reg(i)={'lateralorbitofrontal'}; end
z=strcmp(reg{i},'fp'); if any(z); reg(i)={'frontalpole'}; end
z=strcmp(reg{i},'prec'); if any(z); reg(i)={'precentral'}; end
z=strcmp(reg{i},'postc'); if any(z); reg(i)={'postcentral'}; end
z=strcmp(reg{i},'ip'); if any(z); reg(i)={'inferiorparietal'}; end
z=strcmp(reg{i},'sp'); if any(z); reg(i)={'superiorparietal'}; end
z=strcmp(reg{i},'precun'); if any(z); reg(i)={'precuneus'}; end
z=strcmp(reg{i},'cun'); if any(z); reg(i)={'cuneus'}; end
z=strcmp(reg{i},'lin'); if any(z); reg(i)={'lingual'}; end
end

if isc; reg=reg{1}; end
if ischar(reg); [elecmatrix,eleclabels,anatomy]=getelecs(pt,clin1TDT2); elecs=find(strcmp(anatomy(:,4),reg))';
                                elecmatrix=elecmatrix(elecs,:); eleclabels=eleclabels(elecs,:); anatomy=anatomy(elecs,:); 
else
for i=1:length(reg); 
    [~,~,anatomy]=getelecs(pt,clin1TDT2); elecs=[elecs find(strcmp(anatomy(:,4),reg{i}))'];
end
end

