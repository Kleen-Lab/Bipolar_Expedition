function [d] = ref2bipolar_reg(d,subj)
% d is a matrix of samples by channels by trials consisting of referential intracranial EEG data

[bipolarN,bipolarT]=xlsread(fullfile('/Users/davidcaldwell/Box/KLEENLAB/KLEEN_DRIVE copy/David/Bipolar project','AN_ElectrodeInfoTDT.xlsx'),subj);

% account for case where subject has less channels recorded, chop off last
% part of bipolarN row 

numElecs = size(d,2);
if numElecs < bipolarN(end,2)
   bipolarN(end,2) = numElecs; 
end

for jj = 1:size(d,3)
    for r=1:size(bipolarT,1) %each row of the sheet is a component (grid, strip, or depth)
        if any(strcmpi(bipolarT(r,2),{'grid','minigrid'})) %grids (2-D)
            N=bipolarN(r,3)-1;
            for i=bipolarN(r,1):N+1:bipolarN(r,2)
                d(:,i:i+N-1,jj)=d(:,i:i+N-1,jj)-d(:,i+1:i+N,jj);
                d(:,i+N)=NaN; %last channel in the line will be NaNs
            end
        else; i=bipolarN(r,1); N=bipolarN(r,2)-i; %strips, depths (1-D)
            d(:,i:i+N-1,jj)=d(:,i:i+N-1,jj)-d(:,i+1:i+N,jj);
            d(:,i+N,jj)=NaN; %last channel in the line will be NaNs
        end; clear N i
    end
end

end
