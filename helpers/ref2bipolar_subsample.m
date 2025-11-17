function [d] = ref2bipolar_subsample(d,subj)
% d is a matrix of samples by channels by trials consisting of referential intracranial EEG data
% this sub samples the rereferenced data by skipping every intervening
% electrode

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
            oddRow = true;
            for i=bipolarN(r,1):N+1:bipolarN(r,2)
                if oddRow
                    d(:,i:2:i+N-2,jj)=d(:,i:2:i+N-2,jj)-d(:,i+2:2:i+N,jj);
                    d(:,2:2:i+N,jj) = NaN;% even chans will be NaNs, including last channel
                    oddRow = false;
                else
                    d(:,i:i+N,jj) = NaN;
                    oddRow = true;
                end
            end
        else; i=bipolarN(r,1); N=bipolarN(r,2)-i; %strips, depths (1-D)
            d(:,i:2:i+N-2,jj)=d(:,i:2:i+N-2,jj)-d(:,i+2:2:i+N,jj);
            d(:,2:2:i+N,jj) = NaN;% even chans will be NaNs, including last channel

        end; clear N i
    end
end

end
