function dst=distance3D(e1,e2)
% Finds 3-D euclidean distance for data. 
% a and b should each be [x y z]
s=size(e1); if ndims(e1)>2; disp('inputs incorrect, need [x y z]'); return; end %#ok<ISMAT>
if size(e1)~=size(e2); disp('inputs not same dimensions'); return; end
if (s(1)==1 | s(1)==3) & (s(2)==1 | s(2)==3) & s(1)~=s(2)
if size(e1,1)==3; e1=e1'; e2=e2'; end
dst=sqrt((e1(:,1)-e2(:,1))^2 + (e1(:,2)-e2(:,2))^2 + (e1(:,3)-e2(:,3))^2);
% sqrt((x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2)
else disp('inputs incorrect, need [x y z]'); return; 
end

% end
%     
% % if ndims(a)==2 && any(size(a))==1;
% %     d=
% if any(size(a))== && any(size(a))==3;
%     
%     
% elseif ndims(a)==2 && any(size(a))==1;
%     
% end


