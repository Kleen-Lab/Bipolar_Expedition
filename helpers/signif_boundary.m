function [clusterSig,pValsSig,boundarySig,clustX,clustY] = signif_boundary(cluster,pVals,x,y)
% x,y - size of matrices
% 

numberClusters = length(cluster);
clusterSig = {};
pValsSig = {};
boundarySig = {};
clustX  = {};
clustY = {};

for jj = 1:numberClusters
   if pVals(jj) <= 0.05
      clusterSig{jj} = cluster{jj};
      pValsSig{jj} = pVals(jj);
      [a,b] = ind2sub([x,y],cluster{jj});
      k = boundary(a,b,1);
      clustX{jj} = a;
      clustY{jj} = b;
      boundarySig{jj} = k;
   end
    
end



end

