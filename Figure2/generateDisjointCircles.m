function [xC,r] = generateDisjointCircles(a,b,c,d,minRad,maxRad,DIST,howMany)

% this function generates "howMany" disjoint circles, each of which has
% radius in the interval [minRad,maxRad], and center in the rectangle
% defined by [a,b] x [c,d]. No two circles will be within "DIST" of each
% other. 

xC = [];
r = [];

k = 1;
while k <= 10000 && length(r)<howMany
   x = a + (b-a)*rand(1,1); y = c + (d-c)*rand(1,1);
   RAD = minRad + (maxRad - minRad)*rand(1,1);
   if isempty(r)
      r = RAD;
      xC = [x;y];
   else
       check = 1;
       for m = 1:length(r)
            check = check*(norm(xC(:,m) - [x;y],2)>=r(m)+RAD+DIST); 
       end
       if check
           xC(:,end+1) = [x;y]; r(end+1) = RAD;
       end
   end
   k=k+1;
end

if length(r) < howMany
   fprintf('Failed to generate %i disjoint circles with centers in [%.2f,%.2f] x [%.2f,%.2f] in %i trials.\n',howMany,a,b,c,d,k); 
else
   fprintf('Generated %i disjoint circles with centers in [%.2f,%.2f] x [%.2f,%.2f] in %i trials.\n',howMany,a,b,c,d,k);  
end
end

