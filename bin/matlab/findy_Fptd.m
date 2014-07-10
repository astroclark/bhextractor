function y=findy_Fptd(V, b1, numPCs,d)

y=d*V(:,1:numPCs)*(b1(:,1:numPCs))';


