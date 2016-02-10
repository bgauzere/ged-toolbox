function [linearCost] = linearCost(CostMatrix, PermMatrix)
    linearCost = sum(sum(CostMatrix(PermMatrix==1)));
end
