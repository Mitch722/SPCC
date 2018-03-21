function noViol = countViolations(y, bnds)

noViol = sum((abs(y) > [bnds(1,1); bnds(2,1)]), 2);
noViol = sum(noViol, 1);

end