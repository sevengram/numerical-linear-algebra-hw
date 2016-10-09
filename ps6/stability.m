normA = 0;
normL0 = 0;
normU0 = 0;
normdelta0 = 0;
normL1 = 0;
normU1 = 0;
normdelta1 = 0;
for i=1:50
    A = rand(160);
    A(1,1) = 1e-13;
    normA = normA + norm(A,inf);

    [~,L,U] = myPLU(A,0);
    normL0 = normL0 + norm(L,inf);
    normU0 = normU0 + norm(U,inf);
    normdelta0 = normdelta0 + norm(L*U-A,inf);

    [p,L,U] = myPLU(A,1);
    normL1 = normL1 + norm(L,inf);
    normU1 = normU1 + norm(U,inf);
    normdelta1 = normdelta1 + norm(L*U-A(p,:),inf);
end

normA = normA/50

normL0 = normL0/50
normU0 = normU0/50
normdelta0 = normdelta0/50
normL1 = normL1/50
normU1 = normU1/50
normdelta1 = normdelta1/50

a = normdelta0 / normA
b = normdelta1 / normA
