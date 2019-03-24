# This function will return a SLP decomposition of 2x2 upper/lower triangular matrices in SL(2,p) with the generators of ClassicalStandardGenerators S T D
Triangular2 := function(m)
     local slp;

     #Looking for the cases where [[1,x],[0,1]]
     if IsUpperTriangularMat(m) then
          slp:= StraightLineProgram([[2,IntFFE(m[1][2])]],3);
          return slp;
     #Looking for the cases where [[1,0],[x,1]]
     elif IsLowerTriangularMat(m) then
          slp := StraightLineProgram([[1,-1,2,IntFFE(-m[2][1]),1,1]],3);
          return slp;
     fi;
end;

# This function will return a SLP decomposition of any 2x2 matrix in SL(2,p)
Decompose2 := function(m)
     local d,k,p,slp,m1,m2,m3;

     p := Size(DefaultFieldOfMatrix(m));
     d := Z(p);

     #Algorithim will look for triangular matrices first
     if IntFFE(m[1][1]) = 1 and IntFFE(m[2][2]) = 1 then
          return Triangular2(m);
     #Looking for the cases where [[a,b],[c,d]] with b != 0
     elif IntFFE(m[1][2]) <> 0 then
          m3:=[[1,0],[(m[2][1]+m[2][2]*(d^0-m[1][1])/m[1][2]),1]]*d^0;
          m2:=[[1,m[1][2]],[0,1]]*d^0;
          m1:=[[1,0],[(m[1][1]-d^0)/m[1][2],1]]*d^0;

          slp:=ProductOfStraightLinePrograms(Triangular2(m2),Triangular2(m1));
          slp:=ProductOfStraightLinePrograms(Triangular2(m3),slp);
          return slp;
     #Looking for the cases where [[a,b],[c,d]] with b == 0
     else
          m1:=[[0,1],[-1,0]]*d^0;
          m:=m*m1;
          slp:=ProductOfStraightLinePrograms(Decompose2(m),Triangular2(Inverse(m1)));
          return slp;
     fi;
end;

ClassicalStandardGenerators := function(m,p)
     local d,S,T,D;
     d := Z(p);
     S := [[0,1],[-1,0]]*d^0;
     T := [[1,1],[0,1]]*d^0;
     D := [[d,0],[0,d^-1]]*d^0;
     return [S,T,D];
end;
