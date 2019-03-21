Decompose2 := function(m)
     local d,k,p,slp,m1,m2,m3;

     p := Size(DefaultFieldOfMatrix(m));
     d := Z(p);

     #Algorithim will look for triangular matrices first:
     if IntFFE(m[1][1]) = 1 and IntFFE(m[2][2]) = 1 then
          #Looking for the cases where [[1,d],[0,1]]
          if IsUpperTriangularMat(m) then
               slp:= StraightLineProgram([[2,IntFFE(m[1][2])]],3);
               return slp;
               #k := LogFFE(m[1][2],d);
               #if k mod 2 = 0 then
               #     k := k/2;
               #     slp := StraightLineProgram([[3,k,2,1,3,-k]],3);
               #     return slp;
               #else
               #     k := (k-1)/2;
               #     slp := StraightLineProgram([[3,k,2,IntFFE(d),3,-k]],3);
               #     return slp;
               #fi;
          #Looking for the cases where [[1,0],[d,1]]
          elif IsLowerTriangularMat(m) then
               #k := LogFFE(m[2][1],d);
               slp := StraightLineProgram([[1,-1,2,IntFFE(-m[2][1]),1,1]],3);
               return slp;
          fi;
     #Looking for the cases where [[a,b],[c,d]] with b != 0
     elif IntFFE(m[1][2]) <> 0 then
          #z:=m[2][1]+m[2][2]*(d^0-m[1][1])/m[1][2];
          m3:=[[1,0],[(m[2][1]+m[2][2]*(d^0-m[1][1])/m[1][2]),1]]*d^0;
          m2:=[[1,m[1][2]],[0,1]]*d^0;
          m1:=[[1,0],[(m[1][1]-d^0)/m[1][2],1]]*d^0;
          slp:=ProductOfStraightLinePrograms(Decompose2(m2),Decompose2(m1));
          slp:=ProductOfStraightLinePrograms(Decompose2(m3),slp);
          return slp;
     #Looking for the cases where [[a,b],[c,d]] with b == 0
     else
          m1:=[[0,1],[-1,0]]*d^0;
          m:=m*m1;
          slp:=ProductOfStraightLinePrograms(Decompose2(m),Decompose2(Inverse(m1)));
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
