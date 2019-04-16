Decompose2 := function(m,q)
     local f, d, z, k, m1, m2, m3, basis, coeff, slps, slp, SLP;

     f := GF(q); #DefaultFieldOfMatrix(m); # Finite field of the matrix
     #q := Size(f); # q = p^d which is the size of the finite field
     d := Size(GaloisGroup(f)); # d from q = p^d
     z := Z(q); # Primitive element of the finite field

     # First case: [[1,a],[0,1]]
     if m[1][1] = One(f) and m[2][1] = Zero(f) and m[2][2] = One(f) then
          # For that we're assuming that z^0, z^2 , ... , z^(d-1) is a basis,
          # therefore a = a_0 z^0 + a_2 z^2 + ...

          # Creating the basis
          basis := [];
          k := 0;
          while k < d do
               Append(basis,[z^(2*k)]);
               k := k + 1;
          od;
          basis := Basis(f,basis);

          # Returning the coefficients of m[1][2] in that base
          coeff := Coefficients(basis,m[1][2]);

          # Computing slps to produce the product at the end
          slps := [];
          k := 0;
          while k < d do
               slp := StraightLineProgram([ [[3,k,2,1,3,-k],4], [4,IntFFE(coeff[k+1])] ],3);
               # Each SLP return the operations in the genererators to get a matrix as it follows
               # 1 coeff * Z^(2k)
               # 0 1

               Append(slps,[slp]); # Append slp to array os slps

               k := k + 1;

               # Here we'll start getting slp products. Why does it work?
               # 1 coeff_1 * Z^(2k)  *  1 coeff_2 * Z^(2(k+1))  =>    1 (coeff_1 * Z^(2k) + coeff_2 * Z^(2(k+1)))
               # 0 1                    0 1                           0 1
               # Since we're assuming that z^0 , z^2, ... z^2(d-1) is a basis we get that
               # m[1][2] = coeff_1 * Z^(2k) + coeff_2 * Z^(2(k+1)) + ... + coeff_d-1 * Z^(2(k+d-1))
               if k = 2 then
                    SLP := ProductOfStraightLinePrograms(slps[k-1],slps[k]);
               elif k >= 3 then
                    SLP := ProductOfStraightLinePrograms(SLP,slps[k]);
               else
                    SLP := slps[1];
               fi;
          od;
          return SLP;
     # Second case: [[1,a],[0,1]]
     elif m[1][1] = One(f) and m[1][2] = Zero(f) and m[2][2] = One(f) then
          # For that we're assuming that z^0, z^2 , ... , z^(d-1) is a basis,
          # therefore a = a_0 z^0 + a_2 z^2 + ...

          # Creating the basis
          basis := [];
          k := 0;
          while k < d do
               Append(basis,[z^(2*k)]);
               k := k + 1;
          od;
          basis := Basis(f,basis);

          # Returning the coefficients of m[2][1] in that base
          coeff := Coefficients(basis,m[2][1]);

          # Computing slps to produce the product at the end
          slps := [];
          k := 0;
          while k < d do
               slp := StraightLineProgram([ [[3,k,2,1,3,-k],4], [[4,IntFFE(coeff[k+1])],5], [1,1,5,-1,1,-1] ],3);
               # Each SLP return the operations in the genererators to get a matrix as it follows
               # 1 coeff * Z^(2k)
               # 0 1

               Append(slps,[slp]); # Append slp to array os slps

               k := k + 1;

               # Here we'll start getting slp products. Why does it work?
               # 1                0  *  1                    0  =>    1                                         0
               # coeff_1 * Z^(2k) 1     coeff_2 * Z^(2(k+1)) 1        (coeff_1 * Z^(2k) + coeff_2 * Z^(2(k+1))) 1
               # Since we're assuming that z^0 , z^2, ... z^2(d-1) is a basis we get that
               # m[2][1] = coeff_1 * Z^(2k) + coeff_2 * Z^(2(k+1)) + ... + coeff_d-1 * Z^(2(k+d-1))
               if k = 2 then
                    SLP := ProductOfStraightLinePrograms(slps[k-1],slps[k]);
               elif k >= 3 then
                    SLP := ProductOfStraightLinePrograms(SLP,slps[k]);
               else
                    SLP := slps[1];
               fi;
          od;
          return SLP;
     # Third case: [[a,b],[c,d]] with b != 0
     elif m[1][2] <> Zero(f) then
          m3:=[[1,0],[(m[2][1]+m[2][2]*(One(f)-m[1][1])/m[1][2]),1]]*One(f);
          m2:=[[1,m[1][2]],[0,1]]*One(f);
          m1:=[[1,0],[(m[1][1]-One(f))/m[1][2],1]]*One(f);

          slp:=ProductOfStraightLinePrograms(Decompose2(m2,q),Decompose2(m1,q));
          SLP:=ProductOfStraightLinePrograms(Decompose2(m3,q),slp);

          return SLP;
     # Fourth case: [[a,b],[c,d]] with b == 0
     else
          m1:=[[0,1],[-1,0]]*One(f)^0;
          m:=m*m1;
          slp:=ProductOfStraightLinePrograms(Decompose2(m,q),Decompose2(Inverse(m1),q));
          return slp;
     fi;
end;

ClassicalStandardGenerators := function(m)
     local z,q,S,T,D;
     q := Size(DefaultFieldOfMatrix(m));
     z := Z(q);
     S := [[0,1],[-1,0]]*z^0;
     T := [[1,1],[0,1]]*z^0;
     D := [[z,0],[0,z^-1]]*z^0;
     return [S,T,D];
end;

TestCases := function(l)
     local i, m, q, a, good, bad;

     i := 0;
     q := Primes(Random([1..20]))^Random([1..10]);

     good := 0;
     bad := 0;

     while i < l do
          m := Random(SL(2,q));
          r := Decompose2(m,q);
          gens := ClassicalStandardGenerators(m);
          a := ResultOfStraightLineProgram(r,gens);

          if (a = m) then
               good := good + 1;
          else
               bad := bad + 1;
          fi;

          i := i + 1;

     od;
     
     Print("q = ",q,"\n");
     Print("In a total of ",l," test cases, ",good, " were good and ",bad," were bad");
end;
