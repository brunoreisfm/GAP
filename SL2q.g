Decompose2 := function(m,q)
     local f, d, z, k, basis, coeff, slps, slp, SLP;

     f := GF(q); #DefaultFieldOfMatrix(m); # Finite field of the matrix
     #q := Size(f); # q = p^d which is the size of the finite field
     d := Size(GaloisGroup(f)); # d from q = p^d
     z := Z(q); # Primitive element of the finite field

     if IntFFE(m[1][1]) = 1 and IntFFE(m[2][2]) = 1 and IntFFE(m[2][1]) = 0 then
          # First case: [[1,a],[0,1]]
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
     else
          Display("\nThat decomposition has not been implemented yet\n");
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
