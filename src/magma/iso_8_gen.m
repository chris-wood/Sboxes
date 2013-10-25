n:=2; m:=2;
k:=n*m;
AssertAttribute(FldFin, "PowerPrinting", false);  
              // use polynomial basis for printing finite field elements
//SetLogFile("C:/KUL/doctoraat/RSA_AES/uitvoeri");          // save all output in file "uitvoer"
F:=GF(2);
pol2<Y>:=PolynomialRing(F);    // Polynomial ring over F
Q:=Y^2+Y+1;                     // the polynomial defining GF(2^n)
F4<y>:=ext<F | Q>;             // y is a root of Q, generating GF(2^n)
pol4<X>:=PolynomialRing(F4); // Polynomial ring over GF(2^n)

phi:=y;
P:=X^2 + phi*X + 1;       // the polynomial defining GF((2^n)^m)
F16<x>:=ext<F4 | P>;           // x is a root of P, generating GF(2^k)

Pol16<Z>:=PolynomialRing(F16); // Polynomial ring over GF(2^k)

lambda:=x;
R:=Z^2+lambda*Z+1;
F256<z>:=ext<F16 | R>;           // x is a root of P, generating GF(2^k)

Pol256<T>:=PolynomialRing(F256);


print(Pol256);
print(F);

S:=T^8 + T^4 + T^3 + T + 1;   

for beta1 in Roots(S) do
  beta:=beta1[1];

  n:=8;
  result:=Transpose(Matrix([
  [Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[1],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[1],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[1],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[1]] : i in [0..7]]));

  nt := Matrix([
  [Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[1],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[1],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[1],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[2],
  Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[1]] : i in [0..7]]);

  result_inv:=result^-1;
  print("one basis change");
  print(beta);
  print(result);
  print "";
  print(result_inv);
  print "";
  print(nt);
  print "";
end for;












