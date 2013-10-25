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

phi:=y; // satoh uses phi = {10}_2
P:=X^2 + X + phi;       // the polynomial defining GF((2^n)^m)
F16<x>:=ext<F4 | P>;           // x is a root of P, generating GF(2^k)

Pol16<Z>:=PolynomialRing(F16); // Polynomial ring over GF(2^k)

lambda:=x^3 + x^2; // satoh uses {1100}_2
R:=Z^2+Z+lambda;
F256<z>:=ext<F16 | R>;           // x is a root of P, generating GF(2^k)

// Extend GF(2^8) by degree 2
Pol256<T>:=PolynomialRing(F256);
rho:=z^3 + z; // we use {00001010}_2, 0x0A
S:=T^2 + T + rho;
F65536<v>:=ext<F256 | S>;

print(F65536);

// GF(2^16)/GF(2)
Pol65536<U>:=PolynomialRing(F65536);
T:=U^16 + U^5 + U^3 + U^2 + 1;

// FIND THE ROOTS!
for beta1 in Roots(T) do
  beta:=beta1[1];

  n:=16;
  result:=Transpose(Matrix([[
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[2],

    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[2]] : i in [0..15]]));

    n:=16;
  actualResult:=Transpose(Matrix([[
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[1])[1],

    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[1])[1]] : i in [0..15]]));

  nt := Matrix([[
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[2])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[2])[1])[1])[1],

    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[2])[1])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[2])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[2])[1],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[1])[2],
    Eltseq(Eltseq(Eltseq(Eltseq(beta^(n-i-1))[1])[1])[1])[1]] : i in [0..15]]);

  result_ones:=NumberOfNonZeroEntries(result);
  result_inv:=result^-1;
  print("Result and Inverse");
  print(beta);
  print(result);
  print("");
  print(result_inv);
  print("");
  print(actualResult);
  print("");
  print(actualResult^-1);
  print("");
  print(nt);
  print("");
  print(nt^-1);
end for;












