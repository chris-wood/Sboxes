AssertAttribute(FldFin, "PowerPrinting", false);  
              // use polynomial basis for printing finite field elements
print "INTERPOLATION TEST";
print "";
print "";

// Polynomial ring over GF(2)
F:=GF(2);
pol2<X>:=PolynomialRing(F);    // Polynomial ring over F

// Build GF(2^8)
Q:=X^8+X^4+X^3+X+1;                     // the polynomial defining GF(2^n)
F256<x>:=ext<F | Q>;

// Build up the input and output map
Input:=[];
Output:=[];

// TODO: how to perform the affine transformation?
// TODO: figure out how to do matrix multiplication!
// b = 0x63 = {01100011} = x^6 + x^5 + x + 1

// Create the affine transformation matrix
affine:=Matrix(GF(2),8,8,
    [1,0,0,0,1,1,1,1,
    1,1,0,0,0,1,1,1,
    1,1,1,0,0,0,1,1,
    1,1,1,1,0,0,0,1,
    1,1,1,1,1,0,0,0,
    0,1,1,1,1,1,0,0,
    0,0,1,1,1,1,1,0,
    0,0,0,1,1,1,1,1]);
constant := x^6 + x^5 + x + 1;

// Generate the mapping for interpolation
for e in F256 do
  Input:=Append(Input,e);

  if e eq 0 then
    s:=ElementToSequence(e);
  else
    s:=ElementToSequence(e^(-1));
  end if;

  // perform the matrix computation
  // print "asdasd";
  // s;
  v:=Transpose(Matrix([s]));
  prod:=affine*v;

  // Transform back to the output (we don't transpose...)
  es:=[prod[1][1], prod[2][1], prod[3][1], prod[4][1], prod[5][1], prod[6][1], prod[7][1], prod[8][1]];
  // es;
  elem:=SequenceToElement(es, F256) + constant;
  // e;
  // elem;

  Output:=Append(Output, elem);
end for;

// Just an example...
Fx:=Interpolation(Input,Output);
Fx;
Degree(Fx);
Evaluate(Fx,0);
Evaluate(Fx,x^5+x+1);
Evaluate(Fx,x^7+x^6+x^5+x^4+x^3+x^2+x+1);

//
// For generating the interpolation polynomial...
// Interpolation(I, V) : [ RngElt ], [ RngElt ] -> RngUPolElt
//
// DESCRIPTION:
//
// This function finds a univariate polynomial that evaluates to the values V in the interpolation 
// points I. Let K be a field and n>0 an integer; given sequences I and V, both consisting of n elements of K, 
// return the unique univariate polynomial p over K of degree less than n such that p(I[i]) = V[i] for each 1≤i≤n.
//