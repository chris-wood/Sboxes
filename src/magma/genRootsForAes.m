AssertAttribute(FldFin, "PowerPrinting", false);  
SetQuitOnError(true);

basisMatrix4 := function(binStrings)
  m:=Matrix(GF(2),4,4,
    [
      binStrings[4][4], binStrings[3][4], binStrings[2][4], binStrings[1][4],
      binStrings[4][3], binStrings[3][3], binStrings[2][3], binStrings[1][3],
      binStrings[4][2], binStrings[3][2], binStrings[2][2], binStrings[1][2],
      binStrings[4][1], binStrings[3][1], binStrings[2][1], binStrings[1][1]
    ]);

  // binStrings[8][8], binStrings[7][8], binStrings[6][8], binStrings[5][8], binStrings[4][8], binStrings[3][8], binStrings[2][8], binStrings[1][8];
  // binStrings[8][7], binStrings[7][7], binStrings[6][7], binStrings[5][7], binStrings[4][7], binStrings[3][7], binStrings[2][7], binStrings[1][7];
  // binStrings[8][6], binStrings[7][6], binStrings[6][6], binStrings[5][6], binStrings[4][6], binStrings[3][6], binStrings[2][6], binStrings[1][6];
  // binStrings[8][5], binStrings[7][5], binStrings[6][5], binStrings[5][5], binStrings[4][5], binStrings[3][5], binStrings[2][5], binStrings[1][5];
  // binStrings[8][4], binStrings[7][4], binStrings[6][4], binStrings[5][4], binStrings[4][4], binStrings[3][4], binStrings[2][4], binStrings[1][4];
  // binStrings[8][3], binStrings[7][3], binStrings[6][3], binStrings[5][3], binStrings[4][3], binStrings[3][3], binStrings[2][3], binStrings[1][3];
  // binStrings[8][2], binStrings[7][2], binStrings[6][2], binStrings[5][2], binStrings[4][2], binStrings[3][2], binStrings[2][2], binStrings[1][2];
  // binStrings[8][1], binStrings[7][1], binStrings[6][1], binStrings[5][1], binStrings[4][1], binStrings[3][1], binStrings[2][1], binStrings[1][1];
  return m;
end function;

//////////////////////////////
// toBinary - Print the binary output of an element in GF(2^8)/GF(2^4)/GF(2^2)/GF(2)
//////////////////////////////
toBinary4 := function(elem)
  result:=
    [Eltseq(Eltseq(elem)[1])[1]] cat 
    [Eltseq(Eltseq(elem)[1])[2]] cat 
    [Eltseq(Eltseq(elem)[2])[1]] cat 
    [Eltseq(Eltseq(elem)[2])[2]];
  return result;
end function;

// Build basis binStrings given the polynomial roots
basisElems := function(pr, qr, field)
  basis:=[];
  basis:=Append(basis,SequenceToElement(toBinary4(pr*qr), field));
  basis:=Append(basis,SequenceToElement(toBinary4(pr^2 * qr), field));
  basis:=Append(basis,SequenceToElement(toBinary4(pr * qr^4), field));
  basis:=Append(basis,SequenceToElement(toBinary4(pr^2 * qr^4), field));
  return basis;
end function;

// Build basis binStrings given the polynomial roots
basisElemsPoly := function(pr, qr, field)
  basis:=[];
  basis:=Append(basis,qr^0);
  basis:=Append(basis,pr*(qr^0));
  basis:=Append(basis,pr^0 * qr);
  basis:=Append(basis,pr * qr);
  return basis;
end function;

/// CODE STARTS HERE ///

F:=GF(2);
pol2<W>:=PolynomialRing(F); // polynomial ring over F
S1:=W^2+W+1;
F1<a>:=ext<F|S1>;
S:=W^4 + W + 1;
Fr<z>:=ext<F | S>;
poly256<Z>:=PolynomialRing(Fr);
T:=Z^2+Z+1;

// new way...
// if IsIrreducible(T) then
	T;
  for r1 in Roots(T) do
    // r1[1];
    // Eltseq(r1[1]);
    A:=Z^2+Z+r1[1];
    // A;
	
	try 
		Embed(F1, Fr, r1[1]); 
	catch  e
		print("Embedding already done");
	end try;

    for r2 in Roots(A) do
      // r2[1];
 		// Eltseq(r2[1]);
 		r1[1] * r2[1];
 		r2[1];
 		Eltseq(r2[1], F1);
 		elt<F1 | Eltseq(r2[1], F1)>;


 		// hom<Fr -> F1 | Generator(Fr)>;

 	    //   es:=basisElemsPoly(r1[1], r2[1], Fr);
      // es;
      // binStrings:=[];
      // for i in [1..4] do
      //   binStrings:=Append(binStrings, Eltseq(es[i]));
      // end for;
      // // binStrings;

      // matrix:=basisMatrix4(binStrings);
      // matrix;
      // print "";
      // matrix^-1;
      // print "";
      // Transpose(matrix);
      // print "";
      // Transpose(matrix^(-1));
      // NumberOfNonZeroEntries(matrix);
    end for;
  end for;
// end if;

