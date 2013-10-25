// File: smallNormalVerify.m
// Author: Christopher Wood
// Generate all irreducible polynomials and their roots for 
// polynomials of interest 

// Uncomment for standard basis printing
AssertAttribute(FldFin, "PowerPrinting", false);  
SetQuitOnError(true);

//////////////////////////////
// rootPowers - generate all powers for a given root
//////////////////////////////
rootPowers := function(elem, maxDeg, roots)
  genElems:=[]; 
  for i := 0 to maxDeg do
    e := elem^i;
    // if Position(roots, e) eq 0 then
    if e notin roots then
      genElems:=Append(genElems, elem^i);
    end if; 
  end for;
  return genElems;
end function;

//////////////////////////////
// weakAppend - Try to append an element to a sequence if it's not already present
//////////////////////////////
weakAppend := function(list, elem)
  if elem notin list then
    return Append(list, elem);
  end if;
  return list;
end function;

//////////////////////////////
// toBinary - Print the binary output of an element in GF(2^8)/GF(2^4)/GF(2^2)/GF(2)
//////////////////////////////
toBinary8 := function(elem)
  result:=
    [Eltseq(Eltseq(Eltseq(elem)[1])[1])[1]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[1])[1])[2]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[1])[2])[1]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[1])[2])[2]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[2])[1])[1]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[2])[1])[2]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[2])[2])[1]] cat 
    [Eltseq(Eltseq(Eltseq(elem)[2])[2])[2]];
  return result;
end function;

basisMatrix8 := function(binStrings)
  m:=Matrix(GF(2),8,8,
    [
      binStrings[8][8], binStrings[7][8], binStrings[6][8], binStrings[5][8], binStrings[4][8], binStrings[3][8], binStrings[2][8], binStrings[1][8],
  binStrings[8][7], binStrings[7][7], binStrings[6][7], binStrings[5][7], binStrings[4][7], binStrings[3][7], binStrings[2][7], binStrings[1][7],
  binStrings[8][6], binStrings[7][6], binStrings[6][6], binStrings[5][6], binStrings[4][6], binStrings[3][6], binStrings[2][6], binStrings[1][6],
  binStrings[8][5], binStrings[7][5], binStrings[6][5], binStrings[5][5], binStrings[4][5], binStrings[3][5], binStrings[2][5], binStrings[1][5],
  binStrings[8][4], binStrings[7][4], binStrings[6][4], binStrings[5][4], binStrings[4][4], binStrings[3][4], binStrings[2][4], binStrings[1][4],
  binStrings[8][3], binStrings[7][3], binStrings[6][3], binStrings[5][3], binStrings[4][3], binStrings[3][3], binStrings[2][3], binStrings[1][3],
  binStrings[8][2], binStrings[7][2], binStrings[6][2], binStrings[5][2], binStrings[4][2], binStrings[3][2], binStrings[2][2], binStrings[1][2],
  binStrings[8][1], binStrings[7][1], binStrings[6][1], binStrings[5][1], binStrings[4][1], binStrings[3][1], binStrings[2][1], binStrings[1][1]
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
basisElemsPoly := function(pr, qr)
  basis:=[];
  basis:=Append(basis,qr^0);
  basis:=Append(basis,pr*(qr^0));
  basis:=Append(basis,(pr^0) * qr);
  basis:=Append(basis,pr * qr);
  return basis;
end function;

basisElemsPoly8 := function(pr, qr, rr)
  basis:=[];
  basis:=Append(basis,rr^0);
  basis:=Append(basis,pr*(rr^0));
  basis:=Append(basis,(pr^0) * qr * (rr^0));
  basis:=Append(basis,pr * qr * (rr^0));
  basis:=Append(basis,rr^1);
  basis:=Append(basis,pr*(rr^1));
  basis:=Append(basis,(pr^0) * qr * (rr^1));
  basis:=Append(basis,pr * qr * (rr^1));
  return basis;
end function;

mmult := function(m, s)
  v:=Transpose(Matrix([s])); // create a matrix...
  v;
  m;
  prod:=m*v;
  es:=[prod[1][1], prod[2][1], prod[3][1], prod[4][1]];
  return es;
end function;

////////////////////
//// START HERE ////
////////////////////

F:=GF(2);
pol2<W>:=PolynomialRing(F); // polynomial ring over F
P:=W^2+W+1;                 // the only irreducible polynomial for GF(2^2)
F4<w>:=ext<F | P>;          // create GF(2^2)/GF(2)
poly4<Y>:=PolynomialRing(F4);
S:=W^4 + W + 1;
Fr<z>:=ext<F | S>;

// Generate all roots of P
pRoots:=[];
for e in F4 do
  if Evaluate(P, e) eq 0 then
    pRoots:=Append(pRoots,e);
  end if;
end for;

try 
  Embed(F4, Fr); 
catch  e
  print("Embedding already done. No harm no foul.");
end try;

// Now for each root, plug into the next polynomial higher up and 
// see if we have an irreducible, and if so, find the roots.
// Continue this process higher and higher...
IPs:=[];

// sort the proots
pRoots := Sort(pRoots);

for e in F4 do
  pol4<X>:=PolynomialRing(F4);
  Q:=X^2+X+e;
  if IsIrreducible(Q) then
    print "IS IRREDUCIBLE";
    print e;
  end if;
end for;

for pr in pRoots do
  pol4<X>:=PolynomialRing(F4);
  Q:=X^2+X+pr;

  if IsIrreducible(Q) then
    F16<x>:=ext<F4 | Q>;

    try 
      Embed(F16, Fr); 
    catch  e
      print("Embedding already done. No harm no foul.");
    end try;

    qRoots:=[];
    for e1 in F16 do
      if Evaluate(Q, e1) eq 0 then
        qRoots:=Append(qRoots,e1);
      end if;
    end for;

    for qr in qRoots do
      // Coerce the proot into an element in Freal
      pr;
      qr;
      a:=Fr!pr;
      b:=Fr!qr;
      a;
      b;

      // Printing the polynomials
      P;
      Q;

      // for e in Fr do
      //   if Evaluate(S, e) eq 0 then
      //     e;
      //     print("hells yeah");
      //   end if;
      // end for;

      // Create the basis elements
      es:=basisElemsPoly(Fr!pr, Fr!qr);
      es;
      binStrings:=[];
      for i in [1..4] do
        binStrings:=Append(binStrings, Eltseq(es[i]));
      end for;

      // Now build and output the basis change matrices
      matrix:=basisMatrix4(binStrings);

      // Row order (for VHDL and eyeballs)
      m:=matrix;
      m;
      print "";
      mi:=matrix^-1;
      mi;
      print "";

      // Perform our own checks here...
      // SetPrimitiveElement(F16, qr);
      // for e in Fr do
      //   print "original mapping", e;
      //   mapped := mmult(mi, Reverse(Eltseq(e)));
      //   mapped_l := Seqelt([mapped[4], mapped[3]], F4);
      //   mapped_h := Seqelt([mapped[2], mapped[1]], F4);
      //   together := [mapped_l, mapped_h];
      //   ecomp := Seqelt(together, F16);
      //   mapped;
      //   mapped_l;
      //   mapped_h;

      //   // Compute the addition and square results for both elements
      //   comp_square := ecomp^2;
      //   comp_added := ecomp+1;
      //   ecomp;
      //   comp_square;
      //   comp_added;
      //   square := Eltseq(e^2);
      //   added := Eltseq(e+1);

      //   // map both back and check 
      //   seq_comp_square := [Eltseq(Eltseq(comp_square)[1])[1], Eltseq(Eltseq(comp_square)[1])[2], Eltseq(Eltseq(comp_square)[2])[1], Eltseq(Eltseq(comp_square)[2])[2]];
      //   seq_comp_added := [Eltseq(Eltseq(comp_added)[1])[1], Eltseq(Eltseq(comp_added)[1])[2], Eltseq(Eltseq(comp_added)[2])[1], Eltseq(Eltseq(comp_added)[2])[2]];
      //   mappedSquare:=mmult(m, seq_comp_square);
      //   mappedAdded:=mmult(m, seq_comp_added);
        
      //   // print "squared", e;
      //   // square;
      //   // mappedSquare;

      //   // print "added", e;
      //   // added;
      //   // mappedAdded;
      // end for;  
      
      // Column order (for C code)
      // Transpose(matrix);
      // print "";
      // Transpose(matrix^(-1));
      // NumberOfNonZeroEntries(matrix);
    end for;
  end if;
end for;

// Canright's test vector set in the paper...
S:=W^8 + W^4 + W^3 + W + 1;
Fr<z>:=ext<F | S>;
es:=basisElemsPoly8(z^7+z^5+z^4+z^3+z^2+1, z^6+z^4+z^3+z^2, z^7+z^5+z^5+z^4+z^3+z^2+z+1); //0xBD, 0x5C, 0xFF
es;
binStrings:=[];
for i in [1..8] do
  binStrings:=Append(binStrings, Eltseq(es[i]));
end for;
binStrings;

// Now build and output the basis change matrices
matrix:=basisMatrix8(binStrings);

// Row order (for VHDL and eyeballs)
matrix;
print "";
matrix^-1;
print "";
