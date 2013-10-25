// File: rootGen.m
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
basisElems := function(pr, qr, rr, field)
  basis:=[];
  basis:=Append(basis,SequenceToElement(toBinary8(pr*qr*rr), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr^2*qr*rr), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr*qr^4*rr), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr^2*qr^4*rr), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr*qr*rr^16), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr^2*qr*rr^16), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr*qr^4*rr^16), field));
  basis:=Append(basis,SequenceToElement(toBinary8(pr^2*qr^4*rr^16), field));
  return basis;
end function;

//////////////////////////////
// START ROOT/POLY GENERATION
//////////////////////////////
// Fb:=GF(2);
// p<D>:=PolynomialRing(Fb);
// G:=D^8+D^4+D^3+D+1;
// F:=ext<Fb|G>;
F:=GF(2);
pol2<W>:=PolynomialRing(F); // polynomial ring over F
P:=W^2+W+1;                 // the only irreducible polynomial for GF(2^2)
F4<w>:=ext<F | P>;          // create GF(2^2)/GF(2)

Pfull:=W^8+W^4+W^3+W+1;
Ffull<d>:=ext<F|Pfull>;

// Generate all roots of P
pRoots:=[];
for e in F4 do
  if Evaluate(P, e) eq 0 then
    pRoots:=Append(pRoots,e);
  end if;
end for;

// Now for each root, plug into the next polynomial higher up and 
// see if we have an irreducible, and if so, find the roots.
// Continue this process higher and higher...
IPs:=[];
for pr in pRoots do
  pol4<X>:=PolynomialRing(F4);
  Q:=X^2+X+pr;
  if IsIrreducible(Q) then
    F16<x>:=ext<F4 | Q>;
    qRoots:=[];
    for e1 in F16 do
      if Evaluate(Q, e1) eq 0 then
        qRoots:=Append(qRoots,e1);
      end if;
    end for;
    for qr in qRoots do
      pol16<Y>:=PolynomialRing(F16);
      R:=Y^2+Y+qr; 
      if IsIrreducible(R) then
        F256<y>:=ext<F16 | R>;
        rRoots:=[];
        for e2 in F256 do
          if Evaluate(R, e2) eq 0 then
            rRoots:=Append(rRoots, e2);
          end if;
        end for;

        print "FUCK YES";

        // for rr in rRoots do
        //   print pr,qr,rr;
        //   print P;
        //   print Q;
        //   print R;
        // end for;
      end if;

      // Try the other powers of q, as they might make the polynomial irreducible
      qRootsPowers:=rootPowers(qr, 16, qRoots);
      polys:=[];
      for qr2 in qRootsPowers do
        pol16<Y>:=PolynomialRing(F16);
        R:=Y^2+Y+qr2;

        // qr2;
        if IsIrreducible(R) then
          // print "FUCKYES";

          F256<y>:=ext<F16 | R>;
          rRoots:=[];
          for e2 in F256 do
            if Evaluate(R, e2) eq 0 then
              rRoots:=Append(rRoots, e2);
            end if;
          end for;

          // Testing the output to make sure it matches Canright's results...
          for rr in rRoots do
            // print pr,qr,qr2,rr;
            // print ("displaying...");
            print pr, qr, rr;
            print P;
            print Q;
            print R;

            Pol256<T>:=PolynomialRing(F256);
            S:=W^8 + W^4 + W^3 + W + 1;
            // S:=T^8+T^4+T^3+T+1;
            Freal<t>:=ext<F|S>;

            b1:= pr * qr * rr; // multiply the basis elements together and see if we get what canright had...

            // Eltseq(pr);
            // b1;
            // toBinary8(b1);
            bin:=SequenceToElement(toBinary8(b1), Freal);
            // IsNormal(b1);
            es:=basisElems(pr, qr, rr, Freal);
            // es;
            binStrings:=[];
            for i in [1..8] do
              binStrings:=Append(binStrings, Eltseq(es[i]));
            end for;
            // binStrings;

            matrix:=basisMatrix8(binStrings);
            matrix;
            print "";
            Transpose(matrix);
            print "";
            Transpose(matrix^(-1));
            NumberOfNonZeroEntries(matrix);
          end for;
        end if;
      end for;
    end for;
  end if;
end for;


