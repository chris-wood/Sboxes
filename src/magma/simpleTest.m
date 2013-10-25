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

// Build basis elems given the polynomial roots
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
      // print "starting ", pr, e1, Evaluate(Q, e1);
      if Evaluate(Q, e1) eq 0 then
        qRoots:=Append(qRoots,e1);
      end if;
    end for;
    // pr;
    // Q;
    // qRoots;
  end if;
end for;

Q:=X^2+X+w;
F16<x>:=ext<F4 | Q>;
Pol16<Y>:=PolynomialRing(F16);
R:=Y^2+Y+w*x + w;
F256<y>:=ext<F16|R>; 
for e in F256 do
  if Evaluate(R,e) eq 0 then
    print "start";
    // w+1;
    // w*x + w + 1;
    b:=(w+1)*(w*x + w + 1)*e;
    basisElems(w+1, w*x + w + 1, e, Ffull);
    // IsNormal(b, F16);
    toBinary8(b^(2^0));
    toBinary8(b^(2^1));
    toBinary8(b^(2^2));
    toBinary8(b^(2^3));
    toBinary8(b^(2^4));
    toBinary8(b^(2^5));
    toBinary8(b^(2^6));
    toBinary8(b^(2^7));
  end if;
end for;  

