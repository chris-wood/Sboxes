AssertAttribute(FldFin, "PowerPrinting", false);

buildRoot4Row := function(element, subfield) 
	if element eq 0 then // ADDED
		return [0, 0, 0, 0];
	end if;
	if element eq 1 then
		return [0,0,0,1];
	end if;
	elem := Eltseq(element, subfield);

	// no w term, so insert two leading zeroes
	// handles the case if v or v+1
	if #elem eq 1 then
		tmp := Eltseq(elem[1]);
		row := [0,0,tmp[2],tmp[1]];
		return row;
	end if;

	tmp := [Eltseq(elem[2])[2], Eltseq(elem[2])[1], Eltseq(elem[1])[2], Eltseq(elem[1])[1]];
	return tmp;
end function;

changeSigmaRoot := function(pr1, pr2, F4, sigma) // sigma in GF(2^2)
	if pr2 eq 0 then // ADDED
		r1 := [0, 0];
	end if;
	if pr2 eq 1 then
		r1 := [0, 1];
	else
		r1 := Reverse(Eltseq(pr2));
	end if;
	if pr1 eq 0 then // ADDED
		r2 := [0, 0];
	end if;
	if pr1 eq 1 then
		r2 := [0, 1];
	else
		r2 := Reverse(Eltseq(pr1));
	end if;
	M := Matrix(GF(2), [ r1, r2 ]);

	if sigma eq 0 then
		evector := [0,0];
	elif sigma eq 1 then
		evector := [0,1];
	else
		evector := Reverse(Eltseq(sigma));
	end if;

	ev := Transpose(Matrix(GF(2), [ evector ]));
	prod := Transpose(M)^(-1) * ev;
	newSigma := Seqelt([prod[2][1], prod[1][1]], F4);
	return newSigma;
end function;

changePiRoot := function(pr1, pr2, qr1, qr2, F4, F16, pi) // pi in GF(2^4)/GF(2^2)
	q1p1 := pr1 * qr1;
	q2p1 := pr1 * qr2;
	q1p2 := pr2 * qr1;
	q2p2 := pr2 * qr2;

	r1 := buildRoot4Row(q2p2, F4);
	r2 := buildRoot4Row(q2p1, F4);
	r3 := buildRoot4Row(q1p2, F4);
	r4 := buildRoot4Row(q1p1, F4);

	// r1;
	// r2;
	// r3;
	// r4;

	M := Matrix(GF(2), [  r1, r2, r3, r4 ]);
	evector := [Eltseq(Eltseq(pi)[2])[2], Eltseq(Eltseq(pi)[2])[1], Eltseq(Eltseq(pi)[1])[2], Eltseq(Eltseq(pi)[1])[1]];
	ev := Transpose(Matrix(GF(2), [evector])); // correct
	prod := Transpose(M)^(-1) * ev;
	newPi := Seqelt([Seqelt([prod[4][1],prod[3][1]], F4), Seqelt([prod[2][1],prod[1][1]], F4)], F16);
	return newPi;
end function;

toBinary4 := function(elem, F4) // this puts upper bits in the leftmost index
	return [Eltseq(Eltseq(elem, F4)[2])[2], Eltseq(Eltseq(elem, F4)[2])[1], Eltseq(Eltseq(elem, F4)[1])[2], Eltseq(Eltseq(elem, F4)[1])[1]];
end function;

// Generate all tuples of basis elements (p, r, q, s)
polyGenDegree2_16 := function()
  F:=GF(2);
  pol2<V>:=PolynomialRing(F); // polynomial ring over F
  P:=V^2+V+1;                 // the only irreducible polynomial for GF(2^2)
  F4<v>:=ext<F | P>;          // create GF(2^2)/GF(2)

  // Generate all roots of P
  pRoots:=[];
  for e in F4 do
    if Evaluate(P, e) eq 0 then
      pRoots:=Append(pRoots,e);
    end if;
  end for;

  // Storage containers for the irreducible polynomials
  pIps:=Append([], P);
  qIps:=[];
  rIps:=[];

  // Find all elements in GF(2^2) that make q(x) irreducible
  for pr in F4 do
    pol4<W>:=PolynomialRing(F4);
    Q:=W^2+W+pr;
    if IsIrreducible(Q) and Q notin qIps then
      F16<w>:=ext<F4 | Q>;

      // Find the roots of q(x) - GF(2^4)/GF(2^2)
      qRoots:=[];
      for e1 in F16 do
        if Evaluate(Q, e1) eq 0 then
          qRoots:=Append(qRoots,e1); // there will always be two roots for a degree two extension
        end if;
      end for;

      //field := [1, v, v + 1, w, w + 1, w + v, w + v + 1, v*w, v*w + 1, v*w + v, v*w + v + 1, (v+1)*w, (v+1)*w + 1, (v+1)*w + v, (v+1)*w + v + 1];
      //#field;

      // Find all elements in GF(2^4)/GF(2^2) that make r(y) irreducible
	1, pRoots[1], 1, qRoots[1], pr;
      for e in F16 do
     if e ne 0 then
		newe := changePiRoot(F4!1, pRoots[1], F4!1, qRoots[1], F4, F16, e);
		newe;
		toBinary4(newe, F4);
		inv := e^-1;
		newei := changePiRoot(F4!1, pRoots[1], F4!1, qRoots[1], F4, F16, inv);
		newei;
		toBinary4(newei, F4);
	end if;
	
      end for;

	1, pRoots[1], 1, qRoots[2], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(F4!1, pRoots[1], F4!1, qRoots[2], F4, F16, e);
	toBinary4(newe, F4);
		inv := e^-1;
	newei := changePiRoot(F4!1, pRoots[1], F4!1, qRoots[2], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	1, pRoots[1], qRoots[1], qRoots[2], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(F4!1, pRoots[1], qRoots[1], qRoots[2], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(F4!1, pRoots[1], qRoots[1], qRoots[2], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	1, pRoots[2], 1, qRoots[1], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(F4!1, pRoots[2], F4!1, qRoots[1], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(F4!1, pRoots[2], F4!1, qRoots[1], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	1, pRoots[2], 1, qRoots[2], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(F4!1, pRoots[2], F4!1, qRoots[2], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(F4!1, pRoots[2], F4!1, qRoots[2], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	1, pRoots[2], qRoots[1], qRoots[2], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(F4!1, pRoots[2], qRoots[1], qRoots[2], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(F4!1, pRoots[2], qRoots[1], qRoots[2], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	pRoots[1], pRoots[2], 1, qRoots[1], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(pRoots[1], pRoots[2], F4!1, qRoots[1], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(pRoots[1], pRoots[2], F4!1, qRoots[1], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	pRoots[1], pRoots[2], 1, qRoots[2], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(pRoots[1], pRoots[2], F4!1, qRoots[2], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(pRoots[1], pRoots[2], F4!1, qRoots[2], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

	pRoots[1], pRoots[2], qRoots[1], qRoots[2], pr;
      for e in F16 do
	if e ne 0 then
	newe := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, e);
	toBinary4(newe, F4);
	inv := e^-1;
	newei := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, inv);
	toBinary4(newei, F4);
	end if;
      end for;

    end if;
  end for;
  return [* pIps, qIps, rIps *];
end function;

// polyGenDegree2_16();
F:=GF(2);
pol2<V>:=PolynomialRing(F); // polynomial ring over F
polyGenDegree2_16();

// There's three irreducible polynomials for GF(2^4)

P:=V^4+V+1;                 // the only irreducible polynomial for GF(2^2)
F4<v>:=ext<F | P>;          // create GF(2^2)/GF(2)
P;
for e in F4 do
	if e ne 0 then
		Eltseq(e);
		Eltseq(e^-1);	
	else
		[0,0,0,0];
		[0,0,0,0];
	end if;
end for;

P:=V^4 + V^3 + 1;
F4<v>:=ext<F | P>;          // create GF(2^2)/GF(2)
P;
for e in F4 do
	if e ne 0 then
		Eltseq(e);
		Eltseq(e^-1);	
	else
		[0,0,0,0];
		[0,0,0,0];
	end if;
end for;

P:=V^4 + V^3 + V^2 + V + 1;
F4<v>:=ext<F | P>;          // create GF(2^2)/GF(2)
P;
for e in F4 do
	if e ne 0 then
		Eltseq(e);
		Eltseq(e^-1);	
	else
		[0,0,0,0];
		[0,0,0,0];
	end if;
end for;




