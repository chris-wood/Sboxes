AssertAttribute(FldFin, "PowerPrinting", false);

// CORRECT
buildRoot4Row := function(element, subfield) 
	if element eq 1 then
		return [0,0,0,1];
	end if;
	elem := Eltseq(element, subfield);

	if #elem eq 1 then
		tmp := Eltseq(elem[1]);
		row := [0,0,tmp[2],tmp[1]];
		return row;
	end if;

	tmp := [Eltseq(elem[2])[2], Eltseq(elem[2])[1], Eltseq(elem[1])[2], Eltseq(elem[1])[1]];
	return tmp;
end function;

// TODO: finish this (it'll be used for the gate counting code)
buildRoot8Row := function(element, subfield) 

	// element;

	if element eq 1 then
		return [0,0,0,0,0,0,0,1];
	end if;

	// element;
	// subfield;
	// element in subfield;
	elem := Eltseq(element, subfield);

	// no w term, so insert two leading zeroes
	if #elem eq 1 then
		tmp := Eltseq(elem[1]); // [2] refers to upper bits
		tmp1 := Eltseq(tmp)[1]; // lower
		tmp2 := Eltseq(tmp)[2]; // upper
		row := [0,0,0,0,Eltseq(tmp2)[2], Eltseq(tmp2)[1], Eltseq(tmp1)[2], Eltseq(tmp1)[1]];
		return row;
	end if;

	p4 := Eltseq(Eltseq(elem)[2])[2];
	p3 := Eltseq(Eltseq(elem)[2])[1];
	p2 := Eltseq(Eltseq(elem)[1])[2];
	p1 := Eltseq(Eltseq(elem)[1])[1];

	// element;
	// [ p4, p3, p2, p1 ];

	// tmp := [Eltseq(elem[2])[2], Eltseq(elem[2])[1], Eltseq(elem[1])[2], Eltseq(elem[1])[1]];

	tmp := [Eltseq(p4)[2], Eltseq(p4)[1], Eltseq(p3)[2], Eltseq(p3)[1], Eltseq(p2)[2], Eltseq(p2)[1], Eltseq(p1)[2], Eltseq(p1)[1]];

	return tmp;
end function;

// sigma in GF(2^2)
// CORRECT
changeSigmaRoot := function(pr1, pr2, F4, sigma) // sigma in GF(2^2)
	if pr2 eq 1 then
		r1 := [0, 1];
	else
		r1 := Reverse(Eltseq(pr2));
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

// pi in GF(2^4)/GF(2^2)
// CORRECT
changePiRoot := function(pr1, pr2, qr1, qr2, F4, F16, pi) // pi in GF(2^4)/GF(2^2)
	q1p1 := pr1 * qr1;
	q2p1 := pr1 * qr2;
	q1p2 := pr2 * qr1;
	q2p2 := pr2 * qr2;
	r1 := buildRoot4Row(q2p2, F4);
	r2 := buildRoot4Row(q2p1, F4);
	r3 := buildRoot4Row(q1p2, F4);
	r4 := buildRoot4Row(q1p1, F4);

	// [r1, r2, r3, r4];

	M := Matrix(GF(2), [  r1, r2, r3, r4 ]);
	evector := [Eltseq(Eltseq(pi)[2])[2], Eltseq(Eltseq(pi)[2])[1], Eltseq(Eltseq(pi)[1])[2], Eltseq(Eltseq(pi)[1])[1]];
	ev := Transpose(Matrix(GF(2), [evector])); // correct
	prod := Transpose(M)^(-1) * ev;
	newPi := Seqelt([Seqelt([prod[4][1],prod[3][1]], F4), Seqelt([prod[2][1],prod[1][1]], F4)], F16);
	return newPi;
end function;

// lambda in GF(2^8)/GF(2^4)/GF(2^2)
changeLambdaRoot := function(pr1, pr2, qr1, qr2, rr1, rr2, F4, F16, F256, lambda)
	p1q1r1 := pr1 * qr1 * rr1;
	p1q2r1 := pr1 * qr2 * rr1;
	p2q1r1 := pr2 * qr1 * rr1;
	p2q2r1 := pr2 * qr2 * rr1;
	p1q1r2 := pr1 * qr1 * rr2;
	p1q2r2 := pr1 * qr2 * rr2;
	p2q1r2 := pr2 * qr1 * rr2;
	p2q2r2 := pr2 * qr2 * rr2;
	r1 := buildRoot8Row(F256!p2q2r2, F16);
	r2 := buildRoot8Row(F256!p1q2r2, F16);
	r3 := buildRoot8Row(F256!p2q1r2, F16);
	r4 := buildRoot8Row(F256!p1q1r2, F16);
	r5 := buildRoot8Row(F256!p2q2r1, F16);
	r6 := buildRoot8Row(F256!p1q2r1, F16);
	r7 := buildRoot8Row(F256!p2q1r1, F16);
	r8 := buildRoot8Row(F256!p1q1r1, F16);

	// [r1, r2, r3, r4, r5, r6, r7, r8];

	// Create the inverse basis change matrix...
	M := Matrix(GF(2), [  r1, r2, r3, r4, r5, r6, r7, r8 ]);

	// Create the element sequences
	upper := Eltseq(lambda)[2];
	lower := Eltseq(lambda)[1];
	p4 := Eltseq(upper)[2];
	p3 := Eltseq(upper)[1];
	p2 := Eltseq(lower)[2];
	p1 := Eltseq(lower)[1];

	// evector := [Eltseq(Eltseq(pi)[2])[2], Eltseq(Eltseq(pi)[2])[1], Eltseq(Eltseq(pi)[1])[2], Eltseq(Eltseq(pi)[1])[1]];
	evector := [Eltseq(p4)[2], Eltseq(p4)[1], Eltseq(p3)[2], Eltseq(p3)[1], Eltseq(p2)[2], Eltseq(p2)[1], Eltseq(p1)[2], Eltseq(p1)[1]];
	// lambda;
	// evector;
	ev := Transpose(Matrix(GF(2), [evector])); 
	// Transpose(M)^-1;
	prod := Transpose(M)^(-1) * ev;

	newLambda := Seqelt( 
		[
			Seqelt( [ Seqelt([prod[8][1],prod[7][1]], F4), Seqelt([prod[6][1],prod[5][1]], F4) ], F16), 
			Seqelt( [ Seqelt([prod[4][1],prod[3][1]], F4), Seqelt([prod[2][1],prod[1][1]], F4) ], F16)
		], F256);
	return newLambda;
end function;

// Generate all tuples of basis elements (p, r, q, s)
polyGenDegree2_16 := function()
  F:=GF(2);
  pol2<V>:=PolynomialRing(F); // polynomial ring over F
  P:=V^2+V+1;                 // the only irreducible polynomial for GF(2^2)
  F4<v>:=ext<F | P>;          // create GF(2^2)/GF(2)

  // The base field from which we will perform transformations
  S:=V^8 + V^4 + V^3 + V + 1;
  Fr<z>:=ext<F | S>;

  // Generate all roots of P
  pRoots:=[];
  for e in F4 do
    if Evaluate(P, e) eq 0 then
      pRoots:=Append(pRoots,e);
    end if;
  end for;

  // Embed GF(2^2) in GF(2^8)
  try 
    Embed(F4, Fr); 
  catch  e
    print("Embedding already done. No harm no foul.");
  end try;

  // Storage containers for the irreducible polynomials
  pIps:=Append([], P);
  qIps:=[];
  rIps:=[];
  sIps:=[];

  // Find all elements in GF(2^2) that make q(x) irreducible
  for pr in F4 do
    pol4<W>:=PolynomialRing(F4);
    Q:=W^2+W+pr;
    if IsIrreducible(Q) and Q notin qIps then
      F16<w>:=ext<F4 | Q>;
      qIps:=Append(qIps, Q);

      // Embed GF(2^4)/GF(2^2) in GF(2^8)
      try 
        Embed(F16, Fr); 
      catch  e
        print("Embedding already done. No harm no foul.");
      end try;

      // Find the roots of q(x) - GF(2^4)/GF(2^2)
      qRoots:=[];
      for e1 in F16 do
        if Evaluate(Q, e1) eq 0 then
          qRoots:=Append(qRoots,e1);
        end if;
      end for;

      // Find all elements in GF(2^4)/GF(2^2) that make r(y) irreducible
      for qr in F16 do
        pol8<X>:=PolynomialRing(F16);
        R:=X^2+X+qr;

		if IsIrreducible(R) and R notin rIps then
			F256<x>:=ext<F16 | R>;
			rIps:=Append(rIps, R);

			// Embed GF(2^8)/GF(2^4)/GF(2^2) in GF(2^8)
			try 
				Embed(F256, Fr); 
			catch  e
				print("Embedding already done. No harm no foul.");
			end try;

			// Find the roots of r(y) - GF(2^8)/GF(2^4)/GF(2^2)
			rRoots:=[];
			for e2 in F256 do 
				if Evaluate(R, e2) eq 0 then
					rRoots:=Append(rRoots,e2);
				end if;
			end for;

			// Find all elements in GF(2^8)/GF(2^4)/GF(2^2) that make s(z) irreducible
			for rr in F256 do
				pol16<Y>:=PolynomialRing(F256);
				S:=Y^2+Y+rr;
				if IsIrreducible(S) and S notin sIps then
					// newLambda := changeLambdaRoot(1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2]; 
					// D := Eltseq(newLambda)[1]; 
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;

					// Print out the pairs of roots and whatnot (then count the number of lines to verify correctness)
					// newLambda := changeLambdaRoot(1, pRoots[1], 1, qRoots[1], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					_<v,w,x> := PolynomialRing(F4, 3);
					// _<w> := PolynomialRing(F16, 1);
					// _<x> := PolynomialRing(F256, 1);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					y := (v^2*w^4*x^16 + v^2*w*x^16 + v^2*w^4*x + v*w*x)*(v^2*w^4*x^16 + v*w^4*x^16 + v^2*w*x^16 + v^2*w^4*x + v*w^4*x + v^2*w*x + v*w*v);
					y;
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1];

					// newLambda := changeLambdaRoot(1, pRoots[2], 1, qRoots[1], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1];


					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1];

					// newLambda := changeLambdaRoot(1, pRoots[1], 1, qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1];


					// newLambda := changeLambdaRoot(1, pRoots[2], 1, qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1];

					// newLambda := changeLambdaRoot(1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1];

					// newLambda := changeLambdaRoot(1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1];

					// // ----

					// HERE

					// newLambda := changeLambdaRoot(1, pRoots[1], 1, qRoots[1], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[2], 1, qRoots[1], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[1], 1, qRoots[2], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[2], 1, qRoots[2], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// // y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2];

					// // ----

					// CAW CAW CAW

					// newLambda := changeLambdaRoot(1, pRoots[1], 1, qRoots[1], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[1], 1, qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[1], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(1, pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(1, pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2];

					// newLambda := changeLambdaRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rr);
					// newPi := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
					// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, pr); // sigma in GF(2^2)
					// _<A,B,x> := PolynomialRing(F16, 3);
					// C := Eltseq(newLambda)[2];
					// D := Eltseq(newLambda)[1];
					// N := newPi; 
					// // y := (A^2*(C + newPi*C + D) + B^2 * C)*x + (A^2*newPi*(C + D) + B^2*D); // polynomial basis in GF(2^8) (if rr1 == 1)
					// y := (A^2*(C + D*newPi) + D*B^2*newPi)*x^(16) + (B^2*(D + C*newPi) + C*A^2*newPi)*x; // normal basis in GF(2^8) (if rr1 != 1)
					// y;
					// print P,":",Q,":",R,":",S,":",newLambda,":",newPi,":",newSigma;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2];


					// ONLY GO UP TO HERE!!!! S roots are for GF(2^16)/GF(2^8)...

					// // ----

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];

					// // ----

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];

					// // ----

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];

					// // ----

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];

					// // ----

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];

					// // ----

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];

					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
					// print P,":",Q,":",R,":",S;
					// print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];

				end if;
			end for;
		end if;

   //      if IsIrreducible(R) and R notin rIps then		
			// newQr := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, qr); // pi in GF((2^2)^2)
			// newPr := changeSigmaRoot(1, pRoots[1], F4, pr); // sigma in GF(2^2)
			// pRoots[1], qRoots[1];
			// newQr, newPr; // pr == sigma, qr == pi
			// _<A,B,x> := PolynomialRing(F4, 3);
			// C := Eltseq(newQr)[2]; // qr == pi
			// D := Eltseq(newQr)[1]; // qr == pi
			// N := newPr; // pr == sigma
			// // y := ((C*N^2 + D)*A^2 + C*B^2)*x + ((C + D)*N*A^2 + D*B^2); // for polynomial basis!
			// y := (C*A^2 + D*N*(A^2 + B^2))*x^4 + (C*N*(A^2 + B^2) + D*B^2)*x;
			// y;
   //      end if;
      end for;
    end if;
  end for;
  return [* pIps, qIps, rIps, sIps *];
end function;

polyGenDegree2_16();
