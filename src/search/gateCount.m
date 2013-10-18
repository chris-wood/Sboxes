// File: gateCount.m
// Author: Christopher A. Wood, www.christopher-wood.com
// Description: Count gates for the inverse mapping using 'algebraic' and common subexpression optimizations.

// Uncomment for standard basis printing
AssertAttribute(FldFin, "PowerPrinting", false);  
// SetQuitOnError(true);
//SetLogFile("gateCount.txt");

/// GATE COUNTS FOR ARITHMETIC IN GF(2^8)/GF(2^4)
canrightInv8 := function(P, Q, sigma, pr1, pr2, qr1, qr2, pi)
	if (qr1 eq 1) then // polynomial basis (z)
		case Eltseq(pi):
			when [0, sigma]:
				// print("p case 1");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 67;
					end if; 
				else
					return 67;
				end if;
			when [0, sigma^2]:
				// print("p case 2");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 67;
					end if; 
				else
					return 67;
				end if;
			when [sigma, sigma]:
				// print("p case 3");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 67;
					end if; 
				else
					return 67;
				end if;
			when [sigma^2, sigma^2]:
				// print("p case 4");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 67;
					end if; 
				else
					return 67;
				end if;
			when [1, sigma]:
				// print("p case 5");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 67;
					end if; 
				else
					return 67;
				end if;
			when [sigma, sigma^2]:
				// print("p case 6");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 67;
					end if; 
				else
					return 67;
				end if;
			when [sigma^2, sigma]:
				// print("p case 7");
				// sigma,  pr2,  pi,  Q;
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 68;
					else
						return 68;
					end if; 
				else
					return 67;
				end if;
			when [1, sigma^2]:
				// print("p case 8");
				if (pr1 eq 1) then
					if pr2 eq sigma then
						return 67;
					else
						return 68;
					end if; 
				else
					return 67;
				end if;
		end case;
	else
		return 66; // see Canright's paper - after all optimizations are applied he gets 66 for normal GF(2^4)/GF(2^2) basis
	end if;
end function; 

//////////////// START ROOT GENERATION CODE

pad := function(S, n)
	for i := 1 to n do
		S := Insert(S, 0, 0);
	end for;
	return S;
end function;

// CORRECT
buildRoot4Row := function(element, subfield) 

	// elem := Reverse(Eltseq(element));
	// elem := pad(elem, 4 - #elem);
	// element;
	// elem;
	// return elem;

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

// CORRECT
buildRoot8Row := function(element, subfield) 
	if element eq 1 then
		return [0,0,0,0,0,0,0,1];
	end if;

	// Drop down to a pair of elements in the subfield.
	elem := Eltseq(element, subfield);
	// elem := Reverse(Eltseq(element));
	// elem := pad(elem, 8 - #elem);
	// element;
	// elem;
	// return elem;

	// no w term, so insert two leading zeroes
	if #elem eq 1 then
		tmp := Eltseq(elem[1]); // [2] refers to upper bits, we must have the lower bits here
		tmp1 := Eltseq(tmp)[1]; // lower
		tmp2 := Eltseq(tmp)[2]; // upper
		row := [0,0,0,0,Eltseq(tmp2)[2], Eltseq(tmp2)[1], Eltseq(tmp1)[2], Eltseq(tmp1)[1]];
		return row;
	else
		p4 := Eltseq(Eltseq(elem)[2])[2];
		p3 := Eltseq(Eltseq(elem)[2])[1];
		p2 := Eltseq(Eltseq(elem)[1])[2];
		p1 := Eltseq(Eltseq(elem)[1])[1];
		tmp := 
			[
			Eltseq(p4)[2], Eltseq(p4)[1], 
			Eltseq(p3)[2], Eltseq(p3)[1], 
			Eltseq(p2)[2], Eltseq(p2)[1], 
			Eltseq(p1)[2], Eltseq(p1)[1]
			];
		return tmp;
	end if;
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

// sigma in GF(2^2)
// CORRECT
changeSigmaRoot_bin := function(pr1, pr2, F4, sigma) // sigma in GF(2^2)
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
	return Transpose(prod);
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

	M := Matrix(GF(2), [ r1, r2, r3, r4 ]);
	evector := [Eltseq(Eltseq(pi)[2])[2], Eltseq(Eltseq(pi)[2])[1], Eltseq(Eltseq(pi)[1])[2], Eltseq(Eltseq(pi)[1])[1]];
	ev := Transpose(Matrix(GF(2), [evector])); // correct
	prod := Transpose(M)^(-1) * ev;
	newPi := Seqelt([Seqelt([prod[4][1],prod[3][1]], F4), Seqelt([prod[2][1],prod[1][1]], F4)], F16);
	return newPi;
end function;

changePiRoot_bin := function(pr1, pr2, qr1, qr2, F4, F16, pi) // pi in GF(2^4)/GF(2^2)
	q1p1 := pr1 * qr1;
	q2p1 := pr1 * qr2;
	q1p2 := pr2 * qr1;
	q2p2 := pr2 * qr2;
	r1 := buildRoot4Row(q2p2, F4);
	r2 := buildRoot4Row(q2p1, F4);
	r3 := buildRoot4Row(q1p2, F4);
	r4 := buildRoot4Row(q1p1, F4);

	// [r1, r2, r3, r4];

	M := Matrix(GF(2), [ r1, r2, r3, r4 ]);
	evector := [Eltseq(Eltseq(pi)[2])[2], Eltseq(Eltseq(pi)[2])[1], Eltseq(Eltseq(pi)[1])[2], Eltseq(Eltseq(pi)[1])[1]];
	ev := Transpose(Matrix(GF(2), [evector])); // correct
	prod := Transpose(M)^(-1) * ev;
	return Transpose(prod);
end function;

// lambda in GF(2^8)/GF(2^4)/GF(2^2)
changeLambdaRoot_matrix := function(pr1, pr2, qr1, qr2, rr1, rr2, F4, F16, F256)
	p1q1r1 := F256 ! (pr1 * qr1 * rr1);
	p1q2r1 := F256 ! (pr1 * qr2 * rr1);
	p2q1r1 := F256 ! (pr2 * qr1 * rr1);
	p2q2r1 := F256 ! (pr2 * qr2 * rr1);
	p1q1r2 := F256 ! (pr1 * qr1 * rr2);
	p1q2r2 := F256 ! (pr1 * qr2 * rr2);
	p2q1r2 := F256 ! (pr2 * qr1 * rr2);
	p2q2r2 := F256 ! (pr2 * qr2 * rr2);
	r1 := buildRoot8Row(F256!p2q2r2, F16);
	r2 := buildRoot8Row(F256!p1q2r2, F16);
	r3 := buildRoot8Row(F256!p2q1r2, F16);
	r4 := buildRoot8Row(F256!p1q1r2, F16);
	r5 := buildRoot8Row(F256!p2q2r1, F16);
	r6 := buildRoot8Row(F256!p1q2r1, F16);
	r7 := buildRoot8Row(F256!p2q1r1, F16);
	r8 := buildRoot8Row(F256!p1q1r1, F16);

	// Create the inverse basis change matrix...
	M := Matrix(GF(2), [ r1, r2, r3, r4, r5, r6, r7, r8 ]);

	return Transpose(M);
end function;

// lambda in GF(2^8)/GF(2^4)/GF(2^2)
changeLambdaRoot_bin := function(pr1, pr2, qr1, qr2, rr1, rr2, F4, F16, F256, lambda)
	p1q1r1 := F256 ! (pr1 * qr1 * rr1);
	p1q2r1 := F256 ! (pr1 * qr2 * rr1);
	p2q1r1 := F256 ! (pr2 * qr1 * rr1);
	p2q2r1 := F256 ! (pr2 * qr2 * rr1);
	p1q1r2 := F256 ! (pr1 * qr1 * rr2);
	p1q2r2 := F256 ! (pr1 * qr2 * rr2);
	p2q1r2 := F256 ! (pr2 * qr1 * rr2);
	p2q2r2 := F256 ! (pr2 * qr2 * rr2);
	r1 := buildRoot8Row(F256!p2q2r2, F16);
	r2 := buildRoot8Row(F256!p1q2r2, F16);
	r3 := buildRoot8Row(F256!p2q1r2, F16);
	r4 := buildRoot8Row(F256!p1q1r2, F16);
	r5 := buildRoot8Row(F256!p2q2r1, F16);
	r6 := buildRoot8Row(F256!p1q2r1, F16);
	r7 := buildRoot8Row(F256!p2q1r1, F16);
	r8 := buildRoot8Row(F256!p1q1r1, F16);

	// Create the inverse basis change matrix...
	M := Matrix(GF(2), [ r1, r2, r3, r4, r5, r6, r7, r8 ]);

	// Create the element sequences
	upper := Eltseq(lambda)[2];
	lower := Eltseq(lambda)[1];
	p4 := Eltseq(upper)[2];
	p3 := Eltseq(upper)[1];
	p2 := Eltseq(lower)[2];
	p1 := Eltseq(lower)[1];
	evector := [Eltseq(p4)[2], Eltseq(p4)[1], Eltseq(p3)[2], Eltseq(p3)[1], Eltseq(p2)[2], Eltseq(p2)[1], Eltseq(p1)[2], Eltseq(p1)[1]];
	ev := Transpose(Matrix(GF(2), [evector])); 
	prod := Transpose(M)^(-1) * ev;

	return Transpose(prod);
end function;

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

	// Create the inverse basis change matrix...
	M := Matrix(GF(2), [ r1, r2, r3, r4, r5, r6, r7, r8 ]);

	// Create the element sequences
	upper := Eltseq(lambda)[2];
	lower := Eltseq(lambda)[1];
	p4 := Eltseq(upper)[2];
	p3 := Eltseq(upper)[1];
	p2 := Eltseq(lower)[2];
	p1 := Eltseq(lower)[1];
	evector := [Eltseq(p4)[2], Eltseq(p4)[1], Eltseq(p3)[2], Eltseq(p3)[1], Eltseq(p2)[2], Eltseq(p2)[1], Eltseq(p1)[2], Eltseq(p1)[1]];
	ev := Transpose(Matrix(GF(2), [evector])); 
	prod := Transpose(M)^(-1) * ev;

	newLambda := Seqelt( 
		[
			Seqelt( [ Seqelt([prod[8][1],prod[7][1]], F4), Seqelt([prod[6][1],prod[5][1]], F4) ], F16), 
			Seqelt( [ Seqelt([prod[4][1],prod[3][1]], F4), Seqelt([prod[2][1],prod[1][1]], F4) ], F16)
		], F256);
	return newLambda;
end function;

// TODO: need a function that takes a list of elements and makes a matrix from their binary representation
// 1. For each element in the field, perform mapping and check to see if homomorphism holds
// 2. If so, also check to see if the matrix has an inverse
// 3. If so, do the matrix multiplication

// NOTE: H(x) = Tx, T^-1 is made by inserting 2^i for all i into the matrix (which means that its independent of alpha/beta)
buildIsomorphismMatrix := function(alpha, beta, k, F, map) 
	rows := [];

	// Building the matrix rows (see isoGen.py for an explanation of how this is done)
	// alpha;
	// beta;
	for i := 1 to k do
		// Build the row
		hit := [];
		for j := 1 to (i - 1) do
			hit := Insert(hit, 1, F ! 0);
		end for;
		hit := Append(hit, F ! 1);
		for j := i + 1 to k do
			hit := Append(hit, F ! 0);
		end for;
		// hit;

		index := 0;
		for j := 1 to 2^k do
			if Eltseq(F ! map[j][2]) eq hit then
				index := j;
			end if;
		end for;

		// "checking index ", index;
		rows := Insert(rows, 1, Reverse(Eltseq(F ! map[index][3])));
	end for;
	// rows;

	Ti := Matrix(GF(2), rows);
	// Ti;
	// "";
	Ti := Transpose(Ti);
	// Ti;
	// "";
	T := Ti^(-1);
	return [T, Ti];
end function;

// Function to check and see if the homomorphism is valid
checkHomomorphism := function(alpha, beta, F1, F2, k)
	a_map := [];
	b_map := [];
	map := [];

	// We assume a mapping from alpha^i to beta^i for multiplicative homomorphism
	for i := 0 to 2^k - 1 do 
		a_map := Append(a_map, [i, alpha^i]);
		b_map := Append(b_map, [i, beta^i]);
		map := Append(map, [i, alpha^i, beta^i]);
	end for;

	// Perform the additive homomorphism check:
	// If alpha^r = alpha^i + 1, then beta^r = beta^i + 1
	for i := 0 to 2^k - 1 do
		alpha_r := alpha^i + F1!1;
		beta_r := beta^i + F2!1;
		a_r := 0; // 0 maps to itself
		b_r := 0; // 0 maps to itself
		for r := 0 to 2^k - 1 do
			if alpha_r eq a_map[r + 1][2] then
				a_r := r;
			end if;
			if beta_r eq b_map[r + 1][2] then
				b_r := r;
			end if;
		end for;
		// If this doesn't hold, the mapping isn't valid.
		// if a_r eq 0 then
		// 	alpha^i, beta^i;
		// 	a_r, b_r, i;
		// end if;
		if a_r ne b_r then
			return 0, [];
		end if;
	end for;

	return 1, map;
end function;

// Implementation of Paar's algorithm to find all suitable homomorphic mappings
// F = composite field, F2 = normal field defined by polynomial P
// It's assumed that F1 is embedded in F2... still need to figure out the details of this embedding
findCompositeIsomorphism := function(F1, F2, S, k)
	if IsPrimitive(S) then
		for alpha in F1 do
			if IsPrimitive(alpha) then
				if Evaluate(S, alpha) eq F2!0 then
					"Found one: ", alpha;
				end if;
			end if;
		end for;
	else
		// Field is not primitive, so we're gonna go about this exhaustively
		// by checking all mappings found by generators
		pass := 0;
		fail := 0;
		matrixList := []; // sequence that contains all of the matrices used in this study
		for beta in F2 do
			if beta ne F2!0 and Order(beta) eq (2^k - 1) then
				for alpha in F1 do
					if alpha ne F1!0 and Order(alpha) eq (2^k - 1) then
						res, map := checkHomomorphism(alpha, beta, F1, F2, k);
						if res eq 1 then
							matrices := buildIsomorphismMatrix(alpha, beta, k, F2, map);
							// matrices = [T, Ti]
							// matrices;
							// NumberOfNonZeroEntries(matrices[1]) + NumberOfNonZeroEntries(matrices[2]);

							// Build up the output
							matrixList := Append(matrixList, matrices);

							return matrixList;

							pass := pass + 1;
						else
							// "Found one that didn't work...";
							fail := fail + 1;
						end if;
					end if;
				end for;
			end if;
		end for;

		// Return the final set of matrices that the isomorphism checker should then verify
		return matrixList;
	end if;

	return 0; // default for now
end function;

// TODO: this is not working... really, very big sad face...
testBasisChange8 := function(M, MI, F4, F16, F256, field)

	"HERE WE GO";
	M;
	"";
	MI;
	"";
	test := Random(field);
	test;
	test^-1;
	testseq := Reverse(Eltseq(test));
	testseq := pad(testseq, 8 - #testseq);
	testseq;

	vector := Transpose(Matrix(GF(2),[testseq]));
	sv := MI * vector;
	Transpose(sv);

	subelem := Seqelt( 
	[
		Seqelt( [ Seqelt([sv[8][1],sv[7][1]], F4), Seqelt([sv[6][1],sv[5][1]], F4) ], F16), 
		Seqelt( [ Seqelt([sv[4][1],sv[3][1]], F4), Seqelt([sv[2][1],sv[1][1]], F4) ], F16)
	], F256);

	"non";
	subelem;

	// test^-1;
	testseq := Reverse(Eltseq(test^(-1)));
	testseq := pad(testseq, 8 - #testseq);
	testseq;

	vector := Transpose(Matrix(GF(2),[testseq]));
	sv := MI * vector;
	Transpose(sv);

	subeleminv := Seqelt( 
	[
		Seqelt( [ Seqelt([sv[8][1],sv[7][1]], F4), Seqelt([sv[6][1],sv[5][1]], F4) ], F16), 
		Seqelt( [ Seqelt([sv[4][1],sv[3][1]], F4), Seqelt([sv[2][1],sv[1][1]], F4) ], F16)
	], F256);

	"inv";
	subeleminv;

	// Test that inversion in standard field is the same as inversion in composite field!
	sinverse := subelem^(-1);
	sinverse;

	subelem * subeleminv;

	assert sinverse eq subeleminv;

	// buildRoot8Row(sinverse, F16);
	// prod := M * Transpose(Matrix(GF(2),[buildRoot8Row(sinverse, F16)]));
	// prod;
	// prodseq := [prod[8][1], prod[7][1], prod[6][1], prod[5][1], prod[4][1], prod[3][1], prod[2][1], prod[1][1]];
	// inverse := Seqelt(prodseq, field);

	// // prod := M * Transpose(Matrix(GF(2),[buildRoot8Row(other, F16)]));
	// // prod;
	// // prodseq := [prod[8][1], prod[7][1], prod[6][1], prod[5][1], prod[4][1], prod[3][1], prod[2][1], prod[1][1]];
	// // addee := Seqelt(prodseq, field);

	// inverse; // squared in isomorphic field... 
	// test; // squared in extension field...
	// test^(-1);

	// assert test^(-1) eq inverse;
	quit;

	return 0;
end function;

// basis change matrices using  in GF(2^8)/GF(2^4)/GF(2^2)
totalGateCount8 := function(invCount, v, w, x, M, MI, prr1, prr2, 
	qrr1, qrr2, rrr1, rrr2, F4, F16, F256, field, AFFINE, C)
	// Create the basis change in GF(((2^2)^2)^2)
	M2 := changeLambdaRoot_matrix(prr1, prr2, qrr1, qrr2, rrr1, rrr2, F4, F16, F256);
	M := M * M2;
	MI := M2^(-1) * MI;
	// M * M2;
	// M2^(-1) * MI;

	// Display the elements that are needed for the isomorphism...
	field ! v;
	field ! w;
	field ! x;

	// Display all of the basis change matrices and the important combinations
	// merged encryption sbox
	// (M)^-1; // T^(-1)
	MI;
	mInvCount := NumberOfNonZeroEntries(MI);
	mInvCount;
	prod := (AFFINE * M); // MT
	prod;
	maCount := NumberOfNonZeroEntries(prod);
	maCount;

	// merged decryption sbox
	prod := (AFFINE * M)^(-1); // (MT)^(-1)
	prod;
	maInvCount := NumberOfNonZeroEntries(prod);
	maInvCount;
	M; // T
	mCount := NumberOfNonZeroEntries(M);
	mCount;

	cCount := 2 * NumberOfNonZeroEntries(Matrix(GF(2),[Eltseq(C)]));
	Reverse(Eltseq(C));
	cCount;

	// Display the unoptimized totals 
	// Forward S-box = mInvCount + inv + maCount + cCount
	mInvCount + invCount + maCount + cCount;
	// Inverse S-box = maInvCount + inv + mCount + cCount
	maInvCount + invCount + mCount + cCount;
	// Merged = (mInvCount + maInvCount) + inv + (maCount + mCount) + cCount
	(mInvCount + maInvCount) + invCount + (maCount + mCount) + cCount;
	return 0; // dummy result
end function;

coeffMap2 := function(sigma, F, F4) // these are the only two possibilities for p(v) to be irreducible...
	// sigma;
	if Eltseq(sigma) eq [1,0] then
		return Seqelt([F ! 0,F ! 1], F4); // if sigma = 1, then we want to return v
	elif Eltseq(sigma) eq [0,1] then
		return Seqelt([F ! 1,F ! 1], F4); // if sigma = v, then we want to return v^2 (v + 1)
	end if;
end function;

// Map the normal basis coefficients to those in polynomial
// Magma doesn't use normal bases
coeffMap4 := function(sigma, pi, F4, F16)
	// pi and sigma are in a normal basis representation
	// this function exists because Magma doesn't have support for operations on GF elements in normal bases
	// The mapping is defined as follows:
	// 0                (old pi) -> 0
	// 1                (old pi) -> sigma^2 
	// sigma^2          (old pi) -> sigma
	// sigma^2 + sigma  (old pi) -> 1
	c1 := F4 ! 0; // upper coefficient (c2 w^4 + c1 w) -> [c1, c2]
	c2 := F4 ! 0; // lower coefficient (c2 w^4 + c1 w) -> [c1, c2]
	case Eltseq(sigma):
		when [0, 1]: // sigma (v)
			case Eltseq(Eltseq(pi)[1]): // c1
				when [1,0]: 
					// "c1 v - 1";
					c1 := sigma;
				when [1,1]:
					// "c1 v - v + 1";
					// c1 := 1;
					// c1 := Seqelt([1], F4);
					c1 := F4 ! 1;
				when [0,1]:
					// "c1 v - v";
					c1 := sigma^2;
			end case;
			case Eltseq(Eltseq(pi)[2]): // c2
				when [1,0]:
					// "c2 v - 1";
					c2 := sigma;
				when [1,1]:
					// "c2 v - v + 1";
					// c2 := 1;
					// c2 := Seqelt([1], F4);
					c2 := F4 ! 1;
				when [0,1]:
					// "c2 v - v";
					c2 := sigma^2;
			end case;
		when [1,1]: // sigma^2 (v^2)
			case Eltseq(Eltseq(pi)[1]): // c1
				when [1,0]:
					// "c1 1 - 1";
					c1 := sigma^2;
				when [1,1]:
					// "c1 1 - v + 1";
					// c1 := 1;
					// c1 := Seqelt([1], F4);
					c1 := F4 ! 1;
				when [0,1]:
					// "c1 1 - v";
					c1 := sigma;
			end case;
			case Eltseq(Eltseq(pi)[2]): // c2
				when [1,0]:
					// "c2 1 - 1";
					c2 := sigma^2;
				when [1,1]:
					// "c2 1 - v + 1";
					// c2 := 1;
					// c2 := Seqelt([1], F4);
					c2 := F4 ! 1;
				when [0,1]:
					// "c2 1 - v";
					c2 := sigma;
			end case;
	end case;

	// if (Eltseq(pi)[1] eq 0) then
	// 	c1 := 0;
	// end if;
	// if Eltseq(pi)[2] eq 0 then
	// 	c2 := 0;
	// end if;

	// c1;
	// c2;
	// [c1, c2];

	return Seqelt([c1, c2], F16);
end function;

// Generate all tuples of basis elements (p, r, q, s)
allGen_8 := function(field, S, AFFINE, C)
	F:=GF(2);
	pol2<V>:=PolynomialRing(F); // polynomial ring over F
	P:=V^2+V+1;                 // the only irreducible polynomial for GF(2^2)
	F4<v>:=ext<F | P>;          // create GF(2^2)/GF(2)

	// Embed GF(2^2) in GF(2^8)
	try 
		Embed(F4, field); 
	catch  e
		print("Embedding already done. No harm no foul.");
	end try;

	// Generate all roots of P
	pRoots:=[];
	for e in F4 do
		if Evaluate(P, e) eq 0 then
			pRoots:=Append(pRoots,e);
		end if;
	end for;

	// Find all elements in GF(2^2) that make q(x) irreducible
	for sigma in F4 do
		pol4<W>:=PolynomialRing(F4);
		Q:=W^2 + W + sigma;
		if IsIrreducible(Q) then
			F16<w>:=ext<F4 | Q>;

			// Embed GF(2^4)/GF(2^2) in GF(2^8)
			try 
				Embed(F16, field); 
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
			for pi in F16 do
				pol8<X>:=PolynomialRing(F16);
				R:=X^2 + X + pi;

				if IsIrreducible(R) then
					F256<x>:=ext<F16 | R>;

					// Embed GF(2^8)/GF(2^4)/GF(2^2) in GF(2^8)
					try 
						Embed(F256, field); 
					catch e
						print("Embedding already done. No harm no foul.");
					end try;

					// Find the roots of r(y) - GF(2^8)/GF(2^4)/GF(2^2)
					rRoots:=[];
					for e2 in F256 do 
						if Evaluate(R, e2) eq 0 then
							rRoots:=Append(rRoots,e2);
						end if;
					end for;

					// Find all basis change composite field mappings that work for this polynomial
					matrixList := findCompositeIsomorphism(F256, field, S, 8);
					for mi := 1 to #matrixList do
						T := matrixList[mi][1];
						Ti := matrixList[mi][2];

						// New matrix to check!
						T;
						Ti;

						toss := testBasisChange8(T, Ti, F4, F16, F256, field);

						// // Display the number of gates for this combo...
					    // P, Q, R, S, sigma, pi, 1, pRoots[1], 1, qRoots[1], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[1], 1, qRoots[1], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[1], 1, qRoots[1], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[1], 1, qRoots[1], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[1], 1, qRoots[1], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[1], 1, qRoots[1], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[1], 1, qRoots[2], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[1], 1, qRoots[2], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[1], 1, qRoots[2], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[1], 1, qRoots[2], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[1], 1, qRoots[2], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[1], 1, qRoots[2], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[2], 1, qRoots[1], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[2], 1, qRoots[1], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[2], 1, qRoots[1], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[2], 1, qRoots[1], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[2], 1, qRoots[1], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[2], 1, qRoots[1], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[2], 1, qRoots[2], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[2], 1, qRoots[2], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[2], 1, qRoots[2], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// P, Q, R, S, sigma, pi, 1, pRoots[2], 1, qRoots[2], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[2], 1, qRoots[2], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[2], 1, qRoots[2], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// // USE COEFFICIENT MAPPING
						// newPi := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						// newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// P, Q, R, S, newSigma, newPi, pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// // USE COEFFICIENT MAPPING
						// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[2], F4, F16, pi);
						// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						// newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// P, Q, R, S, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// // USE COEFFICIENT MAPPING
						// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, pi);
						// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						// newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// P, Q, R, S, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// USE COEFFICIENT MAPPING
						newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, pi);
						newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma); 
						newSigma := coeffMap2(newSigma, F, F4); // new
						newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma, so pi is represented in terms of sigma
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						pRoots[1], pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, pRoots[1], pRoots[2], 1, qRoots[1], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, T, Ti, pRoots[1], pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);

						// // USE COEFFICIENT MAPPING
						// newPi := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						// newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// P, Q, R, S, newSigma, newPi, pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// // USE COEFFICIENT MAPPING
						// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[2], F4, F16, pi);
						// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						// newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// P, Q, R, S, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], 1, qRoots[2], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// // USE COEFFICIENT MAPPING
						// newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, pi);
						// newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						// newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// P, Q, R, S, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], 1, qRoots[1], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// USE COEFFICIENT MAPPING
						newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[2], F4, F16, pi);
						newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						newSigma := coeffMap2(newSigma, F, F4); // new
						newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						pRoots[1], pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, pRoots[1], pRoots[2], 1, qRoots[2], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, T, Ti, pRoots[1], pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);

						// USE COEFFICIENT MAPPING
						newPi := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma);
						newSigma := coeffMap2(newSigma, F, F4); // new
						newPi := coeffMap4(newSigma, newPi, F4, F16); // for sigma
						// newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						// newRr := changeLambdaRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, rRoots[1]);
						// newQr := changePiRoot(pRoots[1], pRoots[2], qRoots[1], qRoots[2], F4, F16, qRoots[1]); // pi in GF((2^2)^2)
						// newPr := changeSigmaRoot(pRoots[1], pRoots[2], F4, pRoots[1]); // sigma in GF(2^2)
						// field!newPr, field!newQr, field!newRr;
						// field!pRoots[1], field!qRoots[1], field!rRoots[1];
						// field!pRoots[2], field!qRoots[2], field!rRoots[2];
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						pRoots[1], pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, newSigma, newPi, pRoots[1], pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, pRoots[1], pRoots[2], qRoots[1], qRoots[2], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, T, Ti, pRoots[1], pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);

						// newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, pi);
						// P, Q, R, S, sigma, newPi, 1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, newPi, 1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						// P, Q, R, S, sigma, newPi, 1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, newPi, 1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[2], F4, F16, F256, field, AFFINE, C);

						// newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, pi);
						// P, Q, R, S, sigma, newPi, 1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, sigma, newPi, 1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[1], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						// newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						// P, Q, R, S, sigma, newPi, 1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1];
						// invCount := gatesInv8(P, Q, R, sigma, newPi, 1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1]);
						// toss := totalGateCount8(invCount, v, w, x, 1, pRoots[2], qRoots[1], qRoots[2], 1, rRoots[1], F4, F16, F256, field, AFFINE, C);

						newPi := changePiRoot(1, pRoots[1], qRoots[1], qRoots[2], F4, F16, pi);
						newSigma := changeSigmaRoot(1, pRoots[1], F4, sigma); 
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						1, pRoots[1], qRoots[1], qRoots[2], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, newPi, 1, pRoots[1], qRoots[1], qRoots[2], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, 1, pRoots[1], qRoots[1], qRoots[2], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, 1, T, Ti, pRoots[1], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);



						newSigma := changeSigmaRoot(1, pRoots[2], F4, sigma); 
						newPi := changePiRoot(1, pRoots[2], qRoots[1], qRoots[2], F4, F16, pi);
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						1, pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, newPi, 1, pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, 1, pRoots[2], qRoots[1], qRoots[2], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, 1, T, Ti, pRoots[2], qRoots[1], qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);



						newSigma := changeSigmaRoot(1, pRoots[1], F4, sigma); 
						newPi := changePiRoot(1, pRoots[1], 1, qRoots[1], F4, F16, pi);
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						1, pRoots[1], 1, qRoots[1], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[1], 1, qRoots[1], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, 1, pRoots[1], 1, qRoots[1], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, 1, T, Ti, pRoots[1], 1, qRoots[1], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);


						newSigma := changeSigmaRoot(1, pRoots[2], F4, sigma); 
						newPi := changePiRoot(1, pRoots[2], 1, qRoots[1], F4, F16, pi);
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						1, pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, 1, pRoots[2], 1, qRoots[1], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, 1, T, Ti, pRoots[2], 1, qRoots[1], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);


						newSigma := changeSigmaRoot(1, pRoots[1], F4, sigma); 
						newPi := changePiRoot(1, pRoots[1], 1, qRoots[2], F4, F16, pi);
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						1, pRoots[1], 1, qRoots[2], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[1], 1, qRoots[2], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, 1, pRoots[1], 1, qRoots[2], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, 1, T, Ti, pRoots[1], 1, qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);

						newSigma := changeSigmaRoot(1, pRoots[2], F4, sigma); 
						newPi := changePiRoot(1, pRoots[2], 1, qRoots[2], F4, F16, pi);
						P, Q, R, S; 
						sigma; 
						newSigma;
						pi; 
						newPi;
						PrimitiveElement(field);
						1, pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2];
						// invCount := gatesInv8(P, Q, R, sigma, pi, 1, pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2]);
						invCount := canrightInv8(P, Q, newSigma, 1, pRoots[2], 1, qRoots[2], newPi);
						AFFINE; invCount; 
						toss := totalGateCount8(invCount, v, w, x, 1, T, Ti, pRoots[2], 1, qRoots[2], rRoots[1], rRoots[2], F4, F16, F256, field, AFFINE, C);
					end for;
				end if;
			end for;
		end if;
	end for;
	return 0; // success
end function;

////////////////////////////////////////////////////////////
///////////////    START AES S-BOX SEARCH    ///////////////
////////////////////////////////////////////////////////////

F:=GF(2);
polRing<V>:=PolynomialRing(F);
S := V^8 + V^4 + V^3 + V + 1; // fixed affine for AES polynomial
F256<x>:=ext<F | S>;
affine:=Matrix(GF(2),8,8, // affine == M
	[
		[1,1,1,1,1,0,0,0],
		[0,1,1,1,1,1,0,0],
		[0,0,1,1,1,1,1,0],
		[0,0,0,1,1,1,1,1],
		[1,0,0,0,1,1,1,1],
		[1,1,0,0,0,1,1,1],
		[1,1,1,0,0,0,1,1],
		[1,1,1,1,0,0,0,1]
	]);
constant := x^6 + x^5 + x + 1;
toss := allGen_8(F256, S, affine, constant); // capture output

// buildIsomorphismMatrix(8, F256);

