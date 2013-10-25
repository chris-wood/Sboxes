AssertAttribute(FldFin, "PowerPrinting", false);

interpolate_output := function(F, power, fullElem)
	// Build up the domain/range map for interpolation
	Input1:=[];
	Input2:=[];
	Output1:=[];
	Output2:=[];

	"*** Affine M";
	affine := Matrix(GF(2), 16, 16,
	[
		[0,0,1,0,0,0,0,1,0,0,1,1,1,1,1,0 ],
		[1,1,0,0,0,0,0,1,0,1,1,0,1,0,1,0 ],
		[1,1,0,0,1,0,1,1,0,1,0,1,0,0,1,1 ],
		[1,1,1,0,0,0,1,0,0,1,1,0,0,0,0,0 ],
		[1,1,0,0,0,1,1,0,0,1,1,1,1,0,1,1 ],
		[0,1,0,0,0,0,1,1,0,1,1,1,1,1,0,1 ],
		[0,0,1,0,1,0,1,0,1,1,0,0,1,1,0,0 ],
		[1,0,1,1,1,0,1,1,0,0,0,1,0,1,1,1 ],
		[0,1,0,0,0,0,0,0,1,0,0,1,1,1,0,1 ],
		[1,0,1,1,0,0,0,1,0,0,1,0,1,0,0,0 ],
		[1,0,1,0,0,1,1,1,0,0,1,1,0,1,0,0 ],
		[1,0,1,1,1,0,1,1,1,1,0,1,1,0,0,1 ],
		[1,0,1,0,0,1,0,1,1,0,0,1,0,0,0,1 ],
		[0,1,0,0,0,1,1,1,1,0,0,0,0,0,0,1 ],
		[1,0,0,0,1,1,0,1,0,1,1,1,1,0,0,0 ],
		[1,1,0,1,0,1,1,0,1,0,0,1,1,0,0,0]
	]);
	affine;
	"*** Affine Inverse";
	affine^(-1);

	GF2 := GF(2);
	constant := SequenceToElement([GF2!1,GF2!1,GF2!1,GF2!0,GF2!1,GF2!1,GF2!0,GF2!1,GF2!1,GF2!0,GF2!1,GF2!0,GF2!0,GF2!0,GF2!1,GF2!0], F);
	"*** Constant";
	constant;

	// Compute the inverse (make sure it actually matches what we expect..)
	v := Transpose(Matrix([Reverse(Eltseq(constant))]));
    ci := affine^(-1) * v;
    cinverse := Seqelt(Reverse([ci[1][1], ci[2][1], ci[3][1], ci[4][1], ci[5][1], ci[6][1], ci[7][1], ci[8][1], ci[9][1], ci[10][1], ci[11][1], ci[12][1], ci[13][1], ci[14][1], ci[15][1], ci[16][1]]), F);
    "*** Constant Inverse";
    cinverse;

	for e in F do
		Input1:=Append(Input1,e);

		// Special GF inverse
		if e eq 0 then
			s:=ElementToSequence(e);
		else
			s:=ElementToSequence(e^power);
		end if;
		//e;
		//s;

		// perform the matrix computation
		v:=Transpose(Matrix([Reverse(s)]));
		prod:=affine*v;

		// Transform back to the output (we don't transpose...)
		es:=
		[
			prod[1][1], prod[2][1], prod[3][1], prod[4][1], prod[5][1], prod[6][1], prod[7][1], prod[8][1],
			prod[9][1], prod[10][1], prod[11][1], prod[12][1], prod[13][1], prod[14][1], prod[15][1], prod[16][1]
		];
		// prod;
		// es;
		elem:=SequenceToElement(Reverse(es), F);
		// elem;
		// constant;
		// quit;
		elem := elem + constant;
		//elem;

		if (elem in Output1) then
			fixed := 1;
			"*** fixed found";
			"***  --- THIS SHOULD NOT HAPPEN --- ";
			quit;
			break;
		end if;

		Output1:=Append(Output1, elem);

		// if (e + elem) eq F!0 or (e + elem) eq fullElem then
		// 	e;
		// 	elem;
		// 	e + elem;
		// 	fixed := 1;
		// 	"inverse fixed found?";
		// 	" --- THIS SHOULD NOT HAPPEN --- ";
		// 	quit;
		// 	break;
		// end if;

		// Display the element and its output
		e;
		elem;

		// Now do the inverse S-box
		Input2 := Append(Input2, elem); // elem = y
		s:=Reverse(Eltseq(elem));

		// Perform the matrix computation
		v:=Transpose(Matrix([s]));
		prod:=affine^(-1) * v; // A^(-1) * y

		// Transform back to the output (we don't transpose...)
		es:=Reverse([prod[1][1], prod[2][1], prod[3][1], prod[4][1], prod[5][1], prod[6][1], prod[7][1], prod[8][1], prod[9][1], prod[10][1], prod[11][1], prod[12][1], prod[13][1], prod[14][1], prod[15][1], prod[16][1]]);
		elem := SequenceToElement(es, F) + cinverse; // cinverse = A^(-1) * b
		if elem ne 0 then
			elem := elem^power;
		end if;

		Output2 := Append(Output2, elem);

		if Eltseq(elem) ne Eltseq(e) then
			"*** element did not match it's inverse...";
			"***  --- THIS SHOULD NOT HAPPEN --- ";
			fixed := 1;
			quit;
			break;
		end if;

	end for;

	// This should always pass...
	//if fixed eq 0 then
	"*** Performing interpolation...";
	affine;
	constant;
	NumberOfNonZeroEntries(affine);
	affine^-1;
	cinverse;
	NumberOfNonZeroEntries(affine^-1);
	constant;
	Fx:=Interpolation(Input1,Output1);
	Fx;
	return affine, constant, cinverse;
	//end if;

	return 0;
end function;

// Build up the field
F:=GF(2);
pol2<X>:=PolynomialRing(F);
Q:=X^16+X^5+X^3+X^2+1;
F6K<x>:=ext<F | Q>;
fullElem := x^15 + x^14 + x^13 + x^12 + x^11 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1;
toss := interpolate_output(F6K, -1, fullElem);
