affine_8 := function(F, fullElem)
	// Build up the domain/range map for interpolation
	Input1:=[];
	Input1:=[];
	Output1:=[];
	Output2:=[];

	repeat
	  fixed := 0;
	  affine := RandomMatrix(GF(2), 8, 8);
//		affine := Matrix(GF(2), 8, 8,
//			[
//				[1,1,1,1,1,0,0,0],
//				[0,1,1,1,1,1,0,0],
//				[0,0,1,1,1,1,1,0],
//				[0,0,0,1,1,1,1,1],
//				[1,0,0,0,1,1,1,1],
//				[1,1,0,0,0,1,1,1],
//				[1,1,1,0,0,0,1,1],
//				[1,1,1,1,0,0,0,1]
//			]);

	  try
		toss := affine^(-1);

//	    "Trying new invertible affine matrix";
//	    affine;
	    
	    constant := Random(F);
	    power := -1;
	    v:=Transpose(Matrix([Reverse(Eltseq(constant))]));
	    ci := affine^(-1) * v;
	    // ci;
	    cinverse := Seqelt(Reverse([ci[1][1], ci[2][1], ci[3][1], ci[4][1], ci[5][1], ci[6][1], ci[7][1], ci[8][1]]), F);
	    // cinverse;

	    Input1:=[];
	    Input2:=[];
	    Output1:=[];
	    Output2:=[];

//	    "Got here.";

	    // Generate the mapping for interpolation
	    for e in F do
	      Input1:=Append(Input1,e);

	      // 0 has no inverse, so treat it carefully
	      if e eq 0 then
		s:=ElementToSequence(e);
	      else
		s:=ElementToSequence(e^power);
	      end if;

	      // perform the matrix computation
	      v:=Transpose(Matrix([Reverse(s)]));
	      prod:=affine*v;

	      // Transform back to the output (we don't transpose...)
	      es:=
		[
		  prod[1][1], prod[2][1], prod[3][1], prod[4][1], prod[5][1], prod[6][1], prod[7][1], prod[8][1]
		];
	      elem:=SequenceToElement(Reverse(es), F) + constant;

	      if (elem in Output1) then
//		"Duplicate element found.";
		fixed := 1;
		break;
	      end if;

	      Output1:=Append(Output1, elem);

	      if (e + elem) eq 0 or (e + elem) eq fullElem then
//		"Found a fixed point. Discarding this matrix";
		fixed := 1;
		break;
	      end if;

	      // Now do the inverse S-box
	      Input2 := Append(Input2, elem); // elem = y
	      s:=Reverse(Eltseq(elem));

	      // Perform the matrix computation
	      v:=Transpose(Matrix([s]));
	      prod:=affine^(-1) * v; // A^(-1) * y

	      // Transform back to the output (we don't transpose...)
	      es:=Reverse([prod[1][1], prod[2][1], prod[3][1], prod[4][1], prod[5][1], prod[6][1], prod[7][1], prod[8][1]]);
	      elem := SequenceToElement(es, F) + cinverse; // cinverse = A^(-1) * b
	      if elem ne 0 then
		elem := elem^power;
	      end if;
	      Output2 := Append(Output2, elem);

	      // "output";
	      // elem;

	      if Eltseq(elem) ne Eltseq(e) then
//		"not equal - exiting now.";
		fixed := 1;
		break;
	      end if;

	    end for;

	    if fixed eq 0 then
	      // Perform Lagrangian interpolation and spit out the complexity (# terms)
	      affine;
	      constant;
	      NumberOfNonZeroEntries(affine);
	      affine^-1;
	      cinverse;
	      NumberOfNonZeroEntries(affine^-1);
	      constant;
	      Fx:=Interpolation(Input1,Output1);
	      Fx;
	      //Fx1:=Interpolation(Input2,Output2);
	      //Fx1;
	      #Terms(Fx);
	      //#Terms(Fx1);

		// Return the stuff for the affine transformation
		return affine, constant, affine^(-1), cinverse, Fx, #Terms(Fx);
	    end if;

	  catch e
	  //  e;
	    fixed := 1;
	  end try;
	until fixed eq 0;
end function;

affine_16 := function(F, fullElem)
	// Build up the domain/range map for interpolation
	Input1:=[];
	Input1:=[];
	Output1:=[];
	Output2:=[];

	repeat
	  fixed := 0;
	  affine := RandomMatrix(GF(2), 16, 16);

	  try
	    toss := affine^(-1);

//	    "Trying new invertible affine matrix";
//	    affine;
	//    constant := x^16 + x^5 + x^3 + x + 1;
	    constant := Random(F);
	    power := -1;
	    v:=Transpose(Matrix([Reverse(Eltseq(constant))]));
	    ci := affine^(-1) * v;
	    // ci;
	    cinverse := Seqelt(Reverse([ci[1][1], ci[2][1], ci[3][1], ci[4][1], ci[5][1], ci[6][1], ci[7][1], ci[8][1], ci[9][1], ci[10][1], ci[11][1], ci[12][1], ci[13][1], ci[14][1], ci[15][1], ci[16][1]]), F);
	    // cinverse;

	    Input1:=[];
	    Input2:=[];
	    Output1:=[];
	    Output2:=[];

//	    "Got here.";

	    // Generate the mapping for interpolation
	    for e in F do
	      Input1:=Append(Input1,e);

	      // 0 has no inverse, so treat it carefully
	      if e eq 0 then
		s:=ElementToSequence(e);
	      else
		s:=ElementToSequence(e^power);
	      end if;

	      // perform the matrix computation
	      v:=Transpose(Matrix([Reverse(s)]));
	      prod:=affine*v;

	      // Transform back to the output (we don't transpose...)
	      es:=
		[
		  prod[1][1], prod[2][1], prod[3][1], prod[4][1], prod[5][1], prod[6][1], prod[7][1], prod[8][1],
		  prod[9][1], prod[10][1], prod[11][1], prod[12][1], prod[13][1], prod[14][1], prod[15][1], prod[16][1]
		];
	      elem:=SequenceToElement(Reverse(es), F) + constant;

	      if (elem in Output1) then
//		"Duplicate element found.";
		fixed := 1;
		break;
	      end if;

	      Output1:=Append(Output1, elem);

	      if (e + elem) eq 0 or (e + elem) eq fullElem then
//		"Found a fixed point. Discarding this matrix";
		fixed := 1;
		break;
	      end if;

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

	      // "output";
	      // elem;

	      if Eltseq(elem) ne Eltseq(e) then
		"not equal - exiting now.";
		fixed := 1;
		break;
	      end if;

	    end for;

	    if fixed eq 0 then
	      // Perform Lagrangian interpolation and spit out the complexity (# terms)
//	      "Affine and inverse affine found - performing interpolation.";
	      affine;
	      constant;
	      NumberOfNonZeroEntries(affine);
	      affine^-1;
	      cinverse;
	      NumberOfNonZeroEntries(affine^-1);
	      constant;
	      Fx:=Interpolation(Input1,Output1);
	      Fx;
	      if #Terms(Fx) lt 17 then
	      	return affine, constant, cinverse;
	     end if;
	      //Fx1:=Interpolation(Input2,Output2);
	      //Fx1;
	      // #Terms(Fx);
	      //#Terms(Fx1);

		// Return the stuff for the affine transformation
		// return affine, constant, affine^(-1), cinverse;
	    end if;

	  catch e
	//    e;
	    fixed := 1;
	  end try;
	until fixed eq 0;
end function;

interpolate_output := function(F, fullElem)
	// Build up the domain/range map for interpolation
	Input1:=[];
	Input1:=[];
	Output1:=[];
	Output2:=[];

	"Affine M";
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
	"Inverse";
	affine^(-1);

	GF2 := GF(2);
	c := SequenceToElement([GF2!1,GF2!1,GF2!1,GF2!0,GF2!1,GF2!1,GF2!0,GF2!1,GF2!1,GF2!0,GF2!1,GF2!0,GF2!0,GF2!0,GF2!1,GF2!0], F);
	"Constant";
	c;

	

	return 0;
end function;

// Build the base field stuff
F:=GF(2);
pol2<X>:=PolynomialRing(F);

// // Build GF(2^8)
// Q:=X^8+X^4+X^3+X+1;
// F256<x>:=ext<F | Q>;
// fullElem := x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1;
// toss := affine_8(F256, fullElem);

// Affines for GF(2^16)
Q:=X^16+X^5+X^3+X^2+1;
F6K<x>:=ext<F | Q>;
fullElem := x^15 + x^14 + x^13 + x^12 + x^11 + x^10 + x^9 + x^8 + x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x + 1;
toss := interpolate_output(F6K, fullElem);

