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
// START ROOT/POLY GENERATION
//////////////////////////////
F:=GF(2);
pol2<W>:=PolynomialRing(F); // polynomial ring over F
P:=W^2+W+1;                 // the only irreducible polynomial for GF(2^2)
F4<w>:=ext<F | P>;          // create GF(2^2)/GF(2)

// Generate all roots of P
pRoots:=[];
for e in F4 do
  if Evaluate(P, e) eq 0 then
    pRoots:=Append(pRoots,e);
  end if;
end for;
pRoots;

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

        for rr in rRoots do
          print pr,qr,rr;
          print P;
          print Q;
          print R;
        end for;

        // // Now go up one level higher and generate all irreducible polynomials
        // for rr in rRoots do
        //   pol256<Z>:=PolynomialRing(F256);
        //   S:=Z^2+Z+rr;
        //   if IsIrreducible(S) then
        //     S;
        //     F6K<z>:=ext<F256 | S>;
        //     sRoots:=[];
        //     for e3 in F6K do
        //       if Evaluate(S, e3) eq 0 then
        //         sRoots:=Append(sRoots, e3);
        //       end if;
        //     end for;

        //     print("ROOTS MOFO");
        //     sRoots;
        //     Pol6K<A>:=PolynomialRing(F6K);

        //     // Plug in the irreducible polynomial and solve for the basis change matrix
        //     Q:=A^16+A^5+A^3+A^2+1;
        //     print("ROOTS FOR A^16+A^5+A^3+A^2+1;");
        //     Roots(Q);
        //     for beta1 in Roots(Q) do
        //       beta:=beta1[1];

        //       print "";
        //       print "";

        //       result:=Transpose(Matrix([
        //         [
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[2],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[1],
        //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[2]
        //         ] : i in [0..15]]));
        //       result;

        //       // for i in [0..15] do
        //       //   beta^i;
        //       // end for;
        //     end for;
        //   end if;

        //   // Try the powers of this root...
        //   rRootsPowers:=rootPowers(rr, 65536, rRoots);
        //   for rr2 in rRootsPowers do
        //     pol256<Z>:=PolynomialRing(F256);
        //     S:=Z^2+Z+rr2;
        //     if IsIrreducible(S) then
        //       S;
        //       F6K<z>:=ext<F256 | S>;
        //       sRoots:=[];
        //       for e3 in F6K do
        //         if Evaluate(S, e3) eq 0 then
        //           sRoots:=Append(sRoots, e3);
        //         end if;
        //       end for;

        //       print("ROOTS MOFO");
        //       sRoots;
        //       Pol6K<A>:=PolynomialRing(F6K);

        //       // Plug in the irreducible polynomial and solve for the basis change matrix
        //       Q:=A^16+A^5+A^3+A^2+1;
        //       print("ROOTS FOR A^16+A^5+A^3+A^2+1;");
        //       Roots(Q);
        //       for beta1 in Roots(Q) do
        //         beta:=beta1[1];

        //         print "";
        //         print "";

        //         result:=Transpose(Matrix([
        //           [
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[2],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[1],
        //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[2]
        //           ] : i in [0..15]]));
        //         result;

        //         // for i in [0..15] do
        //         //   beta^i;
        //         // end for;
        //       end for;
        //     end if;
        //   end for;
        // end for;
      end if;

      // Try the other powers of q, as they might make the polynomial irreducible
      qRootsPowers:=rootPowers(qr, 256, qRoots);
      polys:=[];
      for qr2 in qRootsPowers do
        pol16<Y>:=PolynomialRing(F16);
        R:=Y^2+Y+qr2;
        if IsIrreducible(R) then

          // Add R to the list of irreducible polynomials we've seen up to this degree
          // weakAppend([], R);
          // pr;
          // qr2;
          // R;

          F256<y>:=ext<F16 | R>;
          rRoots:=[];
          for e2 in F256 do
            if Evaluate(R, e2) eq 0 then
              rRoots:=Append(rRoots, e2);
            end if;
          end for;

          for rr in rRoots do
            print pr,qr,rr;
            print P;
            print Q;
            print R;
          end for;

          // // Now go up one level higher and generate all irreducible polynomials
          // for rr in rRoots do
          //   pol256<Z>:=PolynomialRing(F256);
          //   S:=Z^2+Z+rr;
          //   if IsIrreducible(S) then
          //     S;
          //     F6K<z>:=ext<F256 | S>;
          //     sRoots:=[];
          //     for e3 in F6K do
          //       if Evaluate(S, e3) eq 0 then
          //         sRoots:=Append(sRoots, e3);
          //       end if;
          //     end for;

          //     print("ROOTS MOFO");
          //     sRoots;
          //     Pol6K<A>:=PolynomialRing(F6K);

          //     // Plug in the irreducible polynomial and solve for the basis change matrix
          //     Q:=A^16+A^5+A^3+A^2+1;
          //     print("ROOTS FOR A^16+A^5+A^3+A^2+1;");
          //     Roots(Q);
          //     for beta1 in Roots(Q) do
          //       beta:=beta1[1];

          //       print "";
          //       print "";

          //       result:=Transpose(Matrix([
          //         [
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[2],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[1],
          //         Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[2]
          //         ] : i in [0..15]]));
          //       result;

          //       // for i in [0..15] do
          //       //   beta^i;
          //       // end for;
          //     end for;
          //   end if;

          //   // Try the powers of this root...
          //   rRootsPowers:=rootPowers(rr, 65536, rRoots);
          //   for rr2 in rRootsPowers do
          //     pol256<Z>:=PolynomialRing(F256);
          //     S:=Z^2+Z+rr2;
          //     if IsIrreducible(S) then
          //       S;
          //       F6K<z>:=ext<F256 | S>;
          //       sRoots:=[];
          //       for e3 in F6K do
          //         if Evaluate(S, e3) eq 0 then
          //           sRoots:=Append(sRoots, e3);
          //         end if;
          //       end for;

          //       print("ROOTS MOFO HERE");
          //       sRoots;
          //       Pol6K<A>:=PolynomialRing(F6K);

          //       // Plug in the irreducible polynomial and solve for the basis change matrix
          //       Q:=A^16+A^5+A^3+A^2+1;
          //       print("ROOTS FOR A^16+A^5+A^3+A^2+1;");
          //       Roots(Q);
          //       for beta1 in Roots(Q) do
          //         beta:=beta1[1];

          //         print "";
          //         print "";

          //         result:=Transpose(Matrix([
          //           [
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[1])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[1])[2])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[1])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[1])[2])[2])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[1])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[1])[2])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[1])[2],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[1],
          //           Eltseq(Eltseq(Eltseq(Eltseq(beta^i)[2])[2])[2])[2]
          //           ] : i in [0..15]]));
          //         result;

          //         // for i in [0..15] do
          //         //   beta^i;
          //         // end for;
          //       end for;
          //     end if;
          //   end for;
          // end for;
        end if;
      end for;
    end for;
  end if;

  // Try other powers, just in case
  pRootPowers:=rootPowers(pr, 2, pRoots);
  for pr2 in pRootPowers do
    Q:=X^2+X+pr2; // use the other roots
    if IsIrreducible(Q) then
      F16<x>:=ext<F4 | Q>;
      qRoots:=[];
      for e in F16 do
        if Evaluate(Q, e) eq 0 then
          qRoots:=Append(qRoots, e);
        end if;
      end for;
      for qr in qRoots do
        pol16<Y>:=PolynomialRing(F16);
        R:=Y^2+Y+qr; 
        if IsIrreducible(R) then
          F256<z>:=ext<F16 | R>;
          rRoots:=[];
          for e2 in F256 do
            if Evaluate(R, e2) eq 0 then
              rRoots:=Append(rRoots,e2);
            end if;
          end for;

          for rr in rRoots do
            print pr,qr,rr;
            print P;
            print Q;
            print R;
          end for;
        end if;


        // Try the other powers of q, as they might make the polynomial irreducible
        qRootsPowers:=rootPowers(qr, 256, qRoots);
        polys:=[];
        for qr2 in qRootsPowers do
          pol16<Y>:=PolynomialRing(F16);
          R:=Y^2+Y+qr2;
          if IsIrreducible(R) then

            // Add R to the list of irreducible polynomials we've seen up to this degree
            // weakAppend([], R);
            // pr;
            // qr2;
            // R;

            F256<y>:=ext<F16 | R>;
            rRoots:=[];
            for e2 in F256 do
              if Evaluate(R, e2) eq 0 then
                rRoots:=Append(rRoots, e2);
              end if;
            end for;
            
            for rr in rRoots do
              print pr,qr,rr;
              print P;
              print Q;
              print R;
            end for;

          end if;
        end for;
      end for;
    end if;
  end for;
end for;


