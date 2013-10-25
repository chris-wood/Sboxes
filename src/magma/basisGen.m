// File: basisGen.m
// Author: Christopher Wood
// Description: A set of helpful functions to generate basis elements for fields of interest.

// Uncomment for standard basis printing
// AssertAttribute(FldFin, "PowerPrinting", false);  
SetQuitOnError(true);
SetLogFile("basisGenOut.txt");

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
              F6K<y>:=ext<F256 | S>;
              sIps:=Append(sIps, S);

              // Embed GF(2^16)/GF(2^8)/GF(2^4)/GF(2^2)
              try 
                Embed(F256, Fr); 
              catch  e
                print("Embedding already done. No harm no foul.");
              end try;

              // Find the roots of s(z) - GF(2^16)/GF(2^8)/GF(2^4)/GF(2^2)
              sRoots:=[];
              for e3 in F6K do 
                if Evaluate(S, e3) eq 0 then
                  sRoots:=Append(sRoots,e3);
                end if;
              end for;

              // TODO: BASIS CHANGE MATRIX WOULD GO HERE
              P,":",Q,":",R,":",S;
              // pRoots;
              // Q;
              // qRoots;
              // R;
              // rRoots;
              // S;
              // sRoots;
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;
  return [* pIps, qIps, rIps, sIps *];
end function;

// Generate all tuples of basis elements (p, r, q, s)
polyGenDegree2_8 := function()
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

          P,":",Q,":",R;
          // pRoots;
          // Q;
          // qRoots;
          // R;
          // rRoots;
        end if;
      end for;
    end if;
  end for;
  return [* pIps, qIps, rIps *];
end function;

// Generate all tuples of basis elements (p, r, q, s)
allGen_8 := function(embed, field)
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
      if embed eq 0 then
        pRoots:=Append(pRoots,e);
      else
        pRoots:=Append(pRoots,field!e);
      end if;
    end if;
  end for;

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
      // qIps:=Append(qIps, Q);

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
          if embed eq 0 then
            qRoots:=Append(qRoots,e1);
          else
            qRoots:=Append(qRoots,field!e1);
          end if;
        end if;
      end for;

      // Find all elements in GF(2^4)/GF(2^2) that make r(y) irreducible
      for qr in F16 do
        pol8<X>:=PolynomialRing(F16);
        R:=X^2+X+qr;

        if IsIrreducible(R) and R notin rIps then
          F256<x>:=ext<F16 | R>;
          // rIps:=Append(rIps, R);

          // Embed GF(2^8)/GF(2^4)/GF(2^2) in GF(2^8)
          try 
            Embed(F256, field); 
          catch  e
            print("Embedding already done. No harm no foul.");
          end try;

          // Find the roots of r(y) - GF(2^8)/GF(2^4)/GF(2^2)
          rRoots:=[];
          for e2 in F256 do 
            if Evaluate(R, e2) eq 0 then
              if embed eq 0 then
                rRoots:=Append(rRoots,e2);
              else
                rRoots:=Append(rRoots,field!e2);
              end if;
            end if;
          end for;

          // Print out the pairs of roots and whatnot (then count the number of lines to verify correctness)
          print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1];
          print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1];
          print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1];
          print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1];
          print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2];
          print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2];
          print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2];
          print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2];

          // ----

          print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1];
          print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1];
          print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1];
          print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2];
          print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2];
          print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2];
          print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2];
          print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2];
          print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2];

          print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2];
          print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2];
          print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1];
          print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1];
          print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2];
          print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2];
          
          print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2];
          print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2];
          print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2];
          print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2];
          
        end if;
      end for;
    end if;
  end for;
  return [* pIps, qIps, rIps *];
end function;

// Generate all tuples of basis elements (p, r, q, s)
allGen_16 := function(embed, field)
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
      if embed eq 0 then
        pRoots:=Append(pRoots,e);
      else
        pRoots:=Append(pRoots,field!e);
      end if;
    end if;
  end for;

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
      // qIps:=Append(qIps, Q);

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
          if embed eq 0 then
            qRoots:=Append(qRoots,e1);
          else
            qRoots:=Append(qRoots,field!e1);
          end if;
        end if;
      end for;

      // Find all elements in GF(2^4)/GF(2^2) that make r(y) irreducible
      for qr in F16 do
        pol8<X>:=PolynomialRing(F16);
        R:=X^2+X+qr;

        if IsIrreducible(R) and R notin rIps then
          F256<x>:=ext<F16 | R>;
          // rIps:=Append(rIps, R);

          // Embed GF(2^8)/GF(2^4)/GF(2^2) in GF(2^8)
          try 
            Embed(F256, field); 
          catch  e
            print("Embedding already done. No harm no foul.");
          end try;

          // Find the roots of r(y) - GF(2^8)/GF(2^4)/GF(2^2)
          rRoots:=[];
          for e2 in F256 do 
            if Evaluate(R, e2) eq 0 then
              if embed eq 0 then
                rRoots:=Append(rRoots,e2);
              else
                rRoots:=Append(rRoots,field!e2);
              end if;
            end if;
          end for;

          // Find all elements in GF(2^8)/GF(2^4)/GF(2^2) that make s(z) irreducible
          for rr in F256 do
            pol16<Y>:=PolynomialRing(F256);
            S:=Y^2+Y+rr;

            if IsIrreducible(S) and S notin sIps then
              F6K<y>:=ext<F256 | S>;
              // sIps:=Append(sIps, S);

              // Embed GF(2^16)/GF(2^8)/GF(2^4)/GF(2^2)
              try 
                Embed(F6K, field); 
              catch  e
                print("Embedding already done. No harm no foul.");
              end try;

              // Find the roots of s(z) - GF(2^16)/GF(2^8)/GF(2^4)/GF(2^2)
              sRoots:=[];
              for e3 in F6K do 
                if Evaluate(S, e3) eq 0 then
                  sRoots:=Append(sRoots,field!e3);
                end if;
              end for;

              // Print out the pairs of roots and whatnot (then count the number of lines to verify correctness)
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[1];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[1];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[1];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[1];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[1];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[1];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[1];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":","1",":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":","1",":",sRoots[2];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":","1",":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":","1",":",sRoots[2];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":","1",":",sRoots[2];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[1],":",sRoots[1],":",sRoots[2];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":","1",":",rRoots[2],":",sRoots[1],":",sRoots[2];

              // ----

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[1],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":","1",":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];

              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[1],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print "1",":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
              print P,":",Q,":",R,":",S;
              print pRoots[1],":",pRoots[2],":",qRoots[1],":",qRoots[2],":",rRoots[1],":",rRoots[2],":",sRoots[1],":",sRoots[2];
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;
  return [* pIps, qIps, rIps, sIps *];
end function;

// Generate all tuples of basis elements (p, r, q, s)
gen16coeffs := function(embed, field)
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
      if embed eq 0 then
        pRoots:=Append(pRoots,e);
      else
        pRoots:=Append(pRoots,field!e);
      end if;
    end if;
  end for;

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
      // qIps:=Append(qIps, Q);

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
          if embed eq 0 then
            qRoots:=Append(qRoots,e1);
          else
            qRoots:=Append(qRoots,field!e1);
          end if;
        end if;
      end for;

      // Find all elements in GF(2^4)/GF(2^2) that make r(y) irreducible
      for qr in F16 do
        pol8<X>:=PolynomialRing(F16);
        R:=X^2+X+qr;

        if IsIrreducible(R) and R notin rIps then
          F256<x>:=ext<F16 | R>;
          // rIps:=Append(rIps, R);

          // Embed GF(2^8)/GF(2^4)/GF(2^2) in GF(2^8)
          try 
            Embed(F256, field); 
          catch  e
            print("Embedding already done. No harm no foul.");
          end try;

          // Find the roots of r(y) - GF(2^8)/GF(2^4)/GF(2^2)
          rRoots:=[];
          for e2 in F256 do 
            if Evaluate(R, e2) eq 0 then
              if embed eq 0 then
                rRoots:=Append(rRoots,e2);
              else
                rRoots:=Append(rRoots,field!e2);
              end if;
            end if;
          end for;

          // Find all elements in GF(2^8)/GF(2^4)/GF(2^2) that make s(z) irreducible
          for rr in F256 do
            pol16<Y>:=PolynomialRing(F256);
            S:=Y^2+Y+rr;

            if IsIrreducible(S) and S notin sIps then
              F6K<y>:=ext<F256 | S>;
              // sIps:=Append(sIps, S);
              Eltseq(rr);
            end if;
          end for;
        end if;
      end for;
    end if;
  end for;
  return [* pIps, qIps, rIps, sIps *];
end function;

////////////////////
//// START HERE ////
////////////////////

AssertAttribute(FldFin, "PowerPrinting", false);
F:=GF(2);
polRing<V>:=PolynomialRing(F);
S:=V^8+V^4+V^3+V+1;
F256<z>:=ext<F|S>;
Q:=V^16+V^5+V^3+V^2+1;
F6K<z> := ext<F|Q>;
// p:=allGen_8(1, F256);
gen16coeffs(0, F6K);

