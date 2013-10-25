// File: gateCount.m
// Author: Christopher Wood
// Description: Count gates for the inverse mapping using 'algebraic' 
// and some common subexpression optimizations.

// Uncomment for standard basis printing
AssertAttribute(FldFin, "PowerPrinting", false);  
SetQuitOnError(true);

////////////////////////////////////////////////////////////////////////
/// GATE COUNTS FOR POLYNOMIAL/NORMAL ARITHMETIC IN GF(2^2)
////////////////////////////////////////////////////////////////////////

gatesPolyInv2 := function()
    return 1; // 1 XOR
end function;

gatesPolyAdd2 := function()
    return 2; // 2 XORs
end function;

gatesPolyMult2 := function()
    // return 7; // 4 XORs, 3 ANDs
    return 4; // ONLY COUNT XOR GATES HERE
end function;

gatesPolySquare2 := function()
    return gatesPolyInv2(); // square is the same as inverse
end function;

gatesPolyScale2 := function()
    return 1; // 1 XOR
end function;

gatesPolySquareScale2 := function(pr, sigma, scalar) 
    if pr eq sigma and pr eq scalar then
        return 0;
    elif pr eq sigma^2 and pr eq scalar then
        return 0;
    else
        return 1;
    end if; 
end function;

gatesNormInv2 := function()
    return 0; // 1 XOR
end function;

gatesNormAdd2 := function()
    return 2; // 2 XORs
end function;

gatesNormMult2 := function()
    // return 7; // 4 XORs, 3 ANDs
    return 4; // only count XOR gates HERE
end function;

// DONE
gatesNormSquare2 := function()
    return gatesNormInv2(); // square is the same as inverse
end function;

// DONE
gatesNormScale2 := function()
    return 1; // 1 XOR
end function;

// DONE
gatesNormSquareScale2 := function() 
    return 1;
end function;

////////////////////////////////////////////////////////////////////////
/// GATE COUNTS FOR POLYNOMIAL ARITHMETIC IN GF(2^4)/GF(2^2)
////////////////////////////////////////////////////////////////////////

gatesInv4 := function(P, Q, sigma, pr1, pr2, qr1, qr2)
    if (qr1 eq 1) then
        if (pr1 eq 1) then
            sum := 2 * gatesPolyAdd2(); 
            sum := sum + gatesPolySquareScale2(pr2, sigma, sigma);
            sum := sum + (3 * gatesPolyMult2());
            sum := sum + gatesPolyInv2();
            return sum;
        else
            sum := 2 * gatesNormAdd2(); 
            sum := sum + gatesNormSquareScale2();
            sum := sum + (3 * gatesNormMult2());
            sum := sum + gatesNormInv2();
            return sum;
        end if; 
    else
        if (pr1 eq 1) then
            sum := 2 * gatesPolyAdd2();
            sum := sum + gatesPolySquareScale2(pr2, sigma, sigma);
            sum := sum + (3 * gatesPolyMult2());
            sum := sum + gatesPolyInv2();
            return sum;
        else
            sum := 2 * gatesNormAdd2();
            sum := sum + gatesNormSquareScale2();
            sum := sum + (3 * gatesNormMult2());
            sum := sum + gatesNormInv2();
            return sum;
        end if; 
    end if;
end function;

gatesAdd4 := function()
    return 4; 
end function;

gatesMult4 := function(P, Q, sigma, pr1, pr2, qr1, qr2)
    if (qr1 eq 1) then
        if (pr1 eq 1) then
            sum := (4 * gatesPolyAdd2());
            sum := sum + (3 * gatesPolyMult2());
            return sum;
        else
            sum := (4 * gatesNormAdd2());
            sum := sum + (3 * gatesNormMult2());
            return sum;
        end if; 
    else
        if (pr1 eq 1) then
            sum := (4 * gatesPolyAdd2());
            sum := sum + (3 * gatesPolyMult2());
            sum := sum + (gatesPolyScale2());
            return sum;
        else
            sum := (4 * gatesNormAdd2());
            sum := sum + (3 * gatesNormMult2());
            sum := sum + (gatesNormScale2());
            return sum;
        end if; 
    end if;
end function;

gatesSquare4 := function(P, Q, sigma, pr1, pr2, qr1, qr2)
    if (qr1 eq 1) then
        if (pr1 eq 1) then
            sum := (2 * gatesPolySquare2()); 
            sum := sum + gatesPolyAdd2();
            sum := sum + gatesPolyScale2();
            return sum;
        else
            sum := (2 * gatesNormSquare2()); 
            sum := sum + gatesNormAdd2();
            sum := sum + gatesNormScale2();
            return sum;
        end if; 
    else
        if (pr1 eq 1) then
            sum := (3 * gatesPolyAdd2());
            sum := sum + (2 * gatesPolySquare2());
            sum := sum + gatesPolyScale2();
            return sum;
        else
            sum := (3 * gatesNormAdd2());
            sum := sum + (2 * gatesNormSquare2());
            sum := sum + gatesNormScale2();
            return sum;
        end if; 
    end if;
end function;

gatesScale4 := function(P, Q, sigma, pr1, pr2, qr1, qr2, pi) 
    if (qr1 eq 1) then // polynomial basis (z)
        case Eltseq(pi):
            when [0, sigma]:
                return 4;
            when [0, sigma^2]:
                return 3;
            when [sigma, sigma]:
                return 4;
            when [sigma^2, sigma^2]:
                return 3;
            when [1, sigma]:
                return 6;
            when [sigma, sigma^2]:
                return 5;
            when [sigma^2, sigma]:
                return 6;
            when [1, sigma^2]:
                return 5;
        end case;
    else // Normal basis 
        case Eltseq(pi):
            when [0, sigma]:
                return 6;
            when [sigma, 0]:
                return 6;
            when [0, sigma^2]:
                return 5;
            when [sigma^2, 0]:
                return 5;
            when [1, sigma]:
                return 5;
            when [sigma, 1]:
                return 5;
            when [1, sigma^2]:
                return 6;
            when [sigma^2, 1]:
                return 6;
        end case;
    end if;
end function;

gatesSquareScale4 := function(P, Q, sigma, pr1, pr2, qr1, qr2, pi) 
    if (qr1 eq 1) then // polynomial basis (z)
        case Eltseq(pi):
            when [0, sigma]:
                if (pr1 eq 1) then
                    return (1 * gatesPolySquare2()) + (1 * gatesPolyAdd2()) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma)) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma^2));
                else
                    return (1 * gatesNormSquare2()) + (1 * gatesNormAdd2()) + 
                        (1 * gatesNormSquareScale2()) + (1 * gatesNormSquareScale2());
                end if;
            when [0, sigma^2]:
                if (pr1 eq 1) then
                    return (1 * gatesPolyAdd2()) + (1 * gatesPolySquare2()) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma)) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma^2));
                else
                    return (1 * gatesNormAdd2()) + (1 * gatesNormSquare2()) + 
                        (1 * gatesNormSquareScale2()) + (1 * gatesNormSquareScale2());
                end if;
            when [sigma, sigma]:
                if (pr1 eq 1) then
                    return (1 * gatesPolySquareScale2(pr2, sigma, sigma)) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma^2)) + (1 * gatesPolyAdd2());
                else
                    return (1 * gatesNormSquareScale2()) + 
                        (1 * gatesNormSquareScale2()) + (1 * gatesNormAdd2());
                end if;
            when [sigma^2, sigma^2]:
                if (pr1 eq 1) then
                    return (1 * gatesPolySquare2()) + 
                        (1 * gatesPolyAdd2()) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma^2));
                else
                    return (1 * gatesNormSquare2()) + 
                        (1 * gatesNormAdd2()) + (1 * gatesNormSquareScale2());
                end if;
            when [1, sigma]:
                if (pr1 eq 1) then
                    return (1 * gatesPolySquareScale2(pr2, sigma, sigma)) + 
                        (1 * gatesPolyAdd2()) + (1 * gatesPolySquare2());
                else
                    return (1 * gatesNormSquareScale2()) + 
                        (1 * gatesNormAdd2()) + (1 * gatesNormSquare2());
                end if;
            when [sigma, sigma^2]:
                if (pr1 eq 1) then
                    return gatesPolyAdd2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2) + 
                        gatesPolySquareScale2(pr2, sigma, sigma);
                else
                    return gatesNormAdd2() + gatesNormSquareScale2() + 
                        gatesNormSquareScale2();
                end if;
            when [sigma^2, sigma]:
                if (pr1 eq 1) then
                    return (2 * gatesPolyAdd2()) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma)) + 
                        gatesPolySquare2();
                else
                    return (2 * gatesNormAdd2()) + 
                        (1 * gatesNormSquareScale2()) + gatesNormSquare2();
                end if;
            when [1, sigma^2]:
                if (pr1 eq 1) then
                    return (2 * gatesPolyAdd2()) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma)) + 
                        (1 * gatesPolySquareScale2(pr2, sigma, sigma^2));
                else
                    return (2 * gatesNormAdd2()) + 
                        (1 * gatesNormSquareScale2()) + (1 * gatesNormSquareScale2());
                end if;
        end case;
    else // Normal basis
        case Eltseq(pi):
            when [0, sigma]:
                if (pr1 eq 1) then
                    return gatesPolyAdd2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma) + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2);
                else
                    return gatesNormAdd2() + gatesNormSquareScale2() + 
                        gatesNormSquareScale2();
                end if;
            when [sigma, 0]:
                if (pr1 eq 1) then
                    return gatesPolyAdd2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma) + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2);
                else
                    return gatesNormAdd2() + gatesNormSquareScale2() + 
                        gatesNormSquareScale2();
                end if;
            when [0, sigma^2]:
                if (pr1 eq 1) then
                    return gatesPolyAdd2() + gatesPolySquare2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2);
                else
                    return gatesNormAdd2() + gatesNormSquare2() + 
                        gatesNormSquareScale2();
                end if;
            when [sigma^2, 0]:
                if (pr1 eq 1) then
                    return gatesPolyAdd2() + gatesPolySquare2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2);
                else
                    return gatesNormAdd2() + gatesNormSquare2() + 
                        gatesNormSquareScale2();
                end if;
            when [1, sigma]:
                if (pr1 eq 1) then
                    return gatesPolySquareScale2(pr2, sigma, sigma) + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2) + 
                        gatesPolyAdd2(); 
                else
                    return gatesNormSquareScale2() + 
                        gatesNormSquareScale2() + gatesPolyAdd2(); 
                end if;
            when [sigma, 1]:
                if (pr1 eq 1) then
                    return gatesPolySquareScale2(pr2, sigma, sigma) + 
                        gatesPolySquareScale2(pr2, sigma, sigma^2) + 
                        gatesPolyAdd2(); 
                else
                    return gatesNormSquareScale2() + 
                        gatesNormSquareScale2() + gatesPolyAdd2(); 
                end if;
            when [1, sigma^2]:
                if (pr1 eq 1) then
                    return gatesPolySquare2() + gatesPolyAdd2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma);
                else
                    return gatesNormSquare2() + gatesNormAdd2() + 
                        gatesNormSquareScale2();
                end if;
            when [sigma^2, 1]:
                if (pr1 eq 1) then
                    return gatesPolySquare2() + gatesPolyAdd2() + 
                        gatesPolySquareScale2(pr2, sigma, sigma);
                else
                    return gatesNormSquare2() + gatesNormAdd2() + 
                        gatesNormSquareScale2();
                end if;
        end case;
    end if;

    print("ERROR - gatesSquareScale4"); // We will not get here.
    quit;
end function;

/// GATE COUNTS FOR ARITHMETIC IN GF(2^8)/GF(2^4)

canrightInv8 := function(P, Q, sigma, pr1, pr2, qr1, qr2, pi)
    if (qr1 eq 1) then // polynomial basis (z)
        case Eltseq(pi):
            when [0, sigma]:
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
        return 66; // see Canright's paper for details
    end if;
end function; 

// DONE
gatesInv8 := function(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2)
    if (rr1 eq 1) then
        sum := 2 * gatesAdd4(); 
        sum := sum + gatesSquareScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi);
        sum := sum + (3 * gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2));
        sum := sum + gatesInv4(P, Q, sigma, pr1, pr2, qr1, qr2);

        return sum - 10; 
    else // normal basis 
        sum := 2 * gatesAdd4(); 
        sum := sum + gatesSquareScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi);
        sum := sum + (3 * gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2));
        sum := sum + gatesInv4(P, Q, sigma, pr1, pr2, qr1, qr2);
        return sum - 15; 
    end if;
end function;

gatesAdd8 := function()
    return 8; 
end function;

gatesMult8 := function(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2)
    if (rr1 eq 1) then // polynomial basis
        sum := (4 * gatesAdd4());
        sum := sum + (3 * gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2));
        sum := sum + (gatesScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi));
        return sum;
    else // normal basis
        sum := (4 * gatesAdd4());
        sum := sum + (3 * gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2));
        sum := sum + (gatesScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi));
        return sum;
    end if;
end function;

gatesSquare8 := function(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2)
    if rr1 eq 1 then // polynomial basis
        sum := (2 * gatesSquare4(P, Q, sigma, pr1, pr2, qr1, qr2));
        sum := sum + gatesAdd4();
        sum := sum + gatesScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi);
        return sum;
    else // normal basis
        sum := (3 * gatesAdd4());
        sum := sum + (2 * gatesSquare4(P, Q, sigma, pr1, pr2, qr1, qr2));
        sum := sum + gatesScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi);
        return sum;
    end if;
end function;

gatesSquareScale8 := function(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2, lambda)
    e3 := Eltseq(lambda)[2];
    e4 := Eltseq(lambda)[1];
    if (rr1 eq 1) then
        if (e4 eq 0) and (e3 eq pi^(-1)) then
            sum := 2 * gatesSquare4(P, Q, sigma, pr1, pr2, qr1, qr2);
            sum := sum + (2 * gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2));
            sum := sum + (2 * gatesAdd4());
            sum := sum + gatesScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi);
            return sum;
        elif (e3 eq pi^(-1)) then
            sum := 5 * gatesAdd4();
            sum := sum + (2 * gatesSquare4(P, Q, sigma, pr1, pr2, qr1, qr2));
            sum := sum + (3 * gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2));
            sum := sum + gatesScale4(P, Q, sigma, pr1, pr2, qr1, qr2, pi);
            return sum;
        elif (e4 eq 0 and e3 eq pi^(-1)) then
            sum := 2 * gatesAdd4();
            sum := sum + (2 * gatesSquare4(P, Q, sigma, pr1, pr2, qr1, qr2));
            sum := sum + gatesMult4(P, Q, sigma, pr1, pr2, qr1, qr2);
            return sum;
        else // regular multiplication
            return gatesMult8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2);
        end if;
    else
        return gatesMult8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2);
    end if;
end function;

/// GATE COUNTS FOR ARITHMETIC IN GF(2^16)/GF(2^8)

gatesInv16 := function(P, Q, R, S, sigma, pi, lambda, pr1, pr2, qr1, qr2, rr1, rr2, sr1, sr2)
    if (sr1 eq 1) then
        sum := 2 * gatesAdd8(); 
        sum := sum + gatesSquareScale8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2, lambda);
        sum := sum + (3 * gatesMult8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2));
        sum := sum + gatesInv8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2);
        return sum - 20;
    else // normal basis 
        sum := 2 * gatesAdd8();
        sum := sum + gatesSquareScale8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2, lambda);
        sum := sum + (3 * gatesMult8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2));
        sum := sum + gatesInv8(P, Q, R, sigma, pi, pr1, pr2, qr1, qr2, rr1, rr2);
        return sum - 30; 
    end if;
end function;

pad := function(S, n)
    for i := 1 to n do
        S := Insert(S, 0, 0);
    end for;
    return S;
end function;

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
    tmp := 
        [
            Eltseq(elem[2])[2], 
            Eltseq(elem[2])[1], 
            Eltseq(elem[1])[2], 
            Eltseq(elem[1])[1]
        ];
    return tmp;
end function;

buildRoot8Row := function(element, subfield) 
    if element eq 1 then
        return [0,0,0,0,0,0,0,1];
    end if;
    elem := Eltseq(element, subfield);
    if #elem eq 1 then
        tmp := Eltseq(elem[1]); 
        tmp1 := Eltseq(tmp)[1]; // lower
        tmp2 := Eltseq(tmp)[2]; // upper
        row := 
            [
                0,0,0,0,
                Eltseq(tmp2)[2], Eltseq(tmp2)[1], 
                Eltseq(tmp1)[2], Eltseq(tmp1)[1]
            ];
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

buildRoot16Row := function(element, F16, F256)
    if element eq 1 then
        return [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1];
    end if;
    elem := Eltseq(element, F256);
    if #elem eq 1 then // no upper 8 bits
        lower := buildRoot8Row(elem[1], F16);
        row := 
            [0,0,0,0,0,0,0,0,
            lower[1], lower[2],
            lower[3], lower[4],
            lower[5], lower[6],
            lower[7], lower[8]
            ];
        return row;
    else
        upper := buildRoot8Row(elem[2], F16);
        lower := buildRoot8Row(elem[1], F16);

        row := 
            [
            upper[1], upper[2],
            upper[3], upper[4],
            upper[5], upper[6],
            upper[7], upper[8],
            lower[1], lower[2],
            lower[3], lower[4],
            lower[5], lower[6],
            lower[7], lower[8]
            ];
        return row;
    end if;
end function;

// sigma in GF(2^2)
changeSigmaRoot := function(pr1, pr2, F4, sigma) 
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
changeSigmaRoot_bin := function(pr1, pr2, F4, sigma) 
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
changePiRoot := function(pr1, pr2, qr1, qr2, F4, F16, pi) 
    q1p1 := pr1 * qr1;
    q2p1 := pr1 * qr2;
    q1p2 := pr2 * qr1;
    q2p2 := pr2 * qr2;
    r1 := buildRoot4Row(q2p2, F4);
    r2 := buildRoot4Row(q2p1, F4);
    r3 := buildRoot4Row(q1p2, F4);
    r4 := buildRoot4Row(q1p1, F4);
    M := Matrix(GF(2), [ r1, r2, r3, r4 ]);
    evector := 
        [
            Eltseq(Eltseq(pi)[2])[2], 
            Eltseq(Eltseq(pi)[2])[1], 
            Eltseq(Eltseq(pi)[1])[2], 
            Eltseq(Eltseq(pi)[1])[1]
        ];
    ev := Transpose(Matrix(GF(2), [evector])); 
    prod := Transpose(M)^(-1) * ev;
    newPi := Seqelt(
        [
            Seqelt([prod[4][1],prod[3][1]], F4), 
            Seqelt([prod[2][1],prod[1][1]], F4)
        ], F16);
    return newPi;
end function;

// pi in GF(2^4)/GF(2^2)
changePiRoot_bin := function(pr1, pr2, qr1, qr2, F4, F16, pi) 
    q1p1 := pr1 * qr1;
    q2p1 := pr1 * qr2;
    q1p2 := pr2 * qr1;
    q2p2 := pr2 * qr2;
    r1 := buildRoot4Row(q2p2, F4);
    r2 := buildRoot4Row(q2p1, F4);
    r3 := buildRoot4Row(q1p2, F4);
    r4 := buildRoot4Row(q1p1, F4);
    M := Matrix(GF(2), [ r1, r2, r3, r4 ]);
    evector := 
        [
            Eltseq(Eltseq(pi)[2])[2], 
            Eltseq(Eltseq(pi)[2])[1], 
            Eltseq(Eltseq(pi)[1])[2], 
            Eltseq(Eltseq(pi)[1])[1]
        ];
    ev := Transpose(Matrix(GF(2), [evector])); 
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

    M := Matrix(GF(2), [ r1, r2, r3, r4, r5, r6, r7, r8 ]);

    upper := Eltseq(lambda)[2];
    lower := Eltseq(lambda)[1];
    p4 := Eltseq(upper)[2];
    p3 := Eltseq(upper)[1];
    p2 := Eltseq(lower)[2];
    p1 := Eltseq(lower)[1];
    evector := 
        [
            Eltseq(p4)[2], Eltseq(p4)[1], 
            Eltseq(p3)[2], Eltseq(p3)[1], 
            Eltseq(p2)[2], Eltseq(p2)[1], 
            Eltseq(p1)[2], Eltseq(p1)[1]
        ];
    ev := Transpose(Matrix(GF(2), [evector])); 
    prod := Transpose(M)^(-1) * ev;

    return Transpose(prod);
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

    M := Matrix(GF(2), [ r1, r2, r3, r4, r5, r6, r7, r8 ]);

    upper := Eltseq(lambda)[2];
    lower := Eltseq(lambda)[1];
    p4 := Eltseq(upper)[2];
    p3 := Eltseq(upper)[1];
    p2 := Eltseq(lower)[2];
    p1 := Eltseq(lower)[1];
    evector := 
        [
            Eltseq(p4)[2], Eltseq(p4)[1], 
            Eltseq(p3)[2], Eltseq(p3)[1], 
            Eltseq(p2)[2], Eltseq(p2)[1], 
            Eltseq(p1)[2], Eltseq(p1)[1]
        ];
    ev := Transpose(Matrix(GF(2), [evector])); 
    prod := Transpose(M)^(-1) * ev;

    newLambda := Seqelt( 
        [
            Seqelt( 
                [ 
                    Seqelt([prod[8][1],prod[7][1]], F4), 
                    Seqelt([prod[6][1],prod[5][1]], F4) 
                ], F16), 
            Seqelt( 
                [ 
                    Seqelt([prod[4][1],prod[3][1]], F4), 
                    Seqelt([prod[2][1],prod[1][1]], F4) 
                ], F16)
        ], F256);
    return newLambda;
end function;

// psi in GF(2^16)/GF(2^8)/GF(2^4)/GF(2^2)
changePsiRoot_matrix := function(pr1, pr2, qr1, qr2, rr1, rr2, sr1, sr2, F4, F16, F256, F6K)
    p1q1r1s1 := pr1 * qr1 * rr1 * sr1;
    p1q2r1s1 := pr1 * qr2 * rr1 * sr1;
    p2q1r1s1 := pr2 * qr1 * rr1 * sr1;
    p2q2r1s1 := pr2 * qr2 * rr1 * sr1;
    p1q1r2s1 := pr1 * qr1 * rr2 * sr1;
    p1q2r2s1 := pr1 * qr2 * rr2 * sr1;
    p2q1r2s1 := pr2 * qr1 * rr2 * sr1;
    p2q2r2s1 := pr2 * qr2 * rr2 * sr1;

    p1q1r1s2 := pr1 * qr1 * rr1 * sr2;
    p1q2r1s2 := pr1 * qr2 * rr1 * sr2;
    p2q1r1s2 := pr2 * qr1 * rr1 * sr2;
    p2q2r1s2 := pr2 * qr2 * rr1 * sr2;
    p1q1r2s2 := pr1 * qr1 * rr2 * sr2;
    p1q2r2s2 := pr1 * qr2 * rr2 * sr2;
    p2q1r2s2 := pr2 * qr1 * rr2 * sr2;
    p2q2r2s2 := pr2 * qr2 * rr2 * sr2;

    r1 := buildRoot16Row(F6K ! p2q2r2s2, F16, F256);
    r2 := buildRoot16Row(F6K ! p1q2r2s2, F16, F256);
    r3 := buildRoot16Row(F6K ! p2q1r2s2, F16, F256);
    r4 := buildRoot16Row(F6K ! p1q1r2s2, F16, F256);
    r5 := buildRoot16Row(F6K ! p2q2r1s2, F16, F256);
    r6 := buildRoot16Row(F6K ! p1q2r1s2, F16, F256);
    r7 := buildRoot16Row(F6K ! p2q1r1s2, F16, F256);
    r8 := buildRoot16Row(F6K ! p1q1r1s2, F16, F256);

    r9  := buildRoot16Row(F6K ! p2q2r2s1, F16, F256);
    r10 := buildRoot16Row(F6K ! p1q2r2s1, F16, F256);
    r11 := buildRoot16Row(F6K ! p2q1r2s1, F16, F256);
    r12 := buildRoot16Row(F6K ! p1q1r2s1, F16, F256);
    r13 := buildRoot16Row(F6K ! p2q2r1s1, F16, F256);
    r14 := buildRoot16Row(F6K ! p1q2r1s1, F16, F256);
    r15 := buildRoot16Row(F6K ! p2q1r1s1, F16, F256);
    r16 := buildRoot16Row(F6K ! p1q1r1s1, F16, F256);

    M := Matrix(GF(2), [ r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16 ]);

    return Transpose(M);
end function;

// basis change matrices using  in GF(2^8)/GF(2^4)/GF(2^2)
totalGateCount8 := function(invCount, v, w, x, prr1, prr2, qrr1, qrr2, rrr1, rrr2, 
    F4, F16, F256, field, AFFINE, C, CINV)

    // Map to isomorphic elements in GF(2^8) 
    pr1 := field ! 1;
    pr2 := field ! v;
    qr1 := field ! 1;
    qr2 := field ! w;
    rr1 := field ! 1;
    rr2 := field ! x;

    // Build the basis change matrix rows
    p1q1r1 := Reverse(Eltseq(pr1 * qr1 * rr1));
    p1q1r1 := pad(p1q1r1, 8 - #p1q1r1);
    p1q2r1 := Reverse(Eltseq(pr1 * qr2 * rr1));
    p1q2r1 := pad(p1q2r1, 8 - #p1q2r1);
    p2q1r1 := Reverse(Eltseq(pr2 * qr1 * rr1));
    p2q1r1 := pad(p2q1r1, 8 - #p2q1r1);
    p2q2r1 := Reverse(Eltseq(pr2 * qr2 * rr1));
    p2q2r1 := pad(p2q2r1, 8 - #p2q2r1);
    p1q1r2 := Reverse(Eltseq(pr1 * qr1 * rr2));
    p1q1r2 := pad(p1q1r2, 8 - #p1q1r2);
    p1q2r2 := Reverse(Eltseq(pr1 * qr2 * rr2));
    p1q2r2 := pad(p1q2r2, 8 - #p1q2r2);
    p2q1r2 := Reverse(Eltseq(pr2 * qr1 * rr2));
    p2q1r2 := pad(p2q1r2, 8 - #p2q1r2);
    p2q2r2 := Reverse(Eltseq(pr2 * qr2 * rr2));
    p2q2r2 := pad(p2q2r2, 8 - #p2q2r2);

    r1 := p2q2r2;
    r2 := p1q2r2;
    r3 := p2q1r2;
    r4 := p1q1r2;
    r5 := p2q2r1;
    r6 := p1q2r1;
    r7 := p2q1r1;
    r8 := p1q1r1;

    // Create the basis change matrix to go from GF(2^8) to GF(((2^2)^2)^2)
    M := Transpose(Matrix(GF(2), [  r1, r2, r3, r4, r5, r6, r7, r8 ])); // M == X
    MI := M^(-1);

    // Create the basis change in GF(((2^2)^2)^2)
    M2 := changeLambdaRoot_matrix(prr1, prr2, qrr1, qrr2, rrr1, rrr2, F4, F16, F256);
    M := M * M2;
    MI := M2^(-1) * MI;

    // Display the elements that are needed for the isomorphism...
    field ! v;
    field ! w;
    field ! x;

    // Display all of the basis change matrices and the 
    // important combinations
    // merged encryption sbox
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
    return 0; 
end function;

totalGateCount16 := function(invCount, v, w, x, y, prr1, prr2, qrr1, qrr2, rrr1, rrr2, srr1, srr2, 
  F4, F16, F256, F6K, field, AFFINE, C, CINV)

    // Map to isomorphic elements in GF(2^8)
    pr1 := field ! 1;
    pr2 := field ! v;
    qr1 := field ! 1;
    qr2 := field ! w;
    rr1 := field ! 1;
    rr2 := field ! x;
    sr1 := field ! 1;
    sr2 := field ! y;

    // Build the basis change matrix rows
    p1q1r1sr1 := Reverse(Eltseq(pr1 * qr1 * rr1 * sr1));
    p1q1r1sr1 := pad(p1q1r1sr1, 16 - #p1q1r1sr1);
    p1q2r1sr1 := Reverse(Eltseq(pr1 * qr2 * rr1 * sr1));
    p1q2r1sr1 := pad(p1q2r1sr1, 16 - #p1q2r1sr1);
    p2q1r1sr1 := Reverse(Eltseq(pr2 * qr1 * rr1 * sr1));
    p2q1r1sr1 := pad(p2q1r1sr1, 16 - #p2q1r1sr1);
    p2q2r1sr1 := Reverse(Eltseq(pr2 * qr2 * rr1 * sr1));
    p2q2r1sr1 := pad(p2q2r1sr1, 16 - #p2q2r1sr1);
    p1q1r2sr1 := Reverse(Eltseq(pr1 * qr1 * rr2 * sr1));
    p1q1r2sr1 := pad(p1q1r2sr1, 16 - #p1q1r2sr1);
    p1q2r2sr1 := Reverse(Eltseq(pr1 * qr2 * rr2 * sr1));
    p1q2r2sr1 := pad(p1q2r2sr1, 16 - #p1q2r2sr1);
    p2q1r2sr1 := Reverse(Eltseq(pr2 * qr1 * rr2 * sr1));
    p2q1r2sr1 := pad(p2q1r2sr1, 16 - #p2q1r2sr1);
    p2q2r2sr1 := Reverse(Eltseq(pr2 * qr2 * rr2 * sr1));
    p2q2r2sr1 := pad(p2q2r2sr1, 16 - #p2q2r2sr1);
    p1q1r1sr2 := Reverse(Eltseq(pr1 * qr1 * rr1 * sr2));
    p1q1r1sr2 := pad(p1q1r1sr2, 16 - #p1q1r1sr2);
    p1q2r1sr2 := Reverse(Eltseq(pr1 * qr2 * rr1 * sr2));
    p1q2r1sr2 := pad(p1q2r1sr2, 16 - #p1q2r1sr2);
    p2q1r1sr2 := Reverse(Eltseq(pr2 * qr1 * rr1 * sr2));
    p2q1r1sr2 := pad(p2q1r1sr2, 16 - #p2q1r1sr2);
    p2q2r1sr2 := Reverse(Eltseq(pr2 * qr2 * rr1 * sr2));
    p2q2r1sr2 := pad(p2q2r1sr2, 16 - #p2q2r1sr2);
    p1q1r2sr2 := Reverse(Eltseq(pr1 * qr1 * rr2 * sr2));
    p1q1r2sr2 := pad(p1q1r2sr2, 16 - #p1q1r2sr2);
    p1q2r2sr2 := Reverse(Eltseq(pr1 * qr2 * rr2 * sr2));
    p1q2r2sr2 := pad(p1q2r2sr2, 16 - #p1q2r2sr2);
    p2q1r2sr2 := Reverse(Eltseq(pr2 * qr1 * rr2 * sr2));
    p2q1r2sr2 := pad(p2q1r2sr2, 16 - #p2q1r2sr2);
    p2q2r2sr2 := Reverse(Eltseq(pr2 * qr2 * rr2 * sr2));
    p2q2r2sr2 := pad(p2q2r2sr2, 16 - #p2q2r2sr2);

    r1 := p2q2r2sr2;
    r2 := p1q2r2sr2;
    r3 := p2q1r2sr2;
    r4 := p1q1r2sr2;
    r5 := p2q2r1sr2;
    r6 := p1q2r1sr2;
    r7 := p2q1r1sr2;
    r8 := p1q1r1sr2;
    r9 := p2q2r2sr1;
    r10 := p1q2r2sr1;
    r11 := p2q1r2sr1;
    r12 := p1q1r2sr1;
    r13 := p2q2r1sr1;
    r14 := p1q2r1sr1;
    r15 := p2q1r1sr1;
    r16 := p1q1r1sr1;

    // Create the inverse basis change matrix
    M := Transpose(Matrix(GF(2), 
      [  
        r1, r2, r3, r4, r5, r6, r7, r8, 
        r9, r10, r11, r12, r13, r14, r15, r16 
      ])); // M == X
    MI := M^(-1);

    // Create the basis change in GF((((2^2)^2)^2)^2)
    M2 := changePsiRoot_matrix(prr1, prr2, qrr1, qrr2, rrr1, rrr2, srr1, srr2, F4, F16, F256, F6K);
    M := M * M2;
    MI := M2^(-1) * MI;

    // Display the elements that are needed for the isomorphism...
    field ! v;
    field ! w;
    field ! x;
    field ! y;

    // Display all of the basis change matrices and the important combinations
    // merged encryption sbox
    MI; // T^(-1)
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

coeffMap2 := function(sigma, F, F4) 
    // these are the only two possibilities for p(v) to be irreducible
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
    // The mapping for each pi coefficient is defined as follows:
    // [0,0] (0)       -> [0,0]
    // [0,1] (v)       -> sigma if sigma = [0,1], else sigma^2
    // [1,0] (v^2)     -> sigma^2 if sigma = [0,1], else sigma
    // [1,1] (v^2+v)   -> 1 (v^2 + v \equiv 1 by p(v))
    c1 := F4 ! 0; // upper coefficient (c2 w^4 + c1 w) -> [c1, c2]
    c2 := F4 ! 0; // lower coefficient (c2 w^4 + c1 w) -> [c1, c2]
    case Eltseq(sigma):
        when [0, 1]: // sigma (v)
            case Eltseq(Eltseq(pi)[1]): // c1
                when [1,0]: 
                    c1 := sigma;
                when [1,1]:
                    c1 := F4 ! 1;
                when [0,1]:
                    c1 := sigma^2;
            end case;
            case Eltseq(Eltseq(pi)[2]): // c2
                when [1,0]:
                    c2 := sigma;
                when [1,1]:
                    c2 := F4 ! 1;
                when [0,1]:
                    c2 := sigma^2;
            end case;
        when [1,1]: // sigma^2 (v^2)
            case Eltseq(Eltseq(pi)[1]): // c1
                when [1,0]:
                    c1 := sigma^2;
                when [1,1]:
                    c1 := F4 ! 1;
                when [0,1]:
                    c1 := sigma;
            end case;
            case Eltseq(Eltseq(pi)[2]): // c2
                when [1,0]:
                    c2 := sigma^2;
                when [1,1]:
                    c2 := F4 ! 1;
                when [0,1]:
                    c2 := sigma;
            end case;
    end case;
    return Seqelt([c1, c2], F16);
end function;

// Generate all tuples of basis elements (p, r, q, s)
allGen_8 := function(embed, field, S, AFFINE, C, CINV)
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

    if pRoots[1] gt pRoots[2] then
        tmp := pRoots[1];
        pRoots[1] := pRoots[2];
        pRoots[2] := tmp;
    end if;

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
                    if embed eq 0 then
                        qRoots:=Append(qRoots,e1);
                    else
                        qRoots:=Append(qRoots,field!e1);
                    end if;
                end if;
            end for;

            if qRoots[1] gt qRoots[2] then
                tmp := qRoots[1];
                qRoots[1] := qRoots[2];
                qRoots[2] := tmp;
            end if;

            // Find all elements in GF(2^4)/GF(2^2) that make r(y) irreducible
            for pi in F16 do
                pol8<X>:=PolynomialRing(F16);
                R:=X^2 + X + pi;

                if IsIrreducible(R) then
                    F256<x>:=ext<F16 | R>;

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

                    if rRoots[1] gt rRoots[2] then
                        tmp := rRoots[1];
                        rRoots[1] := rRoots[2];
                        rRoots[2] := tmp;
                    end if;

                    newPi := changePiRoot(pRoots[1], pRoots[2], 1, qRoots[1], F4, F16, pi);
                    newSigma := changeSigmaRoot(pRoots[1], pRoots[2], F4, sigma); 
                    newSigma := coeffMap2(newSigma, F, F4); 
                    newPi := coeffMap4(newSigma, newPi, F4, F16); 
                    invCount := canrightInv8(P, Q, newSigma, pRoots[1], pRoots[2], 1, qRoots[1], newPi);
                    toss := totalGateCount8(invCount, v, w, x, pRoots[1], pRoots[2], 1, qRoots[1], 

                    // OTHER CASES OMITTED FOR BREVITY

                end if;
            end for;
        end if;
    end for;
    return [* pIps, qIps, rIps *];
end function;

allGen_16 := function(embed, field, T, AFFINE, C, CINV) 
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
    for sigma in F4 do
        pol4<W>:=PolynomialRing(F4);
        Q:=W^2 + W + sigma;
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
            for pi in F16 do
                pol8<X>:=PolynomialRing(F16);
                R:=X^2 + X + pi;

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

                    for lambda in F256 do
                        pol16<Y> := PolynomialRing(F256);
                        S := Y^2 + Y + lambda;

                        if IsIrreducible(S) and S notin sIps then
                            F6K<y> := ext<F256 | S>;
                            // sIps:=Append(sIps, S);

                            try
                                F6K;
                                field;
                                Embed(F6K, field);
                            catch e
                                print e;
                                print("Embedding already done. No harm no foul.");
                            end try;

                            try
                                sRoots := [];
                                for e3 in F6K do 
                                    if Evaluate(S, e3) eq 0 then
                                        if embed eq 0 then
                                            sRoots:=Append(sRoots,e3);
                                        else
                                            sRoots:=Append(sRoots,field!e3);
                                        end if;
                                    end if;
                                end for;
                            catch e
                                e;
                            end try;


                            newPi := changePiRoot(1,pRoots[1],1,qRoots[1], F4, F16, pi);
                            newSigma := changeSigmaRoot(1,pRoots[1], F4, sigma);
                            invCount := gatesInv16(P, Q, R, S, newSigma, newPi, lambda, 1,pRoots[1],1,
                                qRoots[1],1,rRoots[1],1,sRoots[1]);
                            toss := totalGateCount16(invCount,v,w,x,y, 1,pRoots[1],1,qRoots[1],1,
                                rRoots[1],1,sRoots[1], F4, F16, F256, F6K, field, AFFINE, C, CINV);

                            // OTHER CASES OMITTED FOR BREVITY

                        end if;
                    end for;
                end if;
            end for;
        end if;
    end for;

    return [* pIps, qIps, rIps, sIps *];
end function;

////////////////////////////////////////////////////////////
/////////////// START AES ALTERNATIVE SEARCH ///////////////
////////////////////////////////////////////////////////////
F:=GF(2);
polRing<V>:=PolynomialRing(F);
S := V^8 + V^4 + V^3 + V + 1; // fixed affine for AES polynomial
F256<x>:=ext<F | S>;
if PrimitiveElement(F256) ne x + 1 then
    quit;
end if;
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
cinv := constantAffineInv(affine, constant, F256);
p := allGen_8(0, F256, S, affine, constant, cinv); 
////////////////////////////////////////////////////////////
///////////////  END AES ALTERNATIVE SEARCH  ///////////////
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
///////////////   START THESIS TEST CASE 1   ///////////////
////////////////////////////////////////////////////////////
F2 := GF(2);
Poly2<V> := PolynomialRing(F2);
P := V^2 + V + 1;
F4<v> := ext<F2 | P>;
Poly4<W> := PolynomialRing(F4);
Q := W^2 + W + v;
F16<w> := ext<F4 | Q>;
Poly16<X> := PolynomialRing(F16);
R := X^2 + X + (v + 1)*w + v;
F256<x> := ext<F16 | R>;
newSigma := changeSigmaRoot(1, v, F4, v);
newPi := changePiRoot(1, v, w, w^4, F4, F16, (v + 1)*w + v);
gatesInv8(P, Q, R, newSigma, newPi, 1, v, w, w^4, x, x^16);
////////////////////////////////////////////////////////////
///////////////    END THESIS TEST CASE 1    ///////////////
////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////
///////////////   START THESIS TEST CASE 2   ///////////////
////////////////////////////////////////////////////////////
F2 := GF(2);
Poly2<V> := PolynomialRing(F2);
P := V^2 + V + 1;
F4<v> := ext<F2 | P>;
Poly4<W> := PolynomialRing(F4);
Q := W^2 + W + v;
F16<w> := ext<F4 | Q>;
Poly16<X> := PolynomialRing(F16);
R := X^2 + X + ((v + 1)*w + v);
F256<x> := ext<F16 | R>;
Poly256<Y> := PolynomialRing(F256);
S := Y^2 + Y + (v*w + v)*x + w;
F6K<y> := ext<F256 | S>;
gatesInv16(P, Q, R, S, v, ((v + 1)*w + v), (v*w*x + (w + v)), 1, v, 1, w, 1, x, 1, y);
////////////////////////////////////////////////////////////
///////////////    END THESIS TEST CASE 2    ///////////////
////////////////////////////////////////////////////////////
