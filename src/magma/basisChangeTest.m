// basis change test

AssertAttribute(FldFin, "PowerPrinting", false);  
F:=GF(2);
polRing<V>:=PolynomialRing(F);
P:=V^2+V+1;
F4<v>:=ext<F | P>; 
p2<W>:=PolynomialRing(F4);
Q:=W^2 + W + v;
F16<w>:=ext<F4 | Q>;

pi := (v+1)*w + v;
// pi:= v*w;
pi;

// Roots, for testing purposes...
p1:=1;
p2:=v;
q1:=w;
q2:=w^4;

// Compute the four products
q1p1 := p1 * q1;
q2p1 := p1 * q2;
q1p2 := p2 * q1;
q2p2 := p2 * q2;

r1:= [Eltseq(Eltseq(q2p2)[2])[2], Eltseq(Eltseq(q2p2)[2])[1], Eltseq(Eltseq(q2p2)[1])[2], Eltseq(Eltseq(q2p2)[1])[1]];
r2:= [Eltseq(Eltseq(q2p1)[2])[2], Eltseq(Eltseq(q2p1)[2])[1], Eltseq(Eltseq(q2p1)[1])[2], Eltseq(Eltseq(q2p1)[1])[1]];
r3:= [Eltseq(Eltseq(q1p2)[2])[2], Eltseq(Eltseq(q1p2)[2])[1], Eltseq(Eltseq(q1p2)[1])[2], Eltseq(Eltseq(q1p2)[1])[1]];
r4:= [Eltseq(Eltseq(q1p1)[2])[2], Eltseq(Eltseq(q1p1)[2])[1], Eltseq(Eltseq(q1p1)[1])[2], Eltseq(Eltseq(q1p1)[1])[1]];

r1;
r2;
r3;
r4;

q2p1;
q1p2;
q1p1;

M := Matrix(GF(2), [  r1, r2, r3, r4 ]);
M;
"";
Transpose(M); // correct

evector := [Eltseq(Eltseq(pi)[2])[2], Eltseq(Eltseq(pi)[2])[1], Eltseq(Eltseq(pi)[1])[2], Eltseq(Eltseq(pi)[1])[1]];
ev := Transpose(Matrix(GF(2), [evector])); // correct
ev;

prod := Transpose(M)^(-1) * ev;
"";
prod;
[prod[2][1],prod[1][1]];


Seqelt([prod[2][1],prod[1][1]], F4);
Seqelt([prod[4][1],prod[3][1]], F4);
Seqelt([Seqelt([prod[4][1],prod[3][1]], F4), Seqelt([prod[2][1],prod[1][1]], F4)], F16);

// Seqelt([prod[1],prod[1]], F4);
// Seqelt([Seqelt(prod, F4), Seqelt(prod, F4)], F16);


// TODO: need to be able to build this matrix at runtime... but how?

// affine:=Matrix(GF(2),4,4,
//     [1,0,0,0,1,1,1,1,
//     1,1,0,0,0,1,1,1,
//     1,1,1,0,0,0,1,1,
//     1,1,1,1,0,0,0,1,
//     1,1,1,1,1,0,0,0,
//     0,1,1,1,1,1,0,0,
//     0,0,1,1,1,1,1,0,
//     0,0,0,1,1,1,1,1]);
