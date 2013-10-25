/***
 *  Examples of use of the ellbasis package
 *
 *  Distributed under the terms of the GNU Lesser General Public License (L-GPL)
 *                  http://www.gnu.org/licenses/
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 *
 *  Copyright 2007-2008 R. Lercier & J.-M. Couveignes
 */

 /***
 * Changes.
 *
 * 2008-05-27 : R. Lercier
 *              New tests added
 *
 * 2007-12-29 : J.-M. Couveignes, R. Lercier
 *              Initial version.
 *
 ********************************************************************/

// Attach("package/ellbasis.m");

EllBasisTest := function(EB)

    fullret := true;

    d      := EB`d;
    Fq     := EB`Fq;
    Fqd<t> := EB`Fqd;

    ""; " Elliptic Basis Tests :";
    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	A := EllPolynomialToBasis(a, EB);
	ret and:= EllBasisToPolynomial(A, EB) eq a;
    end for;
    " Polynomial Conversion test is OK ?", ret;
    fullret and:= ret;
    
    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	A := EllPolynomialToBasis(a, EB);
	C := EllBasisFrobenius(A, EB);
	ret and:= EllBasisToPolynomial(C, EB) eq a^#Fq;
    end for;
    " Frobenius test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	A := EllPolynomialToBasis(a, EB);
	C := EllBasisFrobeniusInverse(A, EB);
	ret and:= EllBasisToPolynomial(C, EB) eq a^(#Fq^(d-1));
    end for;
    " Inverse Frobenius test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);  b := Random(Fqd);
	A := EllPolynomialToBasis(a, EB);
	B := EllPolynomialToBasis(b, EB);
	C := EllBasisMultiply(A, B, EB);
	ret and:= EllBasisToPolynomial(C, EB) eq a*b;
    end for;
    " Multiplication test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 1 do
	repeat a := Random(Fqd); until a ne 0;
	A := EllPolynomialToBasis(a, EB);
	C := EllBasisInverse(A, EB);
	ret and:= EllBasisToPolynomial(C, EB) eq 1/a;
    end for;
    " Inverse test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 1 do
	a := Random(Fqd); exp := Random(#Fqd);
	A := EllPolynomialToBasis(a, EB);
	C := EllBasisExp(A, exp, EB);
	ret and:= EllBasisToPolynomial(C, EB) eq a^(exp);
    end for;
    " Exponentiation is OK ?", ret;
    fullret and:= ret;

    return fullret;

end function;

EllNormalBasisTest := function(EB)

    fullret := true;
    
    d      := EB`d;
    Fq     := EB`Fq;
    Fqd<t> := EB`Fqd;

    ""; " Normal Elliptic Basis Tests :";

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	A := EllPolynomialToNormalBasis(a, EB);
	ret and:= EllNormalBasisToPolynomial(A, EB) eq a;
    end for;
    " Polynomial Conversion test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	B := EllPolynomialToNormalBasis(a, EB);
	C := EllNormalBasisToBasis(B, EB);
	A := EllBasisToNormalBasis(C, EB);
	ret and:= A eq B;
    end for;
    " Normal to Elliptic conversion test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	A := EllPolynomialToNormalBasis(a, EB);
	C := EllNormalBasisFrobenius(A, EB);
	ret and:= EllNormalBasisToPolynomial(C, EB) eq a^#Fq;
    end for;
    " Frobenius test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);
	A := EllPolynomialToNormalBasis(a, EB);
	C := EllNormalBasisFrobeniusInverse(A, EB);
	ret and:= EllNormalBasisToPolynomial(C, EB) eq a^(#Fq^(d-1));
    end for;
    " Inverse Frobenius test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 10 do
	a := Random(Fqd);  b := Random(Fqd);
	A := EllPolynomialToNormalBasis(a, EB);
	B := EllPolynomialToNormalBasis(b, EB);
	C := EllNormalBasisMultiply(A, B, EB);
	ret and:= EllNormalBasisToPolynomial(C, EB) eq a*b;
    end for;
    " Multiplication test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 1 do
	repeat a := Random(Fqd); until a ne 0;
	A := EllPolynomialToNormalBasis(a, EB);
	C := EllNormalBasisInverse(A, EB);
	ret and:= EllNormalBasisToPolynomial(C, EB) eq 1/a;
    end for;
    " Inverse test is OK ?", ret;
    fullret and:= ret;

    ret := true;
    for i := 1 to 1 do
	a := Random(Fqd); exp := Random(#Fqd);
	A := EllPolynomialToNormalBasis(a, EB);
	C := EllNormalBasisExp(A, exp, EB);
	ret and:= EllNormalBasisToPolynomial(C, EB) eq a^(exp);
    end for;
    " Exponentiation is OK ?", ret;
    fullret and:= ret;

    return fullret;

end function;


EllTest := function(Fq, d)

    printf "Handling an extension of degree %o of a %o\n", d, Fq;
    fullret := true;
    
    EB := EllBasisCreate(Fq, d);
    if EB`Initialised eq false then
	"Sorry, none elliptic basis exists for d =", d;
	return false;
    end if;

    fullret and:= EllBasisTest(EB);

    fullret and:= EllNormalBasisTest(EB);

"";"Current test's ok ?", fullret;

    return fullret;
end function;


fulltest := true;

/* Characteristic 2 examples */
fulltest and:= EllTest(GF(2), 2); "";
fulltest and:= EllTest(GF(2), 3); "";
fulltest and:= EllTest(GF(2), 4); "";
fulltest and:= EllTest(GF(2), 5); "";
fulltest and:= EllTest(GF(2^2), 2); "";
fulltest and:= EllTest(GF(2^2), 3); "";
fulltest and:= EllTest(GF(2^3), 4); "";
fulltest and:= EllTest(GF(2^4), 5); "";

/* Characteristic 3 examples */
fulltest and:= EllTest(GF(3), 2); "";
fulltest and:= EllTest(GF(3), 3); "";
fulltest and:= EllTest(GF(3), 4); "";
fulltest and:= EllTest(GF(3), 5); "";
fulltest and:= EllTest(GF(3), 6); "";
fulltest and:= EllTest(GF(3), 7); "";
fulltest and:= EllTest(GF(3^2), 2); "";
fulltest and:= EllTest(GF(3^2), 3); "";
fulltest and:= EllTest(GF(3^2), 7); "";

/* Characteristic 5 examples */
fulltest and:= EllTest(GF(5), 2); "";
fulltest and:= EllTest(GF(5), 3); "";
fulltest and:= EllTest(GF(5), 10); "";

/* Characteristic 7 examples */
fulltest and:= EllTest(GF(7), 2); "";
fulltest and:= EllTest(GF(7), 10); "";
fulltest and:= EllTest(GF(7), 13); "";

/* Characteristic 53 examples */
fulltest and:= EllTest(GF(53), 10); "";
fulltest and:= EllTest(GF(53), 59); "";

/* A large characteristic example */
fulltest and:= EllTest(GF(100000000000000000039), 7); "";

"Everything's ok ?", fulltest;
