clear;

q := 2^64 - 2^32 + 1;
IsPrime(q);
Fq := GF(q);

// 512-th primitive roots of unity mod q
// needed to multipy in Z[x]/(x^256+1)

omega := Fq ! 237214921853999334;

/* Dilithium params

q := 2^23 - 2^13  + 1;
Fq := GF(q);
omega := Fq ! 1753;

*/


function reverseIndex(i, b)

  a := IntegerToSequence(i, 2);
  if #a lt b then
    a := a cat [Integers() ! 0 : k in [1..(b-#a)]];
  end if;
  
  return SequenceToInteger(Reverse(a), 2);
  
end function;

// bit reverse A for use in in place FFT

function bitReverse(A)  // normally have precomputed bitreverse array
  n := Round(Log(2, #A));
  return [A[reverseIndex(k, n)+1] : k in [0..#A-1]]; 
end function;

function MapTuplesToArray(A)

  B := A[1];
  for i := 2 to #A do
    B  := B cat A[i];
  end for;

  return B;

end function;

zetas := bitReverse([omega^i : i in [0..255]]);

procedure fNTT(~a) 

  k := 0;
  len := 128;
  N := 256;
  while (len gt 0) do
    for start := 0 to N-1 by 2*len do
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      for j := start to (start + len - 1) do
        t := zeta * a[j + len + 1];
        a[j + len + 1] := a[j + 1] - t;
        a[j + 1] := a[j + 1] + t;
      end for;
    end for;
	len := len div 2;
  end while;
end procedure;


// a is now array of tuples
// the original a gets mapped to [[a[i], a[i+128]]

procedure fNTT_2elems(~a) 

  k := 0;
  len := 128;
  N := 256;
  while (len gt 1) do
    for start := 0 to (N div 2)-1 by len do
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      for j := start to (start + (len div 2) - 1) do
        // reading in tuple 1
		A1 := a[j+1];
		// reading in tuple 2
		A2 := a[j+1 + (len div 2)];
		
		t1 := zeta * A1[2];
        t2 := zeta * A2[2];
		
		a[j+1] := [A1[1] + t1, A2[1] + t2];
        a[j+1 + (len div 2)] := [A1[1] - t1, A2[1] - t2];
	    	
      end for;
    end for;
	len := len div 2;
  end while;
  
  // doing len = 1
  
  for j := 0 to (N div 2)-1 do
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      	A1 := a[j+1];
		t1 := zeta * A1[2];
		a[j+1] := [A1[1] + t1, A1[1] - t1];
  end for;
  
end procedure;


printUnrolled := false;
procedure printU(s)

  if printUnrolled then print s; end if;

end procedure;



// a is now array of quadruples
// the original a gets mapped to [[a[i], a[i+128], a[i+1], a[i+1+128]]

procedure fNTT_4elems(~a) 

RW := 0;
AD := 0;
MU := 0;


  k := 0;
  len := 128;
  N := 256;
  while (len gt 2) do
    for start := 0 to (N div 4)-1 by (len div 2) do
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      for j := start to (start + (len div 4) - 1) do
        // reading in tuple 1
		A1 := a[j+1];
		printU("A1 := a[" cat IntegerToString(j+1) cat "];");
		
		// reading in tuple 2
		A2 := a[j+1 + (len div 4)];
		printU("A2 := a[" cat IntegerToString(j+1 + (len div 4)) cat "];");
		
		t1 := zeta * A1[2];
		printU("t1 := " cat IntegerToString(Integers() ! zeta) cat "*A1[2];");
		
        t2 := zeta * A2[2];
		printU("t2 := " cat IntegerToString(Integers() ! zeta) cat "*A2[2];");
		
		t3 := zeta * A1[4];
		printU("t3 := " cat IntegerToString(Integers() ! zeta) cat "*A1[4];");
		
		t4 := zeta * A2[4];
		printU("t4 := " cat IntegerToString(Integers() ! zeta) cat "*A2[4];");
		
		a[j+1] := [A1[1] + t1, A2[1] + t2, A1[3] + t3, A2[3] + t4];
		printU("a[" cat IntegerToString(j+1) cat "] :=  [A1[1] + t1, A2[1] + t2, A1[3] + t3, A2[3] + t4];");
	
        a[j+1 + (len div 4)] := [A1[1] - t1, A2[1] - t2, A1[3] - t3, A2[3] - t4];
	    printU("a[" cat IntegerToString(j+1 + (len div 4)) cat "] :=  [A1[1] - t1, A2[1] - t2, A1[3] - t3, A2[3] - t4];");
		
		RW +:= 4;
		MU +:= 4;
		AD +:= 8;
		
      end for;
    end for;
	len := len div 2;
  end while;

  // doing len = 2 
  
   for j := 0 to (N div 4)-1 do
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      // reading in tuple 1
	  A1 := a[j+1];
	  printU("A1 := a[" cat IntegerToString(j+1) cat "];");
	  
	  t1 := zeta * A1[2];
      t2 := zeta * A1[4];
	  
	  printU("t1 := " cat IntegerToString(Integers() ! zeta) cat "*A1[2];");
	  printU("t2 := " cat IntegerToString(Integers() ! zeta) cat "*A1[4];");
		
	  a[j+1] := [A1[1] + t1, A1[3] + t2, A1[1] - t1, A1[3] - t2];
	  printU("a[" cat IntegerToString(j+1) cat "] :=  [A1[1] + t1, A1[3] + t2, A1[1] - t1, A1[3] - t2];");
	  
	  	RW +:= 2;
		MU +:= 2;
		AD +:= 4;
		
       
    end for;
   
  // doing len = 1

  for j := 0 to (N div 4)-1 do
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      A1 := a[j+1];
	  
	  printU("A1 := a[" cat IntegerToString(j+1) cat "];");
	  
	  t1 := zeta * A1[2];
	  printU("t1 := " cat IntegerToString(Integers() ! zeta) cat "*A1[2];");
	  
	  k := k+1;
      zeta := zetas[k+1];  // index starting at 1
      t2 := zeta * A1[4];
	  printU("t2 := " cat IntegerToString(Integers() ! zeta) cat "*A1[4];");
	 
	  a[j+1] := [A1[1] + t1, A1[1] - t1, A1[3] + t2, A1[3] - t2];
	  
	  printU("a[" cat IntegerToString(j+1) cat "] :=  [A1[1] + t1, A1[1] - t1, A1[3] + t2, A1[3] - t2];");
	  
	  	RW +:= 2;
		MU +:= 2;
		AD +:= 4;
	  
  end for;
  
  print RW, MU, AD;
  
end procedure;



procedure invNTT(~a)
  
  k := 256;
  len := 1;
  N := 256;
  
  while (len lt N) do
    for start := 0 to N-1 by 2*len do
	  k := k-1;
	  zeta := -zetas[k+1];  // indexing from 1
      for j := start to start + len-1 do
        t := a[j + 1];
        a[j + 1] := t + a[j + len + 1];
        a[j + len + 1] := t - a[j + len + 1];
        a[j + len + 1] := zeta * a[j + len + 1];
      end for;
     end for;
    len := 2*len;
  end while;
 
  for j := 0 to N-1 do
	a[j+1] := a[j+1]/N;
  end for;

end procedure;


// input is now array of tuples [a[2k], a[2k+1]]


procedure invNTT_2elem(~a)
  
  k := 256;
  len := 1;
  N := 256;
   
  while (len lt (N div 2)) do
    for start := 0 to (N div 2)-1 by 2*len do
	  k := k-1;
	  zeta1 := -zetas[k+1];  // indexing from 1
	  k := k-1;
	  zeta2 := -zetas[k+1];
      for j := start to start + len-1 do
        A1 := a[j+1];
		A2 := a[j + 1 + len];
		
        u1 := A1[1] + A1[2];
        u2 := zeta1*(A1[1] - A1[2]);
        
		v1 := A2[1] + A2[2];
		v2 := zeta2*(A2[1] - A2[2]);
		
		a[j+1] := [u1, v1];
		a[j+1+len] := [u2, v2];
		
      end for;
     end for;
    len := 2*len;
  end while;
 
 // len = N/2 and scaling by N together
 
 k := k-1;
 zeta := -zetas[k+1];
 
  for j := 0 to (N div 2)-1 do
	 t := a[j + 1][1];
     a[j + 1] := [(t + a[j + 1][2])/N, zeta*(t - a[j+1][2])/N];
  end for;

 // reordering to get linear ordering
 // taking last two elements of 1st row and swapping with first two of second row
 // then interleaving
  
  for j := 0 to (N div 4)-1 do
	 A1 := a[2*j + 1];
     A2 := a[2*j + 2];
	 
	 a[2*j + 1] := [A1[1], A2[1]];
	 a[2*j + 2] := [A1[2], A2[2]];
	
  end for;
  
  // now interleaves so placing in order, could be done using swaps
  
  b := a;
  for j := 0 to (N div 4)-1 do
    b[j + 1] := a[2*j + 1];
	b[j + (N div 4) + 1] := a[2*j + 2];
  end for;
  
  a := b;

end procedure;



procedure invNTT_4elem(~a)
  
  k := 256;
  len := 1;
  N := 256;
   
  // len = 1 separate

  for j := 0 to (N div 4)-1 do
	  k := k-1;
	  zeta1 := -zetas[k+1];  // indexing from 1
	  k := k-1;
	  zeta2 := -zetas[k+1];
      
      A1 := a[j+1];
	
      u1 := A1[1] + A1[2];
      u2 := zeta1*(A1[1] - A1[2]);
        
	  v1 := A1[3] + A1[4];
	  v2 := zeta2*(A1[3] - A1[4]);
		
	  a[j+1] := [u1, v1, u2, v2];
  end for;
  
  len := 2;  
   
  while (len lt (N div 2)) do
    for start := 0 to (N div 4)-1 by len do
	  k := k-1;
	  zeta1 := -zetas[k+1];  // indexing from 1
	  k := k-1;
	  zeta2 := -zetas[k+1];
      for j := start to start + (len div 2)-1 do
        A1 := a[j+1];
		A2 := a[j + 1 + (len div 2)];
		
		// first pairs
		
        u1 := A1[1] + A1[2];
        u2 := zeta1*(A1[1] - A1[2]);
        
		v1 := A2[1] + A2[2];
		v2 := zeta2*(A2[1] - A2[2]);
		
		// second pairs
				
        u3 := A1[3] + A1[4];
        u4 := zeta1*(A1[3] - A1[4]);
        
		v3 := A2[3] + A2[4];
		v4 := zeta2*(A2[3] - A2[4]);
		
		
		a[j+1] := [u1, v1, u3, v3];
		a[j+1+(len div 2)] := [u2, v2, u4, v4];
	
      end for;
     end for;
    len := 2*len;
  end while;
 
 // len = N/2 and scaling by N together
 
 k := k-1;
 zeta := -zetas[k+1];
 
  for j := 0 to (N div 4)-1 do
	 A := a[j + 1];
     a[j + 1] := [(A[1] + A[2])/N, (A[3] + A[4])/N, zeta*(A[1] - A[2])/N, zeta*(A[3] - A[4])/N];
  end for;

 // reordering to get linear ordering
 // taking last two elements of 1st row and swapping with first two of second row
 // then interleaving
  
  for j := 0 to (N div 8)-1 do
	 A1 := a[2*j + 1];
     A2 := a[2*j + 2];
	 
	 a[2*j + 1] := [A1[1], A1[2], A2[1], A2[2]];
	 a[2*j + 2] := [A1[3], A1[4], A2[3], A2[4]];
	
  end for;
  
  // now interleaves so placing in order, could be done using swaps
  
  b := a;
  for j := 0 to (N div 8)-1 do
    b[j + 1] := a[2*j + 1];
	b[j + (N div 8) + 1] := a[2*j + 2];
  end for;
  
  a := b;
  
end procedure;


N := 256;
a := [Random(Fq) : k in [1..N]];

Fqx<x> := PolynomialRing(Fq);
A := Fqx ! a;

Bs := MapTuplesToArray([[Evaluate(A, omega^reverseIndex(128 + i, 8)), Evaluate(A, -omega^reverseIndex(128 + i, 8))]  : i in [0..127]]);

Cs := bitReverse([Evaluate(A, omega^(2*k+1)) : k in [0..255]]);

b := a;


a2s := [ [a[k], a[k+128]] : k in [1..128]];
a4s := [ a2s[2*k + 1] cat a2s[2*k+2] : k in [0..63]];

fNTT_2elems(~a2s); 
fNTT_4elems(~a4s);

fNTT(~a);
print "Forward NTT ", Bs eq a;
print "Forward NTT (2) ", MapTuplesToArray(a2s) eq Bs;
print "Forward NTT (4) ", MapTuplesToArray(a4s) eq Bs;

invNTT(~a);
print "Backward NTT ", a eq b;
invNTT_2elem(~a2s);
print "Backward NTT (2) ", MapTuplesToArray(a2s) eq b;
invNTT_4elem(~a4s);
print "Backward NTT (4) ", MapTuplesToArray(a4s) eq b;



