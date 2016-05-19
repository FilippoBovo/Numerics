# Lieb-Liniger Equations

This is a repository for numerical solutions of mathematical equations:




## Fredholm Equation

To solve the Fredholm integral equation of the second kind, we define the following three functions:

- **BoundFunction:** This module takes the function *f* as input and as output gives the function *f* in the interval *[a,b]* and zero otherwise.

	```
	BoundFunction[f_, a_, b_] :=
		Function[
			Piecewise[{
				{0., # < a},
				{0., # > b},
				{f[#], True}}]
		];
	```

- **Fredholm2ndKind:** Gives the numerical solution of a Fredholm equation of the second kind in the interval *[a,b]*,

	<p align="center">
		<img src="Resources/Fredholm2ndKind.png" style="height: 3em;">
	</p>

	This takes as input the extremes of integration a and b, the kernel *K(x,y)*, *g(x)* and the number, *n*, of subdivision of the integration interval which is used in the numerical solution. This method is a [numerical adaptation](http://mathematica.stackexchange.com/questions/11594/integral-equation-numerical-solution-with-ndsolve) of [S. Rahbar and E. Hashemizadeh, *A Computational Approach to the Fredholm Integral Equation of the Second Kind*, Proceedings of the World Congress on Engineering, 2008](http://www.iaeng.org/publication/WCE2008/WCE2008_pp933-937.pdf).

	```
	Options[Fredholm2ndKind] = {Method -> Automatic};
	Fredholm2ndKind[{a_, b_, k_, g_}, n_?IntegerQ, OptionsPattern[]] :=
		Block[{step, SI, GI, KMatrix, W, DMatrix, f, deltaX, delta, fI, ftemp},
    	step = (b - a)/n;
		SI = Range[a, b, step];
		GI = g /@ SI;
		KMatrix = Outer[k, SI, SI];
		W = {step/2}~Join~ConstantArray[step, n - 1]~Join~{step/2};
		DMatrix = DiagonalMatrix[W];
		deltaX[x_?NumericQ] := 
		W.(k[x, #] & /@ SI) - NIntegrate[k[x, y], {y, a, b}]; 
		delta = deltaX /@ SI;
		fI = LinearSolve[IdentityMatrix[n + 1] + (DiagonalMatrix[delta] - KMatrix.DMatrix), GI];
		f = If[OptionValue[Method] === NoInterpolation,
			fI,
			BoundFunction[Interpolation[Transpose@{SI, fI}], a, b]];
		f];
	```

- **Fredholm2ndKindOut:**  Gives the numerical solution of the Fredholm equation in the interval *[c,d]* outside of *[a,b]*.This takes as input the extremes of integration *a* and *b*, the kernel *K(x,y)*, *g(x)*, the number, *n*, of subdivision of the integration interval *[a,b]* which is used in the numerical solution, the extreme of integration *c* and *d* of the interval outside *[a,b]* and the number of subdivisions of *[c,d]*.

	```
	Fredholm2ndKindOut[{a_, b_, k_, g_}, n_?IntegerQ, {c_, d_}, m_?IntegerQ, fIni_: True] :=
		Block[{fInTempi, stepIn, SIni, stepOut, SOuti, GOuti, KMatrixOut, fOuti},
		
		(* Variable inside the interval [a,b] *)
		stepIn = (b - a)/n;
		SIni = Range[a, b, stepIn]; (*i-th component of the interval*)
		fInTempi = If[fIni === True,
			Fredholm2ndKind[{a, b, k, g}, n, Method -> NoInterpolation],
			fIni];
		
		(* Variable and functions outside the interval [a,b] *)
		stepOut = (d - c)/m;
		SOuti = Range[c, d, stepOut]; (*i-th component of the interval*)
		GOuti = g /@ SOuti; (*i-th component of the g*)
		KMatrixOut = Outer[k, SOuti, SIni]; (*Matrix form of k*)
		fOuti = GOuti + stepIn*(KMatrixOut.fInTempi - (KMatrixOut[[All, 1]]*fInTempi[[1]] + KMatrixOut[[All, n + 1]]*fInTempi[[n + 1]])/2);
		BoundFunction[Interpolation[Transpose[{SOuti, fOuti}]], c, d]];
	```

## Numerical Solutions

```
m = 0.5; (*Renormalized chemical potential*)
n = 1000; (*number of discretization*)

gLL[x_] := 1/(2*Pi);

gSolve[\[Lambda]_] := Block[{KLL, Out},
  KLL[x_, y_] := \[Lambda]/(Pi*(\[Lambda]^2 + (x - y)^2));
  Out = Fredholm2ndKind[{-1., 1., KLL, gLL}, n];
  Out]

eLL[x_] := -m + x^2;

eSolve[\[Lambda]_, a_, b_] := 
	Module[{KLL, fIni, fIn, fLeft, fRight, Out},
		
		KLL[x_, y_] := \[Lambda]/(Pi*(\[Lambda]^2 + (x - y)^2));
	  
		fIni = Fredholm2ndKind[{-1., 1., KLL, eLL}, n, Method -> NoInterpolation];
		fIn = Interpolation[Transpose[{Table[i, {i, -1., 1., 2./n}], fIni}]];
		fLeft = Fredholm2ndKindOut[{-1., 1., KLL, eLL}, n, {a, -1.}, n, fIni];
		fRight = Fredholm2ndKindOut[{-1., 1., KLL, eLL}, n, {1., b}, n, fIni];
	  
		Function[
			Piecewise[{
				{fLeft[#], a <= # < -1.},
				{fIn[#], -1. <= # <= 1.},
				{fRight[#], 1. < # <= b},
				{0., True}}]
		]
	
	]
```

## Plots

### Spectrum

```
PlotSpectrum[a_, b_, \[Lambda]_, t_] := 
  Block[{eLLfunction},
	
   (* Lieb-Liniger *)
   eLLfunction = eSolve[\[Lambda], a, b];
   
   (* Plot *)
   Plot[eLLfunction[x], {x, a, b}, Exclusions -> None, PlotLegends -> {"LL"}]];

PlotSpectrum[-1.5, 1.5, 0.01, 0.01]
```

### Density of States

```
PlotDOS[a_, b_, \[Lambda]_, t_] := Block[{gLLfunction},
	
	(* Lieb-Liniger *)
	gLLfunction = gSolve[\[Lambda]];
	
	(* Plot *)
	Plot[gLLfunction[x], {x, a, b}, Exclusions -> None, PlotLegends -> {"LL"}]];

PlotDOS[-1.5, 1.5, 0.01, 0.01]
```