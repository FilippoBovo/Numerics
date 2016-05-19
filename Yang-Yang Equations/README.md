# Yang-Yang Equations

C. N. Yang and C. P. Yang in 1969 published an exact analysis of the thermodynamics of an interacting Bose gas:

- [C. N. Yang and C. P. Yang, J. Math. Phys. 10, 1115 (1969)](http://dx.doi.org/10.1063/1.1664947)

They found that the spectrum and the density of states (DOS) satisfy respectively these equations:

<p align="center">
	<img src="Resources/YangYangSpectrum.png">
</p>

and

<p align="center">
	<img src="Resources/YangYangDOS.png">
</p>

where *λ* is the contact interaction strength, *m* is the chemical potential and *t* is the temperature. All quantities are in dimensionless units. The relations to the quantities in dimensional units are in Section 3.8 of [my PhD thesis](http://etheses.bham.ac.uk/6320/1/Bovo15PhD.pdf).

The Mathematica notebook of the numerical solution is [available in this repository](Yang-Yang Equations.nb).

## Recursive Solution of Spectrum Equation

We solve the equation for the spectrum with a recursive procedure, until the solution reaches a defined precision.

```
eYY[\[Mu]_, t_, \[Lambda]_, a_, b_, n_, precision_] :=
	Block[{step, xi, \[Epsilon]0, \[Epsilon]0i, kc, kci, kT, \[Epsilon]i, \[Epsilon]iprev, converged, \[Epsilon]converged},
	
	converged = 0;
	step = (b - a)/(n - 1);
	xi = Range[a, b, step];
	\[Epsilon]0[k_] := -\[Mu] + k^2;
	\[Epsilon]0i = \[Epsilon]0 /@ xi;
	kc[x_, y_] := \[Lambda]/(Pi*(\[Lambda]^2 + (x - y)^2));
	kci = Outer[kc, xi, xi];
	kT[\[Epsilon]_] := t*Log[1 + Exp[-\[Epsilon]/t]];
	\[Epsilon]converged[\[CapitalDelta]\[Epsilon]_] := If[Abs[\[CapitalDelta]\[Epsilon]] < precision, 1, 0];
	\[Epsilon]i = \[Epsilon]0i;
	
	While[converged === 0,
		\[Epsilon]iprev = \[Epsilon]i;
		(* Trapezoidal rule sum *)
		\[Epsilon]i = \[Epsilon]0i - step*(kci.(kT /@ \[Epsilon]iprev) - (kci[[All, 1]]*kT[\[Epsilon]iprev[[1]]] + kci[[All, n]]*kT[\[Epsilon]iprev[[n]]])/2); 
		converged = Times @@ (\[Epsilon]converged /@ (\[Epsilon]iprev -\[Epsilon]i));
	];
	
	\[Epsilon]i]
```

## DOS equation: Fredholm approach

The equation for the density of states can be casted as a Fredholm Equation of the Second Kind. Using [this numerical approach](../Lieb-Liniger%20Equations) we solve it in the following way:

```
gYY[\[Mu]_, t_, \[Lambda]_, a_, b_, n_, precision_] := Block[{eYYDiscrete, eYYfunction, k, g},
	
	eYYDiscrete = eYY[\[Mu], t, \[Lambda], a, b, n, precision];
	eYYfunction = Interpolation[Transpose[{Table[i, {i, a, b, (b - a)/(n - 1)}], eYYDiscrete}]];
	g[x_] := 1/(2*Pi*(1 + Exp[eYYfunction[x]/t]));
	k[x_, y_] := 1/(1 + Exp[eYYfunction[x]/t])*\[Lambda]/(Pi*(\[Lambda]^2 + (x - y)^2));
	
	Fredholm2ndKind[{a, b, k, g}, n]
	];
```

## Plots

Finally, we plot the solutions of the spectrum and density of states.

### Spectrum

```
PlotSpectrum[a_, b_] := 
  Block[{eLLfunction1, eLLfunction2, eLLfunction3},

	(* Lieb-Liniger *)
	eLLfunction1 = eSolve[0.01, a, b];
	eLLfunction2 = eSolve[0.1, a, b];
	eLLfunction3 = eSolve[1, a, b];

	(* Plot *)
	Plot[{eLLfunction1[x], eLLfunction2[x], eLLfunction3[x]}, {x, a, b}, Exclusions -> None, PlotLegends -> {0.01, 0.1, 1}]];

PlotSpectrum[-1.5, 1.5]
```

This gives the spectrum for three different values of λ:

<p align="center">
	<img src="Resources/Spectrum.gif">
</p>

### Density of States

```
PlotDOS[a_, b_] := Block[{gLLfunction1, gLLfunction2, gLLfunction3},

	(* Lieb-Liniger *)
	gLLfunction1 = gSolve[0.01];
	gLLfunction2 = gSolve[0.1];
	gLLfunction3 = gSolve[1];
	
	(* Plot *)
	Plot[{gLLfunction1[x], gLLfunction2[x], gLLfunction3[x]}, {x, a, b}, Exclusions -> None, PlotLegends -> {0.01, 0.1, 1}]];

PlotDOS[-1.5, 1.5]
```

This gives the density of states for three different values of λ:

<p align="center">
	<img src="Resources/DOS.png">
</p>
