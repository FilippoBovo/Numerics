<!--- Load Prettify for Mathematica syntax highlight-->
<script src="https://cdn.rawgit.com/google/code-prettify/master/loader/run_prettify.js?lang=css&amp;skin=sunburst"></script>

# Lieb-Liniger Equations

This is a repository for numerical solutions of mathematical equations:


To solve the Fredholm integral equation of the second kind, we define the following three functions:

- **BoundFunction:** This module takes the function f as input and as output gives the function f in the interval [a,b] and zero otherwise.

	<pre class="prettyprint">
	BoundFunction[f_, a_, b_] :=
		Function[
			Piecewise[{
				{0., # < a},
				{0., # > b},
				{f[#], True}}]
		];
	</pre>

- **Fredholm2ndKind:** Gives the numerical solution of the Fredholm equation in the interval [a,b]. This takes as input the extremes of integration a and b, the kernel K(x,y), g(x) and the number, n, of subdivision of the integration interval which is used in the numerical solution. (see this link)

	<pre class="prettyprint">
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
	</pre>

- **Fredholm2ndKindOut:**  Gives the numerical solution of the Fredholm equation in the interval [c,d] outside of [a,b].This takes as input the extremes of integration a and b, the kernel K(x,y), g(x), the number, n, of subdivision of the integration interval [a,b] which is used in the numerical solution, the extreme of integration c and d of the interval outside [a,b] and the number of subdivisions of [c,d].

	<pre class="prettyprint">
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
	</pre>


