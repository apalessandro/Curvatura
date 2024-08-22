BeginPackage["Curvatura`"]

Connection::usage = "Connection[g,x] returns the Christoffel symbols for the metric g written in coordinates x";
RiemannTensor::usage = "RiemannTensor[g,x] returns the Riemann Tensor for the metric g written in coordinates x";
RicciTensor::usage = "RicciTensor[g,x] returns the Ricci Tensor for the metric g written in coordinates x";
RicciScalar::usage = "RicciScalar[g,x] returns the Ricci Scalar for the metric g written in coordinates x";
EinsteinTensor::usage = "EinsteinTensor[g,x] returns the Einstein Tensor for the metric g written in coordinates x";
WeylTensor::usage = "WeylTensor[g,x] returns the Weyl curvature Tensor for the metric g written in coordinates x";
KretschmannScalar::usage = "KretschmannScalar[g,x] returns the Kretschmann Scalar for the metric g written in coordinates x";
Geodesic::usage = "Geodesic[s,g,x] returns the geodesic equation for a test mass of proper time s in a curved spacetime described by the metric g written in coordinates x";
Singularity::usage = "Singularity[g,x] returns the singular points of the spacetime described by the metric g written in coordinates x";

Begin["`Private`"]

gg[g_, x_] := TensorProduct[Inverse[g], Grad[g, x]];

Connection[g_,x_] := 1/2(TensorContract[TensorTranspose[gg[g,x],{1,2,3,4,5}],{2,4}]+
TensorContract[TensorTranspose[gg[g,x],{1,2,5,4,3}],{2,4}]-
TensorContract[TensorTranspose[gg[g,x],{1,2,3,5,4}],{2,4}])//FullSimplify;

DC[g_,x_] := Grad[Connection[g,x],x];

CC[g_,x_] := TensorProduct[Connection[g,x],Connection[g,x]];

RiemannTensor[g_,x_] := (TensorTranspose[DC[g,x],{1,4,2,3}] - TensorTranspose[DC[g,x],{1,3,2,4}] + 
TensorContract[TensorTranspose[CC[g,x],{1,3,5,6,4,2}],{5,6}] - TensorContract[TensorTranspose[CC[g,x],{1,4,5,6,3,2}],{5,6}])//FullSimplify;

RicciTensor[g_,x_] := (TensorContract[TensorTranspose[DC[g,x],{1,2,3,4}],{1,4}]
- TensorContract[TensorTranspose[DC[g,x],{1,2,4,3}],{1,4}]
+ TensorContract[TensorTranspose[CC[g,x],{1,2,3,4,5,6}],{{1,5},{4,6}}]
- TensorContract[TensorTranspose[CC[g,x],{1,2,6,4,3,5}],{{1,5},{4,6}}])//FullSimplify;

RicciScalar[g_,x_] := TensorContract[TensorProduct[Inverse[g], RicciTensor[g,x]],{{1,3},{2,4}}]//FullSimplify;

EinsteinTensor[g_,x_] := (RicciTensor[g,x] - 1/2 g RicciScalar[g,x])//FullSimplify;

WeylTensor[g_,x_] := (TensorContract[TensorProduct[g,RiemannTensor[g,x]],{1,3}]
+1/2*(TensorTranspose[TensorProduct[RicciTensor[g,x],g],{1,4,2,3}] - TensorTranspose[TensorProduct[RicciTensor[g,x],g],{1,3,2,4}] + TensorTranspose[TensorProduct[RicciTensor[g,x],g],{2,3,1,4}] -TensorTranspose[TensorProduct[RicciTensor[g,x],g],{2,4,1,3}])
+1/6*RicciScalar[g,x]*(TensorTranspose[TensorProduct[g,g],{1,3,2,4}] - TensorTranspose[TensorProduct[g,g],{1,4,2,3}]))//FullSimplify;

KretschmannScalar[g_,x_] := TensorContract[
TensorProduct[
TensorContract[TensorProduct[g,RiemannTensor[g,x]],{1,3}],
TensorContract[TensorProduct[RiemannTensor[g,x],Inverse[g],Inverse[g],Inverse[g]],{{2,5},{3,7},{4,9}}]
]
,{{1,5},{2,6},{3,7},{4,8}}]//FullSimplify;

Geodesic[s_,g_,x_] := With[
{y := {Part[x,1][s], Part[x,2][s], Part[x,3][s], Part[x,4][s]}}, 
D[D[y,s],s] + TensorContract[TensorProduct[Connection[g,x], D[y,s], D[y,s]],{{2,4},{3,5}}]]//FullSimplify;

Singularity[g_,x_] := Quiet[
Check[
Reduce[Denominator[KretschmannScalar[g,x]]==0, x, Reals], 
Reduce[Denominator[KretschmannScalar[g,x]]==0, Reals]
]];

End[]
EndPackage[]
