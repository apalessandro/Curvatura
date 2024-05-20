BeginPackage["Einstein`"]

Connection::usage = "Connection[g,x] returns the Christoffel symbols corresponding to the metric g and coordinates x";
RiemannTensor::usage = "RiemannTensor[g,x] returns the Riemann Tensor corresponding to the metric g and coordinates x";
RicciTensor::usage = "RicciTensor[g,x] returns the Ricci Tensor corresponding to the metric g and coordinates x";
RicciScalar::usage = "RicciScalar[g,x] returns the Ricci Scalar corresponding to the metric g and coordinates x";
EinsteinTensor::usage = "EinsteinTensor[g,x] returns the Einstein Tensor corresponding to the metric g and coordinates x";

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

Ricci[g_,x_] := (TensorContract[TensorTranspose[DC[g,x],{1,2,3,4}],{1,4}]
- TensorContract[TensorTranspose[DC[g,x],{1,2,4,3}],{1,4}]
+ TensorContract[TensorTranspose[CC[g,x],{1,2,3,4,5,6}],{{1,5},{4,6}}]
- TensorContract[TensorTranspose[CC[g,x],{1,2,6,4,3,5}],{{1,5},{4,6}}])//FullSimplify;

RicciScalar[g_,x_] := TensorContract[TensorProduct[Inverse[g], Ricci[g,x]],{{1,3},{2,4}}]//FullSimplify;

R[g_, x_] := TensorContract[TensorProduct[Inverse[g], Ricci[g,x]],{{1,3},{2,4}}]//FullSimplify;

EinsteinTensor[g_,x_] := (Ricci[g,x] - 1/2 g R[g,x])//FullSimplify;

End[]
EndPackage[]
