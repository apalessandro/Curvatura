# Curvatura
A Mathematica library for tensor calculus and dynamical evolution in General Relativity.

As an example, let's try to calculate the curvature tensors in an expanding universe.

In order to do that, we need to define a metric and a set of coordinates. 

An isotropic and homogeneous universe is described by the Friedmann-Robertson-Walker metric. In spherical coordinates, this is
```math
ds^2 = -dt^2 + a^2(t)\left[\frac{dr^2}{1-kr^2} + r^2 (d\theta^2 + \sin^2\theta d\phi^2) \right]
```

We first load the Curvatura Package:
```Mathematica
Needs["Curvatura`"]
```

Then, we define the FRW metric as a (4,4) tensor
```Mathematica
g = {
{-1,0,0,0},
{0,a[t]^2/(1-k r^2),0,0},
{0,0,a[t]^2r^2,0},
{0,0,0,a[t]^2r^2 Sin[\[Theta]]^2}
}
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/400158ae-f3e3-4ba4-b577-bb6394a6d2bc)

and the spherical coordinates the metric is written in
```Mathematica
x = {t,r,\[Theta],\[Phi]}
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/b27f0e21-9933-4504-8b10-f489e1ba14f4)

We can now compute the connection coefficients (Christoffel symbols) using the pre-defined function Connection[g,x]:
```Mathematica
Connection[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/3072f400-abe5-4821-ba7d-cac63f260097)

We do the same for the Ricci tensor
```Mathematica
RicciTensor[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/a97b734b-cd5f-4b56-802a-99707542add1)

the Ricci scalar
```Mathematica
RicciScalar[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/bd9325c4-40df-49e1-8cc3-1b6a4202ea28)

and the Einstein tensor
```Mathematica
EinsteinTensor[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/50387ad2-60cc-464c-acc1-d73f78d52dc3)

We can also check explicitly that the contraction of the Riemann curvature tensor gives the Ricci tensor ($R_{ab} = R^c_{acb}$):
```Mathematica
TensorContract[RiemannTensor[g,x],{1,3}]==RicciTensor[g,x]//FullSimplify
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/87a7a13f-dda7-42ab-ae1e-f890f792d173)

The Weyl curvature tensor vanishes identically in any conformally flat metric, and we can check this explicitly for FRW by using the relevant function:
```Mathematica
WeylTensor[g,x] == ConstantArray[0,{4,4,4,4}]//FullSimplify
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/9ada3762-575e-4a9b-bc33-1b5e6327e9c4)

We can also calculate the Kretschmann scalar for the FRW metric:
```Mathematica
KretschmannScalar[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/0e2740ee-1877-4a30-95e7-89e1d6450d9c)

The Singularity function returns all the singular points of the spacetime described by the metric $g$. In the FRW spacetime, there is one gravitational singularity at $a(t) = 0$, corresponding to the Big Bang:
```Mathematica
Singularity[g,x]
```
![image](https://github.com/user-attachments/assets/042c67fd-9adb-487c-9659-a5df83f30876)

