# Curvatura
A Tensor Calculus Mathematica repository for calculations in General Relativity.

As an example, let's try to calculate the Ricci and Einstein tensors in an expanding universe.

In order to compute the curvature tensors, we need to define a metric and a set of coordinates. 

An isotropic and homogeneous universe is described by the Friedmann-Robertson-Walker metric. In spherical coordinates, this is
```math
ds^2 = -dt^2 + a^2(t)\left[\frac{dr^2}{1-kr^2} + r^2 (d\theta^2 + \sin^2\theta d\phi^2) \right]
```

We first load the Einstein Package:
```Mathematica
Needs["Einstein`"]
```

Then, we define the FRW metric 
```Mathematica
g = {
{-1,0,0,0},
{0,a[t]^2/(1-k r^2),0,0},
{0,0,a[t]^2r^2,0},
{0,0,0,a[t]^2r^2 Sin[\[Theta]]^2}
}
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/400158ae-f3e3-4ba4-b577-bb6394a6d2bc)

and the spherical coordinates
```Mathematica
x = {t,r,\[Theta],\[Phi]}
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/b27f0e21-9933-4504-8b10-f489e1ba14f4)

We can now compute the connection coefficients by using the pre-defined function Connection[g,x]:
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
