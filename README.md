# Curvatura
A Tensor Calculus Mathematica package for calculations in General Relativity.

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

Then, we define the metric 
```Mathematica
g = {
{-1,0,0,0},
{0,a[t]^2/(1-k r^2),0,0},
{0,0,a[t]^2r^2,0},
{0,0,0,a[t]^2r^2 Sin[\[Theta]]^2}
}
```
and the spherical coordinates
```Mathematica
x = {t,r,\[Theta],\[Phi]}
```

We can now compute the connection coefficients by using the pre-defined function Connection[g,x]:
```Mathematica
Connection[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/9d5df438-d71d-4874-8d0e-0ca5240a907f)

We do the same for the Ricci tensor
```Mathematica
RicciTensor[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/006ac216-a150-49ed-986b-4961d94dc4a9)
the Ricci scalar
```Mathematica
RicciScalar[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/aeee824a-f141-4ecb-8f0c-600169e42486)

and the Einstein tensor
```Mathematica
EinsteinTensor[g,x]
```
![image](https://github.com/apalessandro/Curvatura/assets/48097299/04595de1-308e-4c5c-880a-9c23789a1722)



