A few years later (in fact in 2002, so 17 years after Perlin published his first paper on noise), Perlin proposed a new improved version of his noise function. The original noise function has two problems.

## The smoothstep function is not ideal and here is why...

First, remember that the interpolant we used to interpolate the result of the dot product between the gradients and the vectors pointing toward the point within the cell, is the smoothstep function. This function has the form \(3t^2 - 2t^3\). Remember as well that we used the derivative of this function (\(6t - 6t^2\)) to compute the Noise function partial derivatives. As explained in the previous chapter, these derivatives are useful to compute the "true" normal (true as in mathematically correct) of vertices displaced by a noise function. Now the problem with this is that the derivative of the smoothstep function derivative (what we call in mathematics the second order derivative or \(f''(x)\)) is not continuous when \(t=0\) or \(t=1\). The second order derivative of the smoothstep function is \(6 - 12t\). So when \(t=0\) \(f''(x)=6\) and when \(t=1\), \(f''(x)=-6\). This is a problem because if the derivative of the function used to compute the normal of the displaced mesh is not continuous then it will introduce a discontinuity in these normals wherever \(t=0\) and \(t=1\). That happens when the vertices of the displaced mesh are aligned with the edges of the noise grid cells (i.e. when the vertices coordinates have integer values). These visual artefacts are clearly visible in the image below (left).

It is important to understand that while the analytical solution to compute the noise function derivative is correct mathematically, this doesn't mean that the result is necessarily visually what we look for. This is nothing else than a mathematical construct and this construct may have properties we don't like. So while the results are mathematically correct, they create visual artefacts we would rather prefer to avoid.

![](/images/noise-part-2/noise-discontinuity2.png?)

Furthermore most documents mentioning the problem of this function don't really explain why we have this problem and provide a way of visualising it. To better understand why the problem is there in the first place, you need to remember a detail about functions and function derivatives. In physics, a common function is the function that allows you to compute the **position** of a moving object as a function of a time. Let's call this function \(p(t)\). If you recall your lessons on physics, you can compute the **speed** of this object as the derivative of this function with respect to time, in other words \(p'(t)=p(t)/\partial t\). And finally you can compute the object **acceleration** as the derivative of the speed function with respect to time again or \(p''(t)=p'(t)/\partial t\). So in essence you can see a function \(f(x)\) and its derivatives \(f'(x)\) and \(f''(x)\) as respectively the position, the speed and the acceleration of an object. In our example the smoothstep function is by analogy, the position, the first order derivative of the smoothstep function that we used to compute the tangent is the speed and finally the second order derivative of the smoothstep function is the acceleration. Of course these are just analogies. But the point is, if the acceleration has discontinuities, then we can expect the speed to change abruptly from time to time. In our particular case, this means that if the smoothstep function second order derivative has discontinuities, then we can expect the tangent to vary abruptly from time to time. And that's precisely what causes the visual artefacts at the cell boundaries. The normals in these regions change abruptly.

<details>
In computer graphics and mathematics, when the original function is continuous we speak of geometry continuity and we say the function is C0 (or G0 when we deal with geometry). If the function derivative is continuous we speak of tangent continuity (C1 or G1). And when the function second order derivative is also continuous we speak of curvature continuity, in other word, the curvature where two curves or two surfaces join is the same for the curves or two surfaces at the point of contact).
</details>

The image below shows a slice of the noise function and its analytically computed normals. It is pretty clear that the orientation of the normals and the direction of these normals is not the same on each side of the cell boundary.

![](/images/noise-part-2/noise-discontinuity.png?)

![Figure 1: bicubic (smoothstep) and quintic interpolant.](/images/noise-part-2/noise-curves.png?)

Hopefully there is a solution. Well, if the only problem is that the second order derivative of the smoothstep function is discontinuous then all we need to do to fix the problem, is to choose a function for the interpolant that has a second order derivative that is, obviously, continuous. There are quite a few of them that have a similar shape to the smoothstep function but Perlin choose this quintic function (figure 1):

$$6t^5-15t^4+10t^3.$$

Its first order derivative is:

$$30t^4-60t^3+30t^2.$$

And its second order derivative is:

$$120t^3-180t^2+60t.$$

As you can see this function is continuous (when t=0 and t=1 the function is equal to 0).

![](/images/noise-part-2/noise-discontinuity1.png?)

## Choosing random gradients is not ideal and here is why...

The next problem is that because the gradients are chosen randomly, sometimes the gradient directions are aligned along the same axis or along the diagonal of the cell as showed in figure 2 [missing]. The noise function tends to take on very high values in these regions (close to 1) causing in parts the "splotchy appearance" of the original Noise implementation.

Hopefully this problem can be solved easily by replacing the 256 random directions by a set of 12 directions that are different from each other (and excluding situations in which these gradients would point along the same direction). The permutation table provides plenty of randomness so in practice it doesn't matter if these 12 vectors are predefined. Perlin in his 2002 paper decided to choose the following vectors:

$$
\begin{array}{l}
(1,1,0),(-1,1,0),(1,-1,0),(-1,-1,0),\\
(1,0,1),(-1,0,1),(1,0,-1),(-1,0,-1),\\
(0,1,1),(0,-1,1),(0,1,-1),(0,-1,-1)
\end{array}
$$

![Figure 3: the 12 directions chosen by Perlin as gradients.](/images/noise-part-2/noise-12directions.png?)

Though keep in mind that our permutation table returns integers in the range [0,255], so in order to choose one of these vectors randomly we would need to do something like <span class="code-inline">P[Hash(x,y,z)] % 12</span> though as you know we don't like to use the modulo operator too much because it's slower than bitwise operators and it requires to handle negative values differently than positive values (bitwise operators fix this problem). To replace the modulo operator with a bitwise operator, Perlin proposes to extend the array of 12 vectors to 16 vectors, adding to the first 12 the following 4 directions:

$$
\begin{array}{l}
(1,1,0),(-1,1,0),(0,-1,1), (0,-1,-1).
\end{array}
$$ 

We already used these directions but Perlin claims that because the first 12 directions form a regular tetrahedron, adding some of these directions redundantly doesn't introduce any bias in the texture.

> The final result has the same non-directional appearance as the original distribution but less clumping.

The splotchy like appearance is not as visible anymore (check out the comparaison image below). Implementing this change is actually quite simple. Our hash function returns numbers in the range [0,255] so we need to use another bitwise operator to set this number to the range [0,15].

```
uint8_t p = hash(xi, yi, zi) & 15;
```

Then note that because all directions are defined with 0 and 1, the dot product simplified into a sum of the coordinates of the right-hand vector.

```
// float a = dot(c000, p000);
// float a = c000.x * p000.x + c000.y * p000.y + c000.z * p000.z;
// for the gradient (1, 1, 0) for instance, this simplifies to:
float a = 1 * p000.x + 1 * p000.y + 0 * p000.z = p000.x + p000.y; 
```

We can move this code into its own function:

```
float gradientDotV( 
    // a value between 0 and 255 
    uint8_t perm, 
    // vector from one of the corners of the cell to the point where the noise function is computed
    float x, float y, float z) 
{ 
    switch(perm & 15) { 
        case  0: return  x + y;  //(1,1,0) 
        case  1: return -x + y;  //(-1,1,0) 
        case  2: return  x - y;  //(1,-1,0) 
        case  3: return -x - y;  //(-1,-1,0) 
        case  4: return  x + z;  //(1,0,1) 
        case  5: return -x + z;  //(-1,0,1) 
        case  6: return  x - z;  //(1,0,-1) 
        case  7: return -x - z;  //(-1,0,-1) 
        case  8: return  y + z;  //(0,1,1), 
        case  9: return -y + z;  //(0,-1,1), 
        case 10: return  y - z;  //(0,1,-1), 
        case 11: return -y - z;  //(0,-1,-1) 
        case 12: return  y + x;  //(1,1,0) 
        case 13: return -x + y;  //(-1,1,0) 
        case 14: return -y + z;  //(0,-1,1) 
        case 15: return -y - z;  //(0,-1,-1) 
    } 
} 
 
float eval(...) 
{ 
    ... 
    float a = gradientDotV(hash(xi0, yi0, zi0), x0, y0, z0); 
    float b = gradientDotV(hash(xi1, yi0, zi0), x1, y0, z0); 
    float c = gradientDotV(hash(xi0, yi1, zi0), x0, y1, z0); 
    float d = gradientDotV(hash(xi1, yi1, zi0), x1, y1, z0); 
    float e = gradientDotV(hash(xi0, yi0, zi1), x0, y0, z1); 
    float f = gradientDotV(hash(xi1, yi0, zi1), x1, y0, z1); 
    float g = gradientDotV(hash(xi0, yi1, zi1), x0, y1, z1); 
    float h = gradientDotV(hash(xi1, yi1, zi1), x1, y1, z1); 
    ... 
}
```

The image below shows a comparison between the original Perlin noise implementation (left) and the improved version (right):

![](/images/noise-part-2/noise-compare.png?)

While not being significantly different the two implementations give distinct results. The improved version creates vertical and horizontal lines though that you don't see in the orignal version. This due to the fact that in the improved version, the gradients are aligned along the xy, xz and yz planes. While the choice of a quintic function as interpolant definitely reduces the visual artefacts and can be called an improvement, we don't think using the predefined gradients over random gradients is necessarily a good thing. You can find the complete implementation of the improved version in the last chapter of this lesson.

```
inline 
float quintic(const float &t) 
{ 
    return t * t * t * (t * (t * 6 - 15) + 10); 
} 
 
inline 
float quinticDeriv(const float &t) 
{ 
    return 30 * t * t * (t * (t - 2) + 1); 
} 
 
class PerlinNoise 
{ 
    ... 
    /* inline */ 
    uint8_t hash(const int &x, const int &y, const int &z) const 
    { return permutationTable[permutationTable[permutationTable[x] + y] + z]; } 
 
    float gradientDotV( 
        // a value between 0 and 255
        uint8_t perm, 
        // vector from one of the corners of the cell to the point where the noise function is computed
        float x, float y, float z) const 
    { 
        switch (perm & 15) { 
            case  0: return  x + y;  //(1,1,0) 
            case  1: return -x + y;  //(-1,1,0) 
            case  2: return  x - y;  //(1,-1,0) 
            case  3: return -x - y;  //(-1,-1,0) 
            case  4: return  x + z;  //(1,0,1) 
            case  5: return -x + z;  //(-1,0,1) 
            case  6: return  x - z;  //(1,0,-1) 
            case  7: return -x - z;  //(-1,0,-1) 
            case  8: return  y + z;  //(0,1,1), 
            case  9: return -y + z;  //(0,-1,1), 
            case 10: return  y - z;  //(0,1,-1), 
            case 11: return -y - z;  //(0,-1,-1) 
            case 12: return  y + x;  //(1,1,0) 
            case 13: return -x + y;  //(-1,1,0) 
            case 14: return -y + z;  //(0,-1,1) 
            case 15: return -y - z;  //(0,-1,-1) 
        } 
    } 
 
    float eval(const Vec3f &p, Vec3f& derivs) const 
    { 
        int xi0 = ((int)std::floor(p.x)) & tableSizeMask; 
        int yi0 = ((int)std::floor(p.y)) & tableSizeMask; 
        int zi0 = ((int)std::floor(p.z)) & tableSizeMask; 
 
        int xi1 = (xi0 + 1) & tableSizeMask; 
        int yi1 = (yi0 + 1) & tableSizeMask; 
        int zi1 = (zi0 + 1) & tableSizeMask; 
 
        float tx = p.x - ((int)std::floor(p.x)); 
        float ty = p.y - ((int)std::floor(p.y)); 
        float tz = p.z - ((int)std::floor(p.z)); 
 
        float u = quintic(tx); 
        float v = quintic(ty); 
        float w = quintic(tz); 
 
        // generate vectors going from the grid points to p
        float x0 = tx, x1 = tx - 1; 
        float y0 = ty, y1 = ty - 1; 
        float z0 = tz, z1 = tz - 1; 
 
        float a = gradientDotV(hash(xi0, yi0, zi0), x0, y0, z0); 
        float b = gradientDotV(hash(xi1, yi0, zi0), x1, y0, z0); 
        float c = gradientDotV(hash(xi0, yi1, zi0), x0, y1, z0); 
        float d = gradientDotV(hash(xi1, yi1, zi0), x1, y1, z0); 
        float e = gradientDotV(hash(xi0, yi0, zi1), x0, y0, z1); 
        float f = gradientDotV(hash(xi1, yi0, zi1), x1, y0, z1); 
        float g = gradientDotV(hash(xi0, yi1, zi1), x0, y1, z1); 
        float h = gradientDotV(hash(xi1, yi1, zi1), x1, y1, z1); 
 
        float du = quinticDeriv(tx); 
        float dv = quinticDeriv(ty); 
        float dw = quinticDeriv(tz); 
 
        float k0 = a; 
        float k1 = (b - a); 
        float k2 = (c - a); 
        float k3 = (e - a); 
        float k4 = (a + d - b - c); 
        float k5 = (a + f - b - e); 
        float k6 = (a + g - c - e); 
        float k7 = (b + c + e + h - a - d - f - g); 
 
        derivs.x = du *(k1 + k4 * v + k5 * w + k7 * v * w); 
        derivs.y = dv *(k2 + k4 * u + k6 * w + k7 * v * w); 
        derivs.z = dw *(k3 + k5 * u + k6 * v + k7 * v * w); 
 
        return k0 + k1 * u + k2 * v + k3 * w + k4 * u * v + k5 * u * w + k6 * v * w + k7 * u * v * w; 
    } 
    ... 
};
```

## What's Next?

This concludes our lesson on Perlin noise. The lesson is quite extensive already but there is certainly a lot more to say about noise functions in general. First, there exists quite a few more ways of generating a procedural noise. Some other common methods are the **wavelet noise** (developed by Pixar) the **Gabor noise** and the **simplex noise** which is also a kind of noise developed by Ken Perlin. We will write lessons on each one of these techniques in the future.

Finally it is important to be able to judge the quality of noise implementation against others. How do we that? We already mentioned the qualities that a noise function should have although. Among these qualities the noise should have a distribution of frequencies as uniform as possible and yet still appear random and smooth. One way of finding out whether a given noise function fits this criteria is to analyse its output in frequency domain. We will also learn about this technique in a future lesson.

## How to Use the Noise Function?

A procedural noise function in general can be used for many things:

- Terrain: in this lesson, we just studied the implementation of the function itself, though [we also showed how to use it to create something like a terrain](/lessons/procedural-generation-vritual-worlds%20/perlin-noise-part-2/perlin-noise-terrain-mesh). We haven't shown yet how to create terrains with really complex and realistic surface details using the noise function. We will show this technique in a separate lesson yet to be published.
- Water surfaces: it can also be used to simulate a water surface. By offsetting the the y coordinative of the input point, it is possible to animate the displaced surface like waves. This is simple technique that works well though there are better methods to simulate water surfaces.
- Cloud: it can also be used to add surface details to volumes such as cloud. We use the noise function in the lesson on volume rendering to simulate the appearance of various cloudy shapes.
- Texturing: finally it can be used of couse to add texture details to surfaces. It can be used to modulate any input of the shader such as the color, the specular, or add some bump to the surface of an object.
- Animation: the noise function can also be used to add some noise to animations. A typical use case for example is to add some noise to a camera animation to simulate some camera shakes.

## References

_Improving Noise_. Ken Perlin (2002)

_State of the Art in Procedural Noise Functions_. A. Lagae and al. (2010)