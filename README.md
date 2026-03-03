# SIS_ANUPAMA_KAR_202603xx.zip 
Usage of AI: Used AI to understand and learn how to set everything up on GitHub, to help in markdown MathJax formatting for equations, and to learn how to create .gif from MATLAB, to polish my base code, to define some of the functions, and to add some fancy items like indicator plots. Also used it to research and read about different integrators to enhance my 
knowledge base.  

# Table of Contents
- [1. Singularity for asymmetric Euler angle sets](#1-singularity-for-asymmetric-euler-angle-sets-page-117)
- [2. Ambiguity of Euler parameters (quaternions)](#2-ambiguity-of-euler-parameters-quaternion-page-127)
- [3. Classical Rodrigues Parameters (CRP)](#3-classical-rodrigues-parameters-page-135)
- [4. Comparison of different numerical integrators](#4-comparison-of-different-numerical-integrators)
- [5. 3D animation explanation](#5-3d-animation-explanation)
- [References](#references)


## 1.	Singularity for asymmetric Euler angle sets
### 1.1. Background and set-up of the EOM  
Euler angles describe the orientation of a rigid body with respect to a fixed coordinate system using three successive rotations. Let $(\phi, \theta, \psi)$ be the three Euler angles (often called roll, pitch, and yaw).

![Figure 1.1: Successive yaw, pitch, and roll rotations](Successive%20yaw,%20pitch,%20and%20rotations.png)

**Figure 1.1:** Successive yaw, pitch, and roll rotations. This figure illustrates the sequence of rotations about the principal axes, corresponding to the Euler angles $(\psi, \theta, \phi)$.

The transformation from the body frame to the inertial frame can be represented by the rotation matrix $R$:

$$
R = M_3(\psi)\,M_2(\theta)\,M_1(\phi)
$$

$$
M_1(\phi) =
\begin{bmatrix}
1 & 0 & 0 \\
0 & \cos\phi & \sin\phi \\
0 & \sin\phi & \cos\phi
\end{bmatrix}
$$

$$
M_2(\theta) =
\begin{bmatrix}
\cos\theta & 0 & \sin\theta \\
0 & 1 & 0 \\
\sin\theta & 0 & \cos\theta
\end{bmatrix}
$$

$$
M_3(\psi) =
\begin{bmatrix}
\cos\psi & \sin\psi & 0 \\
\sin\psi & \cos\psi & 0 \\
0 & 0 & 1
\end{bmatrix}
$$
\
The body angular velocity vector in general form is:

$$
\omega_B = \begin{bmatrix} \omega_x \\ \omega_y \\ \omega_z \end{bmatrix}
$$

In the MATLAB script, the specific angular velocity vector is defined as a function of time $t$:

```matlab
w_fun = @(t) [0.2*sin(0.02*t);
				0.15;
				0.2*cos(0.02*t)];
```

Designate the Direction Cosine Matrix (DCM) relating Euler angle rates to body angular velocity as matrix $B$:

$$
B = \begin{bmatrix}
1 & 0 & -\sin\theta \\
0 & \cos\phi & \sin\phi \cos\theta \\
0 & -\sin\phi & \cos\phi \cos\theta
\end{bmatrix}
$$

Assign the Euler angle rates as:

$$
\dot{y} = \begin{bmatrix} \dot{\psi} \\ \dot{\theta} \\ \dot{\phi} \end{bmatrix}
$$
**Equation (3.57)** [2] shows the Euler Angles kinematic equation of motion as reproduced below:
<div style="border:2px solid black; padding:10px; display:inline-block;">
The Euler Angles kinematic equation of motion is:

$$
\dot{y} = B^{-1} \omega_B
$$
</div>

where 

$$
B^{-1} = \frac{1}{\cos \theta} 
\begin{bmatrix}
0 & \sin \phi & \cos \phi \\
0 & \cos \theta \cos \phi & -\cos \theta \sin \phi \\
\cos \theta & \sin \theta \sin \phi & \sin \theta \cos \phi
\end{bmatrix}
$$

**Note: Here B in the code is equal to $B^{-1}$ in the derivation.**
 
In the *main.m* script, it is coded as a function 

```matlab
function dydt = EulerODE(t,y,w_fun)
% y = [psi theta phi]'
psi = y(1);
theta = y(2);
phi = y(3);

w = w_fun(t);

s2 = sin(theta);
c2 = cos(theta);
s3 = sin(phi);
c3 = cos(phi);

B = 1/c2*[0 s3 c3;
          0 c2*c3 -c2*s3;
          c2 s2*s3 s2*c3];

dydt = B*w;

end
```
Looking at the inverse B matrix, the **singularity occurs when $\cos \theta = 0$, i.e., $\theta = \pi/2$ or 90°**. At this point, the factor $1/\cos\theta$ becomes **infinite**, making the matrix undefined. 

Physically, this corresponds to **gimbal lock**, where the Euler angles lose one degree of freedom: the rotation axes align such that you cannot uniquely determine all three angular rates from the Euler angle rates. Near $\theta = 90^\circ$, small changes in angular velocity cause large changes in the Euler rates, making control or simulation unstable.

### 1.2. How-To

### 1.3. Results

![Figure 1.2: 3-2-1 Euler Angles vs Time (ode45)](3-2-1%20Euler%20Angles%20ve%20Time%20(ode45).png)

**Figure 1.2:** 3-2-1 Euler Angles vs Time (ode45). This figure shows the time evolution of the Euler angles using the ode45 integrator.


### 1.4. Analysis
## 2.	Ambiguity of Euler parameters (quaternions)  
### 2.1. Background and set-up of the EOM 
To curb the singularity issue in Euler angles, four Euler parameters are used, mainly $\beta_0$, $\beta_1$, $\beta_2$, and $\beta_3$ (or $q_0$, $q_1$, $q_2$, and $q_3$).

<div style="border:2px solid black; padding:10px; display:inline-block;">
The Euler Parameters kinematic equation of motion is:

$$
\dot{\mathbf{q}} = \frac{1}{2} \, \Omega \, \mathbf{q}, \quad  
\Omega =
\begin{bmatrix}
0 & \omega_x & \omega_y & \omega_z \\
-\omega_x & 0 & -\omega_z & \omega_y \\
-\omega_y & \omega_z & 0 & -\omega_x \\
-\omega_z & -\omega_y & \omega_x & 0
\end{bmatrix}
$$
</div>  

Here, $\dot{\mathbf{q}}$ = $\dot{\boldsymbol{\beta}}$, $\mathbf{q}$ = $\boldsymbol{\omega}_B$, and $B(\boldsymbol{\beta})$ = $\Omega$.

The Junkins and Schaub textbook gives a thorough explanation and derivation of the quarternions. **Equation (3.104)** [2] is shown below:

$$
\dot{\boldsymbol{\beta}} =
\frac{1}{2}
\begin{bmatrix}
0 & -\omega_1 & -\omega_2 & -\omega_3 \\
\omega_1 & 0 & \omega_3 & -\omega_2 \\
\omega_2 & -\omega_3 & 0 & \omega_1 \\
\omega_3 & \omega_2 & -\omega_1 & 0
\end{bmatrix}
\begin{bmatrix}
\beta_0 \\ \beta_1 \\ \beta_2 \\ \beta_3
\end{bmatrix}
$$

where 
The quaternion column vector is:

$$
\boldsymbol{\beta} =
\begin{bmatrix}
\beta_0 \\
\beta_1 \\
\beta_2 \\
\beta_3
\end{bmatrix}
$$
or 
$$
\mathbf{q} =
\begin{bmatrix}
q_0 \\
q_1 \\
q_2 \\
q_3
\end{bmatrix}
$$
This is defined as a function on the *main.m* file as 

````matlab
function dqdt = QuatODE(t,q,w_fun)
% q = [b0 b1 b2 b3]'
w = w_fun(t);

Omega = [ 0    -w(1) -w(2) -w(3);
          w(1)  0     w(3) -w(2);
          w(2) -w(3)  0     w(1);
          w(3)  w(2) -w(1)  0];

dqdt = 0.5 * Omega * q;

end
````
Quaternions have a **sign ambiguity** because a quaternion and its negative represent the same physical rotation. That is, $\mathbf{q} = [q_0, q_1, q_2, q_3]^T$ and $-\mathbf{q} = [-q_0, -q_1, -q_2, -q_3]^T$ produce identical rotations in 3D space. This can cause apparent “jumps” in quaternion values when integrating or interpolating rotations if the sign flips, even though the actual orientation does not change. Care must be taken in numerical algorithms to maintain sign consistency to avoid discontinuities.

### 2.2. How-To  

![Figure 2.1: Quaternion Components](Quaternion%20Components.png)

**Figure 2.1:** Quaternion components. This figure shows the individual components of the quaternion over time.

![Figure 2.2: Quaternion Ambiguity Detection](Quaternion%20Ambiguity%20Detection.png)

**Figure 2.2:** Quaternion ambiguity detection. This figure illustrates the detection of ambiguity in quaternion representation.

![Figure 2.3: Quaternion Norm](Quaternion%20Norm.png)

**Figure 2.3:** Quaternion norm. This figure shows the norm of the quaternion, which should remain close to 1 for a valid rotation.
### 2.3. Results   
### 2.4. Analysis  
## 3.	Classical Rodrigues Parameters (CRP)  
### 3.1. Background and set-up of the EOM  
<div style="border:2px solid black; padding:10px; display:inline-block;">
The standard CRP kinematic equation is:

$$
\dot{\mathbf{r}} = \frac{1}{2} \left( \mathbf{I}_3 + [\mathbf{r}]_\times + \mathbf{r} \mathbf{r}^T \right) \boldsymbol{\omega}_B
$$
</div>

where

$$
[\mathbf{r}]_\times =
\begin{bmatrix}
0 & -r_3 & r_2 \\
r_3 & 0 & -r_1 \\
-r_2 & r_1 & 0
\end{bmatrix}.
$$

This matches the form given on **Equation (3.131)** [2]. 

On the MATLAB file, this CRP function is defined as below

````matlab
function drdt = CRPODE(t,r,w_fun)

w = w_fun(t);

r1 = r(1);
r2 = r(2);
r3 = r(3);

r_tilde = [  0   -r3   r2;
             r3    0   -r1;
            -r2   r1    0];

drdt = 0.5 * ( eye(3) + r_tilde + r*r' ) * w;

end
````

CRP was developed to mitigate the ambiguity issue arising from the quaternions. 
Classical Rodrigues Parameters (CRP) are defined from a quaternion 
$\mathbf{q} = [q_0, \mathbf{q}_v]^T$, where $\mathbf{q}_v = [q_1, q_2, q_3]^T$, as:

$$
\mathbf{r} = \frac{\mathbf{q}_v}{q_0} =
\frac{1}{q_0}
\begin{bmatrix}
q_1 \\ q_2 \\ q_3
\end{bmatrix}, \quad q_0 \neq 0 \ (\text{i.e., rotation angle } \theta \neq 180^\circ)
$$

The mapping from a rotation to CRP is **unique**, except at $\theta = 180^\circ$.

While this solves the ambiguity issue, it reintroduces the singularity issue at $\theta = 180^\circ$.

### 3.2. How-To  
### 3.3. Results  
### 3.4. Analysis  
## 4.	Comparison of different numerical integrators    
Singularity was hit. So, a tolerance was added for ode45 and ode15s.

![Figure 1.4: Integrator Step Size vs Time vs Pitch](Integrator%20Step%20Size%20vs%20Time%20vs%20Pitch.png)

**Figure 1.4:** Integrator step size vs time vs pitch. This figure illustrates how the integrator step size varies with time and pitch angle during the simulation.

## 5. 3D animation explanation  
Once MATLAB runs successfully, it will generate a .gif file and save it in the designated folder as *Aircraft_Attitude.gif*. The legends of the .gif are same as the legends of the Euler Angles plot (figure 1.2.). Since the .gif updates itiratively, putting a legend that would stay fixed was rendering bad graphics.  

![Figure 1.3: 3D Animation Frame](3D%20animation%20frame.png)

**Figure 1.3:** 3D animation frame. This figure shows a sample frame from the MATLAB-generated aircraft attitude animation.



It is to be noted that the animation is a rendering of the computation and in no way represents the actual physical motion of an aircraft. So, singularity in Euler angles or CRPs results in tumbling or rapid spinning. Quarternion ambiguity does not show any symptoms. Only singularity affects the actual computer attitude. 
## References
[1] Ross Dynamics Lab, “Euler Angle Simulation with MATLAB – Integrating the Rotational Kinematic Differential Equations,” YouTube video, 30 Jul. 2021. https://www.youtube.com/watch?v=vwn_JT0SDXQ. 

[2] Schaub, H., and Junkins, J. L., Analytical Mechanics of Space Systems, 2nd ed., American Institute of Aeronautics and Astronautics, Reston, VA, 2009.  

[3] MathWorks, “Choose an ODE Solver,” MATLAB & Simulink Documentation, The MathWorks, Inc., Natick, MA, accessed Mar. 2, 2026. Available: https://www.mathworks.com/help/matlab/math/choose-an-ode-solver.html