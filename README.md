[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://github.com/dindagustiayu/First-Order-Differential-Equations-in-Chemistry/blob/main/First%20ODE%20in%20Chemistry.py)


# Ordinary Differential Equation of Chemistry

A differential equation (DE) is an equation that relates some derivatives of an unkniwn function to each other and to this function itself. An ordinary differential equation is a (ODE), where the unknown function depends on one variable only. We call this variable the independent variable and the unknown function the dependent variable. If the unknown function depends on several independent variables and the equation involves partial derivatives, then we have a partial differential equation. A chemicallly reacting system without space variations can be described by a set of ordinary with time $t$ as independent variable.

The mathematical description of various processes in chemistry and physics is possible by describing then with the help of differential equations which are based on simple model assumptions and defining the boundary conditions.

## First-order Differential Equations
Many processes and phenomena in chemistry, and generally in sciences, can be described by first-order differential equations. These equations are the most important and most frequently used to describe natural laws. 
The general form of a first order equation can be defined as:
<p align='center'>
    $$\frac{dy}{dt} = f(t, \ y) \quad(1)$$
</p>

where $f(t, \ y)$ is some functions of the independent variable (e.g., time $t$) and $y$ itself, the solution is sought, $y(t)$, in general, varies with time in some way determined by an _initial condition_: a known value of $y$ at some particular time (e.g., $t =0$).

## Chemical reaction kinetics

Chemical reaction kinetics is the study of rates of chemical processes (reactions). The goal is to find the relations between the concetrations of educts or products of a chemical reaction (as depending variable) and the time (as independent variable). In general, _all_ chemical reactions can be described mathematically by first-order differential equations. Their solutions, however, depend directly on the nature of the chemical reaction itself. The latter is characterized  by the so-called _reaction order_, which has nothing what so ever to do with the _order of a differential equation_. The reaction order of a chemical reaction is simply defined by the sum of exponents of concentrations occuring in th rate law.

Simple reactions like the transformation of $A \rightarrow P$ can be described by the differential equation:
<p align='center'>
    $$-d[A] = k \cdot A \cdot dt \quad(2)$$
</p>

The first order kinetics are summarized by the differential equation:
<p align='center'>
    $$\frac{d[A]}{dt} = -k [A]$$
</p>

where $[A]_0$ is the initial concentration of $[A]$.


The _first-order differential equation_ describes also a _first-order reaction_ in chemical kinetics, due to the exponent 1 of the concentration $A$. Following a separation of variables, the integration results in:

<p align='center'>
    $$ \int_{A_0}^{A_t} \frac{d[A]}{[A]} = -k \int_0^{t} \ dt \quad(3)$$
</p>
The result

<p align='center'>
    $$ \ln \frac{A_t}{A_0} = -k t \quad(4)$$
</p>

can be also written in the exponential form:
<p align='center'>
    $$[A](t) = [A]_0 e^{-k t} \quad(5)$$
</p>


## Preliminaries

Ordinary differential equations can be solved numerically with,
- `scipy.integrate.solve_ivp`: solve an initial value problem.
- `scipy.integrate.odeint`: solve first-order differential equations.

This function is based on the well-tested Fortran LSODA routine, which can automatically switch between stiff anf nostiff algorithms.

`solve_ivp` or `odeint` take three arguments:
- A function object returning $dy/dt$ as $f(t, \ y)$, given $t$ and $y$.
- An initial condition, $y_0$.
- A sequence of $t$, the initial and final time values at which to calculate the solution, $y(t)$, as a tuple ($t0, \ tf$).

The derivative function is simply:
```
def dydt(y, t):
    return -k * y
```
the order of the argument is important. Or 

```
def dydt(y, t):
    return -k * y

soln = solve_ivp(dydt, (t0, tf), [y0])
```
The returned object, `soln`, contains information about the solution.

### Example - The decomposition of $N_2O_5$
The decomposition of $N_2O_5$ is first order with a rate constant, $k = 5 \times 10^{-4} \ s^{-1}$. Given an initial consentration of $[N_2O_5] = 0.02 \ mol \ dm^{-3}$, what is the concentration of $N_2O_5$ after 2 hours?.

### Solution
1. Write a single dependent variable, $y(t) \equiv [A]$, which is a function of the independent variable, $t$ (time). We have:
<p align='center'>
    $$ \frac{dy}{dt} = -ky$$
</p>

2. We need to provide a function returning $dy/dt$ as $f(y, \ t)$ (in general a function of both $y$ and $t$), an intitial condition, $y(0)$ and a sequence of time points upon which to calculate the solution.
   
3. Calculate concentration through analytical equation, $\[ A \](t)= [A]_0 \ e^{-kt}$.

4. Determinate the value of $y(t_f)$ by the ODE solver algorithm.

5. Set the `dense_Output` argument to `True` to define an `Odesolution`.


```Python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# First-order reaction rate constant, s-1
k = 5.0e-4

# Initial condition on y: 100% of N2O5 is present at t = 0
y0 = 0.02

# Initial and final time points for the integration (s).
t0, tf = 0, 2*60 * 60

# 1. Return dy/dt = f(y, t) at t time
def dydt(t, y):
    return -k * y

# 2. Integrate the differential equation
soln = solve_ivp(dydt, (t0, tf), [y0])
print(soln)
```
```
  message: The solver successfully reached the end of the integration interval.
  success: True
   status: 0
        t: [ 0.000e+00  4.618e-01  5.080e+00  5.126e+01  5.130e+02
             2.376e+03  4.149e+03  5.960e+03  7.200e+03]
        y: [[ 2.000e-02  2.000e-02  1.995e-02  1.949e-02  1.547e-02
              6.101e-03  2.515e-03  1.018e-03  5.475e-04]]
      sol: None
 t_events: None
 y_events: None
     nfev: 50
     njev: 0
      nlu: 0
```

The solution was obtained succesfully in 50 function evaluations and the initial, intermediate and final values of and are reported. Next, we need calculate consentration of  after 2 hours.

```Python
y0 = y0 * np.exp(-k * tf)
print(f'N205:{y0}')
```
```
N205:0.0005464744489458512
```
The concentration of $N_2O_5$ after 7200 secs is $5.46 \times 10^{-4} \ mol \ dm^{-3}$. 

The twenty-one time points, given as the array $t$ in the returned solution, that were used to determine the value of $y(t_f)$ were chosen by the ODE solver algorithm, To follow the reactant concentration as a function of time in higher resolution, provide a suitable sequence of time points as the argument `t_eval`.
```Python
# 4. Integrate the differential equation, report at a 21 time points.
t_eval = np.linspace(tf, t0, 21)
soln = solve_ivp(dydt, (tf, t0), [y0], t_eval = t_eval)

t, y = soln.t, soln.y[0]

print(f't: {t}')
print(f'y: {y}')
```
```

t: [7200. 6840. 6480. 6120. 5760. 5400. 5040. 4680. 4320. 3960. 3600. 3240.
 2880. 2520. 2160. 1800. 1440. 1080.  720.  360.    0.]
y: [0.00054647 0.00065425 0.00078329 0.00093786 0.00112272 0.00134387
 0.00160875 0.00192618 0.00230643 0.00276162 0.00330668 0.00395789
 0.00473649 0.00566897 0.00678709 0.00812786 0.00973354 0.01165358
 0.01395182 0.01670321 0.01999755]
```
```Python
# 5. Solve the ODE
soln = solve_ivp(dydt, (t0, tf), [y0], dense_output=True)

# Evaluate the solution, y(t), at 21 time points
t= np.linspace(t0, tf, 21)
y = soln.sol(t)[0]

# Plot and compare the numerical and extract solutions
plt.plot(t, y, 'o', color='k', label=r'solve_ivp')
plt.plot(t, y0 * np.exp(-k * t), color='gray', label='exact')
plt.xlabel('t/s')
plt.ylabel(r'$[N_2 O_5]$(t)  / mol.dm-3')
plt.legend
plt.savefig('Numerical and Analytical ODE.svg', bbox_inches='tight')
plt.show()
```
<div align='center'>
<img src='https://github.com/dindagustiayu/First-Order-Differential-Equations-in-Chemistry/blob/main/First%20ODE%20Kinetic%20reactions/Numerical%20and%20Analytical%20ODE.svg'>
</div>

<p align='center'>
Figure 1. Numerical and analytical solutions to the ordinary differential equation governing the decomposition of $[N_2 O_5]$.
</p>

The resulting plot demonstrates that the numerical algorithm was able to follow the true solution accurately in this case.

## Coupled First-Order Differential Equations

In chemistry we are usually concerned with systems of DE of the first order in more than one dependent variable: $y_1(t), \ y_2(t), \ldots, \ y_n(t)$:

<p align='center'>
    $$\begin{align} \frac{dy_1}{ft} &=f_1 (y_1, y_2, \ldots, y_n;t), \\ \frac{dy_2}{dt} &= f_2(y_1, y_2, \ldots, y_n;t), \\ \vdots \\ \frac{dy_n}{dt} &= f_n (y_1, y_2, \ldots, y_n;t) \end{align}$$
</p>

In this case, the function passed to `solve_ivp` must return a sequence of derivatives, $dy_1/dt, \ dy_2/dt, \ldots, \ dy_n/dt$ for each of the dependent variables, that is, it evaluates the earlier mentioned functions $f_i(y_1, y_2, \ldots, y_n;t)$ for each of the y_i passed to it in a sequence, y. 

In majority of chemical reactions, however, more than one educt is involve. A reaction proceeds via two first-order reaction steps:
<p align='center'>
    $$ A \longrightarrow B \longrightarrow P $$
</p>

with rate constants  $k_1$ and $k_2$. 

<p align='center'>
    $$ \begin{align} A \rightarrow B \quad k_1 \\ B \rightarrow P \quad \end{align}$$
</p>

The equations governing the rate of change of A and B are
<p align='center'>
    $$\begin{align} \frac{d[A]}{dt} &= -k_1 [A] \quad (6) \\ \frac{d[B]}{dt} &= k_1[A] -k_2 [B] \quad(7) \end{align}$$
</p>

We can solve this pair of coupled equations analytically, but in our numerical solution, let $y_1 \equiv [A]$ and $y_2 \equiv [B]$:
<p align='center'>
    $$\begin{align} \frac{dy_1}{dt} &= -k_1 y_1 \quad (8) \\ \frac{dy_2}{dt} &= k_1 y_1 - k_2 y_2 \quad(9) \end{align}$$ 
</p>

The derivative function would then be defined as:
```
def deriv(t, y):
    """return dy_i/dt for each y_i at time t,"""
    y1, y2 = y
    dy1dt = -k1 * y1
    dy2dt = k1 * y1 - k2* y2
    return dy1dt, dy2dt
```

The code mentioned here integrates equations above for $k_1 = 0.2 \ s^{-1}$, $k_2 = 0.8 \ s^{-1}$ and initial conditions $y_1(0) = 100, \ y_2(0) = 0 $, and compares with the analytical result (__Figure 1__).

```Python
# First-order reaction rate constant, s-1
k1, k2 = 0.2, 0.8

# Initial concentration of species A and B (mol. dm-3)
y0 = A0, B0 = 100, 0

# Initial and final time points for the integration (s)
t0, tf = 0, 40

def deriv(t, y, k1, k2):
    y1, y2 = y
    dy1dt = -k1 * y1
    dy2dt = k1 * y1 - k2 * y2
    return dy1dt, dy2dt

# Integrate the differential equation.
soln = solve_ivp(deriv, (t0, tf), y0, dense_output=True, args=(k1, k2))
print(soln.message)

t = np.linspace(t0, tf, 200)
A, B = soln.sol(t)

# The concentration, [P], is determined by conservation
P = A0- A -B

plt.plot(t, A, ls='-', label=r'$ \mathrm{[A]}$')
plt.plot(t, B, ls=':' , label=r'$ \mathrm{[B]}$')
plt.plot(t, P, ls = '-.', label=r'$ \mathrm{[P]}$')
plt.xlabel(r'$ \mathrm{t / s}$')
plt.ylabel(r'$ \mathrm{conc / mol \ dm^{-3}}$')
plt.legend()
plt.savefig('ODE in more than on dependent variables.svg', bbox_inches='tight')
plt.show()
```
The solver successfully reached the end of the integration interval.
<div align='center'>
<img src='https://github.com/dindagustiayu/First-Order-Differential-Equations-in-Chemistry/blob/main/First%20ODE%20Kinetic%20reactions/ODE%20in%20more%20than%20on%20dependent%20variables.svg'>
</div>

<p align='center'>
Figure 2. Numerical solution to the ODE governing the reaction the sequence $A \rightarrow B \rightarrow P$ for rate constant $k_1 = 0.2 \ s^{-1}, \ k_2 = 0.8 \ s^{-1}$ and taking $[A]_0 = 100 \ mol \ dm^{-3}, \ [B]_0 = [P]_0 = 0$.
</p>

## Example - The Chapman Mechanism
Ozone destruction and regeneration occur naturally in the stratosphere. A simple mechanism for the formation consists of following four reactions (known as the _[Chapman cycle](https://en.wikipedia.org/wiki/Herbert_Chapman)_ in 1930):

<div align='center'>
    $$\mbox{Chapman Mechanism of Ozone Production and Destriction}$$
</div>

<div align='center'>

|Reaction | Rate| Kinetics order |
|--------------------------|--------------------------|----------|
|1) $O_2 + hv \longrightarrow O + O $ | Slow (creates odd oxygen)| $k_1 = 3 \times 10^{-12} \ s^{-1}$ |
|2) $O + O_2 + M \longrightarrow O_3 + M$ | Fast (creates ozone)| $k_2 = 1.2 \times 10^{-33} \ cm^6 molec^{-2} s^{-1}$ |
|3) $O_3 + hv \longrightarrow O_2 + O$ | Fast (destroy ozone) | $k_3 = 5.5 \times 10^{-4} \ s^{-1}$ |
|4) $O + O_3 \longrightarrow 2O_2$ | Slow (destroys odd oxygen)| $k_4 = 6.9 \times 10^{-16} \ cm^3 molec^{-1} s^{-1}$|

</div>

where M is any non-reactive species that can take up the energy released in reaction (2) to stabilize $O_3$. $O_3$ is not a very stable molecule and (without the presence of M) the $O_3$ formed by the collision of $O_2$ and $O$ would immediately fall apart to give back $O$ and $O_2$. Given that $N_2$ and $O_2$ are the major components in the atmosphere, M is either $O_2$ or $N_2$. These reactions lead to the following rate equations for $[O], \ [O_3]$ and $[ O_2]$:

<p align='center'>
    $$\begin{align} \frac{d[O_2]}{dt} &= -k_1 [O_2] - k_2 [O_2] [O] [M] + k_3 [O_3] + 2k_4 [O] [O_3] \quad(10 \\ \frac{d[O]}{dt} &= 2k_1[O_2] - k_2 [O_2] [O] [M] + k_3 [O_3] - k_4 [O] [O_3] \quad(11) \\ \frac{d[O_3]}{dt} &= k_2[O_2][O][M] -k_3[O_3] -k_4[O][O_3] \quad(12) \end{align}$$ 
</p>

the rate constants $k_1$ and $k_3$ depend on light intensity, which in this case is the light intensity of the sun. The rate constants apply at an altitude of 25 km, where $[M] = 9 \times 10^{17} \ molec \ cm^{-3}$.

We can use the __steady-state approximation__ to solve for the concentration of $O$ and $O_3$. The steady-state approximation assumes that after an initial time period, the concentration of the reaction intermediates remain a constant with time, i.e the rate of change of the intermediate's concentration with time is zero. Hence, using the steady state approximation

<p align='center'>
    $$\frac{d[O]}{dt} =\frac{d[O_3]}{dt} = 0 \quad(13)$$
</p>

we can solve for $[O]$ and $[O_3]$

<p align='center'>
    $$\begin{align} [O] &=\frac{2k_1 [O_2] + k_3 [O_3]}{k_2 [O_2] [M]+ k_4[O_3]} \quad(14) \\ [O_3] &=\sqrt{\frac{k_1 k_2}{k_3 k_4}} [O_2] [M]^{\frac{1}{2}} \quad(15) \\ \frac{[O]}{[O_3]} &= \frac{k_3}{k_2 [O_2][M]} \quad(16) \end{align}$$
</p>

The Chapman mechanism have been very well studied and the rate constant for these reactions at different altitudes are known. For example, this code performs numerical integration of the above ODE using `scipy.integrate.solve_ivp` and compare the result with the steady-state values above.

```Python
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Reaction rate constant
k1 = 3.e-12
k2 = 1.2e-33
k3 = 5.5e-4
k4 = 6.9e-16

# Return the dx/dt for each reactions
def deriv(t, c, M):
    O2, O, O3 = c # ozone formation reaction O + O2 + M => O3 + M
    dO2dt = -k1*O2 - k2*O2*O*M + k3 * O3 + 2 *k4 * O * O3
    dOdt = 2*k1*O2 - k2 * O2*O*M + k3*O3 - k4*O*O3
    dO3dt = k2*O2*O*M -k3*O3 -k4*O*O3
    return dO2dt, dOdt, dO3dt

# Total molecule concentration, M, and O2 concentration, cO2
M = 9.e17
cO2 = 0.21 * M #Referring to the concentration or partial pressure of oxygen in the air 20.95% OR 0.21

# Initial conditions for [O2], [O], [O3]
c0 = [cO2, 0, 0]

# Integrate the differential equations over a suitable time grid (s).
ti, tf = 0, 5.e7 # approximately 1.5 years typical

# NB We need a solver that is robust to stiff ODE
soln = solve_ivp(deriv, (ti, tf), c0, args=(M,), method='LSODA')
t, c = soln.t, soln.y

# Steady state approximation solution for comparison
cO3ss = np.sqrt(k1 * k2 / k3 / k4 * M) * cO2
cOss = cO3ss * k3 / k2 /cO2 / M


print('Numerical values:\n[O] = {:g} molec/cm3, [O3] = {:g} molec/cm3'.format(*c[1:, -1]))
print('Steady-state values:\n[O]ss={:g} molec/cm3, [O3]ss= {:g} molec/cm3'.format(cOss, cO3ss))
```
The output shows that the steady-state approximation works well at this altitude:
```
Numerical values:
[O] = 4.69124e+07 molec/cm3, [O3] = 1.7407e+13 molec/cm3
Steady-state values:
[O]ss=4.7055e+07 molec/cm3, [O3]ss= 1.74634e+13 molec/cm3
```

```Python
# Plot the evolution of [O3] and [O] with time
plt.plot(t, c[2],'D', label=r'$\mathrm{[O_3]}$')
plt.plot(t, c[1],'^', label=r'$\mathrm{[O]}$')
plt.yscale('log')
plt.ylim(1.e5, 1.e14)
plt.xlabel(r'$t / \mathrm{s}$')
plt.ylabel(r'$[\cdot] / \ \mathrm{molec \ cm^{-3}}$')
plt.legend()
plt.title('The Chapman cycle for generating ozone')
plt.savefig('The Chapman cycle for generating ozone.svg', bbox_inches ='tight')
plt.show()
```
Clearly, the steady-state approximation works well for this system (Figure 3)

<div align='center'>
<img src="https://github.com/dindagustiayu/First-Order-Differential-Equations-in-Chemistry/blob/main/First%20ODE%20Kinetic%20reactions/The%20Chapman%20cycle%20for%20generating%20ozone.svg">
</div>

<p align='center'>
$$\mbox{Figure 3. Ozone and Oxygen atom concentrations rapidly reach a steady-state under the Chapman Cycle reaction mechanism}$$
</p>

## Conclusion
1. The first-order equations is  a simple modle assumption and defining the boundary conditions to describe many processes in chemistry and physics.
2. In chemical reaction kinetics, first-order differential equations are completely describing the variation $dy$ of a function $y(x)$ and other quantities.
3. In computational, the first-order equations is essential to understand algorithm and numerical method to solve problems. The default algorithm used by `solve_ivp` is `RK45`, an explicit Runge-Kutta method, which is general-purpose approach for many problems.



 
