{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "b1331b43-73b5-4a6a-8621-1a2466940a78",
   "metadata": {},
   "outputs": [],
   "source": [
    " ---\n",
    "title: \"Ordinary Differential Equation of Chemistry\"\n",
    "date: \"2026-4-26\"\n",
    "categories: [Kinetic energy]\n",
    "---"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "d8ecf196-3050-4cd7-97f3-c6a447cde473",
   "metadata": {},
   "source": [
    "[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)]()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "0c5d0413-14a1-4fa2-99d2-c6616af18f11",
   "metadata": {},
   "source": [
    "# Ordinary Differential Equation of Chemistry\n",
    "\n",
    "A differential equation (DE) is an equation that relates some derivatives of an unkniwn function to each other and to this function itself. An ordinary differential equation is a (ODE), where the unknown function depends on one variable only. We call this variable the independent variable and the unknown function the dependent variable. If the unknown function depends on several independent variables and the equation involves partial derivatives, then we have a partial differential equation. A chemicallly reacting system without space variations can be described by a set of ordinary with time $t$ as independent variable.\n",
    "\n",
    "The mathematical description of various processes in chemistry and physics is possible by describing then with the help of differential equations which are based on simple model assumptions and defining the boundary conditions.\n",
    "\n",
    "## First-order Differential Equations\n",
    "Many processes and phenomena in chemistry, and generally in sciences, can be described by first-order differential equations. These equations are the most important and most frequently used to describe natural laws. \n",
    "The general form of a first order equation can be defined as:\n",
    "<p align='center'>\n",
    "    $$\\frac{dy}{dt} = f(t, \\ y) \\tag{1}$$\n",
    "</p>\n",
    "\n",
    "where $f(t, \\ y)$ is some functions of the independent variable (e.g., time $t$) and $y$ itself, the solution is sought, $y(t)$, in general, varies with time in some way determined by an _initial condition_: a known value of $y$ at some particular time (e.g., $t =0$).\n",
    "\n",
    "## Chemical reaction kinetics\n",
    "\n",
    "Chemical reaction kinetics is the study of rates of chemical processes (reactions). The goal is to find the relations between the concetrations of educts or products of a chemical reaction (as depending variable) and the time (as independent variable). In general, _all_ chemical reactions can be described mathematically by first-order differential equations. Their solutions, however, depend directly on the nature of the chemical reaction itself. The latter is characterized  by the so-called _reaction order_, which has nothing what so ever to do with the _order of a differential equation_. The reaction order of a chemical reaction is simply defined by the sum of exponents of concentrations occuring in th rate law.\n",
    "\n",
    "Simple reactions like the transformation of $A \\rightarrow P$ can be described by the differential equation:\n",
    "<p align='center'>\n",
    "    $$-d[A] = k \\cdot A \\cdot dt \\tag{2}$$\n",
    "</p>\n",
    "\n",
    "The first order kinetics are summarized by the differential equation:\n",
    "<p align='center'>\n",
    "    $$\\frac{d[A]}{dt} = -k [A]$$\n",
    "</p>\n",
    "\n",
    "where $[A]_0$ is the initial concentration of $[A]$.\n",
    "\n",
    "\n",
    "The _first-order differential equation_ describes also a _first-order reaction_ in chemical kinetics, due to the exponent 1 of the concentration $A$. Following a separation of variables, the integration results in:\n",
    "\n",
    "<p align='center'>\n",
    "    $$ \\int_{A_0}^{A_t} \\frac{d[A]}{[A]} = -k \\int_0^{t} \\ dt \\tag{3}$$\n",
    "</p>\n",
    "The result\n",
    "\n",
    "<p align='center'>\n",
    "    $$ \\ln \\frac{A_t}{A_0} = -k t \\tag{4}$$\n",
    "</p>\n",
    "\n",
    "can be also written in the exponential form:\n",
    "<p align='center'>\n",
    "    $$[A](t) = [A]_0 e^{-k t} \\tag{5}$$\n",
    "</p>\n",
    "\n",
    "\n",
    "## Preliminaries\n",
    "\n",
    "Ordinary differential equations can be solved numerically with,\n",
    "- `scipy.integrate.solve_ivp`: solve an initial value problem.\n",
    "- `scipy.integrate.odeint`: solve first-order differential equations.\n",
    "\n",
    "This function is based on the well-tested Fortran LSODA routine, which can automatically switch between stiff anf nostiff algorithms.\n",
    "\n",
    "`solve_ivp` or `odeint` take three arguments:\n",
    "- A function object returning $dy/dt$ as $f(t, \\ y)$, given $t$ and $y$.\n",
    "- An initial condition, $y_0$.\n",
    "- A sequence of $t$, the initial and final time values at which to calculate the solution, $y(t)$, as a tuple ($t0, \\ tf$).\n",
    "\n",
    "The derivative function is simply:\n",
    "```\n",
    "def dydt(y, t):\n",
    "    return -k * y\n",
    "```\n",
    "the order of the argument is important. Or \n",
    "\n",
    "```\n",
    "def dydt(y, t):\n",
    "    return -k * y\n",
    "\n",
    "soln = solve_ivp(dydt, (t0, tf), [y0])\n",
    "```\n",
    "The returned object, `soln`, contains information about the solution.\n",
    "\n",
    "### Example - The decomposition of $N_2O_5$\n",
    "The decomposition of $N_2O_5$ is first order with a rate constant, $k = 5 \\times 10^{-4} \\ s^{-1}$. Given an initial consentration of $[N_2O_5] = 0.02 \\ mol \\ dm^{-3}$, what is the concentration of $N_2O_5$ after 2 hours?.\n",
    "\n",
    "### Solution\n",
    "1. Write a single dependent variable, $y(t) \\equiv [A]$, which is a function of the independent variable, $t$ (time). We have:\n",
    "<p align='center'>\n",
    "    $$ \\frac{dy}{dt} = -ky$$\n",
    "</p>\n",
    "\n",
    "2. We need to provide a function returning $dy/dt$ as $f(y, \\ t)$ (in general a function of both $y$ and $t$), an intitial condition, $y(0)$ and a sequence of time points upon which to calculate the solution.\n",
    "   \n",
    "3. Calculate concentration through analytical equation, $[A](t) = [A]_0 e^{-k t}$.\n",
    "\n",
    "4. Determinate the value of $y(t_f)$ by the ODE solver algorithm.\n",
    "\n",
    "5. Set the `dense_Output` argument to `True` defines an `Odesolution`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 14,
   "id": "319d6ee8-2712-4d29-9e59-f68c069f884c",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "  message: The solver successfully reached the end of the integration interval.\n",
      "  success: True\n",
      "   status: 0\n",
      "        t: [ 0.000e+00  4.618e-01  5.080e+00  5.126e+01  5.130e+02\n",
      "             2.376e+03  4.149e+03  5.960e+03  7.200e+03]\n",
      "        y: [[ 2.000e-02  2.000e-02  1.995e-02  1.949e-02  1.547e-02\n",
      "              6.101e-03  2.515e-03  1.018e-03  5.475e-04]]\n",
      "      sol: None\n",
      " t_events: None\n",
      " y_events: None\n",
      "     nfev: 50\n",
      "     njev: 0\n",
      "      nlu: 0\n"
     ]
    }
   ],
   "source": [
    "import numpy as np\n",
    "from scipy.integrate import solve_ivp\n",
    "import matplotlib.pyplot as plt\n",
    "\n",
    "# First-order reaction rate constant, s-1\n",
    "k = 5.0e-4\n",
    "\n",
    "# Initial condition on y: 100% of N2O5 is present at t = 0\n",
    "y0 = 0.02\n",
    "\n",
    "# Initial and final time points for the integration (s).\n",
    "t0, tf = 0, 2*60 * 60\n",
    "\n",
    "# 1. Return dy/dt = f(y, t) at t time\n",
    "def dydt(t, y):\n",
    "    return -k * y\n",
    "\n",
    "# 2. Integrate the differential equation\n",
    "soln = solve_ivp(dydt, (t0, tf), [y0])\n",
    "print(soln)"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "c8cab1c5-22f5-48a8-8bbf-b70be0c1e312",
   "metadata": {},
   "source": [
    "The solution was obtained succesfully in 50 function evaluations and the initial, intermediate and final values of $t$ and $y$ are reported. Next, we need calculate consentration of $N_2O_5$ after 2 hours."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 15,
   "id": "364c6ed1-9582-4ae2-bf95-b15e06e93f57",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "N205:0.0005464744489458512\n"
     ]
    }
   ],
   "source": [
    "y0 = y0 * np.exp(-k * tf)\n",
    "print(f'N205:{y0}')"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4a4f6757-9ba0-486e-b019-33f04f36d551",
   "metadata": {},
   "source": [
    "The concentration of $N_2O_5$ after 7200 secs is $5.46 \\times 10^{-4} \\ mol \\ dm^{-3}$. "
   ]
  },
  {
   "cell_type": "markdown",
   "id": "039bb5e2-625a-4019-ba3a-818a705f1d20",
   "metadata": {},
   "source": [
    "The twenty one time points, given as the array $t$ in the returned solution, that were used to determine the value of $y(t_f)$ were chosen by the ODE solver algorithm, To follow the reactant concentration as a function of time in higher resolution, provide a suitable sequence of time points as the argument `t_eval`."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 16,
   "id": "06b88236-adf2-43de-9a9b-122bfca47efd",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "t: [7200. 6840. 6480. 6120. 5760. 5400. 5040. 4680. 4320. 3960. 3600. 3240.\n",
      " 2880. 2520. 2160. 1800. 1440. 1080.  720.  360.    0.]\n",
      "y: [0.00054647 0.00065425 0.00078329 0.00093786 0.00112272 0.00134387\n",
      " 0.00160875 0.00192618 0.00230643 0.00276162 0.00330668 0.00395789\n",
      " 0.00473649 0.00566897 0.00678709 0.00812786 0.00973354 0.01165358\n",
      " 0.01395182 0.01670321 0.01999755]\n"
     ]
    }
   ],
   "source": [
    "# 4. Integrate the differential equation, report at a 21 time points.\n",
    "t_eval = np.linspace(tf, t0, 21)\n",
    "soln = solve_ivp(dydt, (tf, t0), [y0], t_eval = t_eval)\n",
    "\n",
    "t, y = soln.t, soln.y[0]\n",
    "\n",
    "print(f't: {t}')\n",
    "print(f'y: {y}')"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 17,
   "id": "4740d1e7-689d-4d17-981a-cab5b01bb945",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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21+2bFlrBwAA/oshP/gpV2rw4MF2XcHTxYsX3cKFC+22F0S9//771KMCACAOAirE0qVLFzd27FhXunRpu71o0SILrEqUKOE+++wzux8AAMRGpXTEo6CpU6dONptPCega9tPsP690AgAAiI2ACokO/7Vq1cqunzp1yoYCd+7c6Xbs2OHKlSvHqwYAQACG/HBNgaUT5syZwysGAEAcBFRIFq90wubNm6lDBQBAHARUSJZChQq5mjVr2vW5c+fyqgEAEICACsl266232uXatWvdkSNHeOUAAPj/EVAh2YoXL+4qV65s1+mlAgDgvwiokCJegc+VK1faWn8AAICACimkpWdUNuHq1atu3rx5vH4AABBQIS25VMuWLbOinwAARDuG/JBilSpVsqVoLl265F/rDwCAaEZAhRRTPSovl0pr/V24cIFXEQAQ1QiokCrVqlVzhQsXdufPn3dLlizhVQQARDUCKqTuFydLFqueLgsWLHCXL1/mlQQARC0CKqRanTp13PXXX+9Onz7tVqxYwSsJAIhaBFRItaxZs7omTZr4C32qlAIAANGIgApp0qBBA5c7d253/PhxW5IGAIBoRECFNImJiXGNGze263PmzHE+n49XFAAQdQiokGY33XSTBVYHDx50v//+O68oACDqEFAhzXLlyuUaNWpk1+mlAgBEIwIqpAslpytJfffu3W7Hjh28qgCAqEJAhXSRN29eV79+fX8vFQAA0YSACummadOmtizNli1b3N69e3llAQBRg4AK6aZgwYKudu3adp1eKgBANCGgQrrylqNZv369O3z4MK8uACAqEFAhXRUrVsxVrVrVXz0dAIBoQECFdHfrrbfa5apVq9yJEyd4hQEAES+kAqqhQ4e6ChUquJw5c7qGDRu63377Lcn2s2fPtnZqX7FiRTd8+PB4bcaNG+dq1KjhcuTIYZcTJkxI0/M+9dRTlnj9/vvvp/IsI1+ZMmXs9dTafvPmzQv24QAAED0B1ZgxY1z//v3dSy+95JYvX+6aN2/u2rdv73bu3Jlg+23btrm77rrL2qn9oEGDXN++fS2A8syfP99169bN9ejRw61cudIuH3jgAbdw4cJUPe93331nP1uqVKkMehUir5dq2bJl7syZM8E+HAAAMtR1vhBZfE3rwWmh3WHDhvn3Va9e3XXu3Nm99dZb8dq/8MILbuLEiZb87OnT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22/ujjvu4AUHAIQkMsoRsrJmzWoJ6rJgwQJ37NixYB8SAABp76HSUjPJRQ4V0oNKKFSsWNFm+/38889WTBYAgLAOqJQcHEhrrik5vWrVqnZb09zVq0ACMdKL8vHUS/XRRx+5devWuR07drhy5crxAgMAwjegCsybUg9Uvnz53Oeff+4KFixo+zQk89hjj7nmzZun/5EiahUvXtw1aNDAAnglqffu3ZuJDwCAyFjLr3Tp0lbJumbNmrH2a8afehT27t2bXscYEVjLL23OnDnj/vWvf7kLFy64Tp06uXr16qXTOwMAQNo/v7Ok5QkOHDgQb7/qBp06dSq1DwskKE+ePP6ez19++cUq9QMAECpSHVDde++9Nrw3duxYt3v3btt0/Q9/+IPr0qVL+h4l4Jxr3LixDS+fPn3azZ07l9cEABD+AdXw4cNdhw4d3COPPGJJwtoefvhh1759ezd06ND0PUpACX/ZsvlrUc2bN8+6XwEACOscqsDcli1btjg9zI033mhDM4iPHKr0od8zTYTQbL/atWvTGwoACO8cKo8CKK3lV7duXYIpZEoZhbZt2/qXP9JQMwAAEVcpXUMyKsQIZJSSJUv6Z/mpjEIaO1kBAAi9gErJ6j179kzvhwViue2221z27Nmthyru4twAAITN4siJeeaZZ9L7IYF4VFT21ltvtWKzWpKmWrVqFmABABAMLI6MsNWkSRNLEFTC4Pz584N9OACAKMbiyAhb6pFSzt748ePdnDlzXP369a3nCgCAsFocOamZWEBmqFWrllu0aJHlUqmCeoECBdy+ffsscV2V1bVYNwAAIV+HCslDHaqMo2BqxIgRNtvv448/toBKypQp4wYPHkytKgBAaNehOn78uPvHP/7hnnjiCde7d2/3z3/+k+rVyHTqoVq1alWsGlWyZ88e17VrVxsSBAAgI6U6oFqyZImrVKmSBVFHjx51hw8fdu+9957tW7ZsWfoeJZCIK1euuH79+tlMv0uXLrny5cu76tWr231e52v//v2tHQAAIRdQDRgwwHXs2NFt377degAmTJjgtm3b5u6++277AAMyw2+//WZDfuqS1fp+0qZNG38JBQVVu3btsnYAAIRkD9ULL7xgC9Z6dP0vf/mL3QdkBi9fSubOnWtDzgULFnS33357ou0AAAiZgEqJWTt37oy3X70BTF1HZtFsPs/Fixfd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qq+jCtFxe4U0oylu75FHie8Jtc+WLZsl8CbVxnvMlDyvgikltOvDa/r06desJaXeA7UJ3ID0pF5P9ZRmz57dbd682dbWVI8rACAMAyot8FqlShWbcad1/AJn0CnXIzn0LVvlChSoBNLtpk2bJvgzWng5bvtp06a5Ro0a2QdMUm28x0zu83rB1KZNmyx3zAvYgGArV66cLVujLxKafapZrgRVABAEvjSqW7eu7+DBg3Z98eLFvho1avi+/vpru12vXr1kP87o0aN92bNn940YMcK3bt06X//+/X158uTxbd++3e5/8cUXfT169PC337p1qy937ty+AQMGWHv9nH5+7Nix/jZz5871Zc2a1ff222/71q9fb5fZsmXzLViwINnPe+nSJV/Hjh19ZcqU8a1YscK3b98+/3bhwoVkn9+JEyc0HmmXQHrbtGmT73//9399r7/+uv0fuHLliu/y5cu+mTNn+r755hu71G0AQMok9/M7zQGVAqhAhw8f9rVo0cL317/+1Ve/fv0UPdaQIUN85cqV88XExPgaNGjgmz17tv++nj17+lq2bBmr/axZs+w51L58+fK+YcOGxXvMb7/91le1alULmqpVq+YbN25cip5327Zt9kImtOlDKrkIqJDRNm7c6HvjjTcsqHr//fd9ZcuWjfX7qi8FCf3+AwDS/vmdqjpUgVq3bu0GDx5sw30eTePWchnK6SBR9v9QhwqZQXXZlKAumoH7ww8/xJusoWH6Ll268IYAQDAXR45r9+7dlr+RUHX0uXPnumbNmqXl4SMGARUyg/Kn9CVHW5YsWawI7aRJk/z3K6gqU6aMTa7QhAwAQJAKe2r2mxLQtRSL1sVTCYHElpohmAIyl9ah1Pbdd99Zz5RKmdxzzz0WXIn2qYTHtdbJBACkTLYUtnevvfaarXmnb73vvvuuDempSnr9+vVtxpyWcdFWqVKllD40gDTat2+fXa5atcqCqE6dOtn/S3270hC8KvgHtgMABCmg+uMf/xgrV0pr52n5FgVZKjfw/vvvW5kBcqeA4FRR92glA1X2V77UjTfe6B5//HH3zTffWLd1YDsAQNqlOYdKFED99NNPluyqBYuVx6E/2vgvcqiQGfR/r3z58lZo1/uvreBJBUDz5cvnTp065WbMmGG5VeRQAUCQF0cWffNVDpWWvyhatKh9+9UQw5dffukOHTqU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Hn3vO+959L3n0fy+R+t5h9q5ayBPwZRWTcmSJfG5fPRQZRK9CWXKlMmwx9cvXLB/6YKFc4++9533PPre82h+36P1vEPp3JPqmQqpsgkAAADhjIAKAAAgjQiowpyWt3nttdeicpkbzj363nfe8+h7z6P5fY/W8w7XcycpHQAAII3ooQIAAEgjAioAAIA0IqACAABIIwIqAACANCKgCnNDhw51FSpUcDlz5nQNGzZ0v/32mwsnv/76q7vnnnusAq0qyH/33XfxKtS+/vrrdn+uXLlcq1at3Nq1a2O1uXDhgvvTn/7kihQp4vLkyeM6duzodu/eHavNsWPHXI8ePaw4mzZd17pMwfLWW2+5m266ySrnFytWzHXu3Nlt3LgxKs592LBhrk6dOv6CfU2aNHFTpkyJ+PNO6HdAv/P9+/eP+HPXOelcA7cSJUpE/HnLnj173COPPOIKFy7scufO7erVq+eWLl0a8edevnz5eO+5tmeeeSZyzzuzVmtG+hs9erQve/bsvk8++cS3bt06X79+/Xx58uTx7dixI2xe7smTJ/teeukl37hx42w17wkTJsS6/+233/bly5fP7l+9erWvW7duvpIlS/pOnjzpb9OnTx9f6dKlfdOnT/ctW7bM17p1a1/dunV9ly9f9rdp166dr1atWr558+bZput33323L1jatm3rGzlypG/NmjW+FStW+Dp06OC74YYbfKdPn474c584caJv0qRJvo0bN9o2aNAg+z3WaxHJ5x1o0aJFvvLly/vq1Klj/289kXrur732mq9mzZq+ffv2+beDBw9G/HkfPXrUV65cOV+vXr18Cxcu9G3bts33888/+zZv3hzx537w4MFY77eOXX/jZ86cGbHnTUAVxm6++Wb7hQtUrVo134svvugLR3EDqqtXr/pKlChh//E858+f9+XPn983fPhwu338+HH7MFZw6dmzZ48vS5YsvqlTp9ptBZt67AULFvjbzJ8/3/Zt2LDBFwr0x0fHM3v27Kg7dylYsKDv008/jYrzPnXqlK9y5cr2IdGyZUt/QBXJ566ASh+ECYnk837hhRd8t956a6L3R/K5x6Xf80qVKtk5R+p5M+QXpi5evGjdxm3atIm1X7fnzZvnIsG2bdvc/v37Y52jiry1bNnSf456DS5duhSrjbqQa9Wq5W8zf/586wpu3Lixv80tt9xi+0Ll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",
      "text/plain": [
       "<Figure size 640x480 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# 5. Solve the ODE\n",
    "soln = solve_ivp(dydt, (t0, tf), [y0], dense_output=True)\n",
    "\n",
    "# Evaluate the solution, y(t), at 21 time points\n",
    "t= np.linspace(t0, tf, 21)\n",
    "y = soln.sol(t)[0]\n",
    "\n",
    "# Plot and compare the numerical and extract solutions\n",
    "plt.plot(t, y, 'o', color='k', label=r'solve_ivp')\n",
    "plt.plot(t, y0 * np.exp(-k * t), color='gray', label='exact')\n",
    "plt.xlabel('t/s')\n",
    "plt.ylabel(r'$[N_2 O_5]$(t)  / mol.dm-3')\n",
    "plt.legend\n",
    "plt.savefig('Numerical and Analytical ODE.svg', bbox_inches='tight')\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "544a8a0c-97c9-4928-8610-b1b1bd62996b",
   "metadata": {},
   "source": [
    "__Figure 1__. Numerical anad analytical solutions to the ordinary differential equation governing the decomposition of $[N_2 O_5]$. \n",
    "\n",
    "The resulting plot demonstrates that the numerical algorithm was able to follow the true solution accurately in this case."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "6941d795-c70b-4d44-9055-ce1a5fd7756a",
   "metadata": {},
   "source": [
    "## Coupled First-Order Differential Equations\n",
    "\n",
    "In chemistry we are usually concerned with systems of DE of the first order in more than one dependent variable: $y_1(t), \\ y_2(t), \\ldots, \\ y_n(t)$:\n",
    "\n",
    "<p align='center'>\n",
    "    $$\\begin{align} \\frac{dy_1}{ft} &=f_1 (y_1, y_2, \\ldots, y_n;t), \\\\ \\frac{dy_2}{dt} &= f_2(y_1, y_2, \\ldots, y_n;t), \\\\ \\vdots \\\\ \\frac{dy_n}{dt} &= f_n (y_1, y_2, \\ldots, y_n;t) \\end{align}$$\n",
    "</p>\n",
    "\n",
    "In this case, the function passed to `solve_ivp` must return a sequence of derivatives, $dy_1/dt, \\ dy_2/dt, \\ldots, \\ dy_n/dt$ for each of the dependent variables, that is, it evaluates the earlier mentioned functions $f_i(y_1, y_2, \\ldots, y_n;t)$ for each of the y_i passed to it in a sequence, y. \n",
    "\n",
    "In majority of chemical reactions, however, more than one educt is involve. A reaction proceeds via two first-order reaction steps:\n",
    "<p align='center'>\n",
    "    $$ A \\longrightarrow B \\longrightarrow P $$\n",
    "</p>\n",
    "\n",
    "with rate constants  $k_1$ and $k_2$. \n",
    "\n",
    "<p align='center'>\n",
    "    $$ \\begin{align} A \\rightarrow B \\quad k_1 \\\\ B \\rightarrow P \\quad \\end{align}$$\n",
    "</p>\n",
    "\n",
    "The equations governing the rate of change of A and B are\n",
    "<p align='center'>\n",
    "    $$\\begin{align} \\frac{d[A]}{dt} &= -k_1 [A] \\\\ \\frac{d[B]}{dt} &= k_1[A] -k_2 [B] \\end{align}$$\n",
    "</p>\n",
    "\n",
    "We can solve this pair of coupled equations analytically, but in our numerical solution, let $y_1 \\equiv [A]$ and $y_2 \\equiv [B]$:\n",
    "<p align='center'>\n",
    "    $$\\begin{align} \\frac{dy_1}{dt} &= -k_1 y_1 \\\\ \\frac{dy_2}{dt} &= k_1 y_1 - k_2 y_2 \\end{align}$$ \n",
    "</p>\n",
    "\n",
    "The derivative function would then be difined as:\n",
    "```\n",
    "def deriv(t, y):\n",
    "    \"\"\"return dy_i/dt for each y_i at time t,\"\"\"\n",
    "    y1, y2 = y\n",
    "    dy1dt = -k1 * y1\n",
    "    dy2dt = k1 * y1 - k2* y2\n",
    "    return dy1dt, dy2dt\n",
    "```\n",
    "\n",
    "The code mentioned here integrates equations above for $k_1 = 0.2 \\ s^{-1}$, $k_2 = 0.8 \\ s^{-1}$ and initial conditions $y_1(0) = 100, \\ y_2(0) = 0 $, and compares with the analytical result (__Figure 1__)."
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 18,
   "id": "55fa5b03-17ed-4983-8b14-f122d8cf40af",
   "metadata": {},
   "outputs": [
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "The solver successfully reached the end of the integration interval.\n"
     ]
    },
    {
     "data": {
      "image/png": "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dpfyHtV+/S6m+uwwpa8jg6AG9uGHH5piYmJMLi4upk6dOtWYyqp248aNk7+4ItCLiIgwjR492nTgwAGTWom/cERgcO4mpqxapnaKsgBieqerq6upX79+8g+wvbxH8UU6dOhQU3BwsLym0dHR8vipU6dManGh9ya2+fPnW89R83X8r/en9mt41113Wf++FO9h0KBB1gBI7dfuUt6j2q/fpQZBSl5HJ/F/9d/fRERERGRbmBNEREREDolBEBERETkkBkFERETkkBgEERERkUNiEEREREQOiUEQEREROSQGQUREROSQGAQRERGRQ2IQRERERA6JQRARERE5JAZBRGS3BgwYgClTplzwsTvuuAPTp09v8DYRke3QKd0AIqKGZjQasXz5cixdupQfPpEDY08QEdkl0dOzdu1azJ49G05OTnJLTEyUj23cuBEajQbdu3eX93/66Se0bdsW7u7uCAwMxODBg1FcXKzwOyCi+saeICKySyL4OXr0KOLj4/Hyyy/LY8HBwfJW9ACNHDlSBkKpqakYP348Zs6ciRtuuAGFhYVYv349TCaTwu+AiOobgyAisku+vr5wcXGBh4cHwsLCajwmgqC3335b7osgSK/XY/To0YiJiZHHRK8QEdk/DocRkUM5dOgQUlJS5JCX0L59ewwaNEgGPjfddBM+/fRT5ObmKt1MImoADIKIyKGIXqAhQ4bI/B9Bq9Vi5cqV+O2339C6dWt88MEHaNGiBRISEpRuKhHVMwZBRGS3xHCYwWCoceyXX37BddddV+OYSJru3bs3XnrpJezevVs+b/HixQ3cWiJqaMwJIiK71bhxY2zdulXOCvPy8pJT47dv344lS5ZYzxGPr1q1CkOHDkVISIi8n5mZiVatWinadiKqf+wJIiK79cQTT8jhLjHMJWaG/frrr3JavAh2LHx8fLBu3Tpcc801aN68OZ599lm88847GD58uKJtJ6L652TiPFAichBiGKxPnz6YNm2a0k0hIhvAniAichgiABI1gYiIBPYEERERkUNiTxARERE5JAZBRERE5JAYBBEREZFDYhBEREREDolBEBERETkkBkFERETkkBgEERERkUNiEEREREQOiUEQERERwRH9Py0sJi2RJhW7AAAAAElFTkSuQmCC",
      "text/plain": [
       "<Figure size 640x480 with 1 Axes>"
      ]
     },
     "metadata": {},
     "output_type": "display_data"
    }
   ],
   "source": [
    "# First-order reaction rate constant, s-1\n",
    "k1, k2 = 0.2, 0.8\n",
    "\n",
    "# Initial concentration of species A and B (mol. dm-3)\n",
    "y0 = A0, B0 = 100, 0\n",
    "\n",
    "# Initial and final time points for the integration (s)\n",
    "t0, tf = 0, 40\n",
    "\n",
    "def deriv(t, y, k1, k2):\n",
    "    y1, y2 = y\n",
    "    dy1dt = -k1 * y1\n",
    "    dy2dt = k1 * y1 - k2 * y2\n",
    "    return dy1dt, dy2dt\n",
    "\n",
    "# Integrate the differential equation.\n",
    "soln = solve_ivp(deriv, (t0, tf), y0, dense_output=True, args=(k1, k2))\n",
    "print(soln.message)\n",
    "\n",
    "t = np.linspace(t0, tf, 200)\n",
    "A, B = soln.sol(t)\n",
    "\n",
    "# The concentration, [P], is determined by conservation\n",
    "P = A0- A -B\n",
    "\n",
    "plt.plot(t, A, ls='-', label=r'$ \\mathrm{[A]}$')\n",
    "plt.plot(t, B, ls=':' , label=r'$ \\mathrm{[B]}$')\n",
    "plt.plot(t, P, ls = '-.', label=r'$ \\mathrm{[P]}$')\n",
    "plt.xlabel(r'$ \\mathrm{t / s}$')\n",
    "plt.ylabel(r'$ \\mathrm{conc / mol \\ dm^{-3}}$')\n",
    "plt.legend()\n",
    "plt.savefig('ODE in more than on dependent variables.svg', bbox_inches='tight')\n",
    "plt.show()"
   ]
  },
  {
   "cell_type": "markdown",
   "id": "4887ca1f-c084-44b1-be1e-5288fe4c3619",
   "metadata": {},
   "source": [
    "__Figure 2__. Numerical solution to the ODE governing the reaction the sequence $A \\rightarrow B \\rightarrow P$ for rate constant $k_1 = 0.2 \\ s^{-1}, \\ k_2 = 0.8 \\ s^{-1}$ and taking $[A]_0 = 100 \\ mol \\ dm^{-3}, \\ [B]_0 = [P]_0 = 0$."
   ]
  },
  {
   "cell_type": "markdown",
   "id": "b397b462-8a2b-41b7-8a11-a0685a940bdb",
   "metadata": {},
   "source": [
    "## Conclusion\n",
    "1. The first-order equations is  a simple modle assumption and defining the boundary conditions to describe many processes in chemistry and physics.\n",
    "2. In chemical reaction kinetics, first-order differential equations are completely describing the variation $dy$ of a function $y(x)$ and other quantities.\n",
    "3. In computational, the first-order equations is essential to understand algorithm and numerical method to solve problems. The default algorithm used by `solve_ivp` is `RK45`, an explicit Runge-Kutta method, which is general-purpose approach for many problems.\n"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "e8135c40-2855-4d1c-aa20-3a32c3c78410",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.14.2"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
