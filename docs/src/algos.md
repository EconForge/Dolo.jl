# Solution Algorithms

## Value iteration

Value function algorithms require a Bellman representation of the model:

> $s_t = {\color{red} g} \left( m_{t-1}, s_{t-1}, x_{t-1}, m_t \right)$
> $V(m_t, s_t) = \max_{\underline{x}(m_t,s_t) \leq x_t \leq \overline{x}(m_t,s_t)} {\color{red} u}\left(m_t, s_t, x_t\right) + {\color{blue} \beta} E_{m_t} V\left(m_{t+1}, s_{t+1}\right)$

where $g$ is the transition function, $u$ the instantaneous felicity function and $\beta$ the time-discount parameter.

The solution of this problem produces naturally two functions, $x=\varphi(m,s)$ and  $V=\varphi(m,s)$, for the controls and the value function respectively.

Given an initial value function $V^n(m,s)$ applying the maximum operator produces a new and improved value function $\tilde{V}^n(m,s)$ and a corresponding policy rule $x=\varphi^n(m,s)$. This is an *improvement* step.

Note that at this stage, $\tilde{V}^n$ is not the value of following $\varphi^n$ forever. An evaluation step performs the recursion (note the absence of a $\max$):

> $x_t=\varphi^n(m_s, s_t)$
> $s_t = {\color{red} g} \left( m_{t-1}, s_{t-1}, x_{t-1}, m_t \right)$
> $V^{n+1}_{k+1}(m_t, s_t) = {\color{red} u}\left(m_t, s_t, x_t\right) + {\color{blue} \beta} E_{m_t} V^{n+1}_{k}\left(m_{t+1}, s_{t+1}\right)$

Starting from the initial guess, $\tilde{V}^{n+1}$, this recursion converges to the value $V^{n+1}()$ associated to $\varphi^{n+1}()$. If high accuracy is not required, it is common to restrict the number of steps to perform a *partial evaluation*.

The value function algorithm implemented in Dolo follows the following scheme:

1. Given an initial guess for the policy rule $\varphi^0()$, evaluate the corresponding value function $V^0()$.
2. Then given $\varphi_n$:
    - perform one improvement step to get $\varphi^{n+1}()$, $\tilde{V}^{n+1}()$
    - perform a partial evaluation of $\varphi^{n+1}()$
    starting from $\tilde{V}^{n+1}()$ to get $V^{n+1}()$ using at most $K$ steps.
3. Compute $\eta_{n+1}=|\varphi^n-\varphi^{n+1}|$ and $\epsilon_{n+1}=|V^n-V^{n+1}|$
    - if $\eta_{n+1}<\tau_{\eta}$ and $\epsilon_{n+1}<\tau_{\epsilon}$, return
    - else return to step 2.

This algorithm is controlled by the precision parameters $\tau_{\eta}$, $\tau_{\epsilon}$ and the evaluation length $K$. Note that setting $K=0$ corresponds to the naive (but robust ?) VFI algorithm while $K$ high corresponds to Howard improvements which converge faster: convergence of the outer loop $\varphi^n$ is quadratic instead of geometric.

```@docs
value_iteration
```

## Time iteration


We consider a model with the form:

> $s_t = g\left(m_{t-1}, s_{t-1}, x_{t-1}, m_t \right)$
>
> $0 = E_t \left[ f\left(m_t, s_{t}, x_{t}, m_{t+1}, s_{t+1}, x_{t+1} \right) \right]$

where $g$ is the state transition function, and $f$ is the arbitrage
equation.

The time iteration algorithm consists in approximating the optimal
controls as a function of exogenous and endogenous controls
$x_t = \varphi(m_t,s_t)$. At step $n$, the current guess for the control, $x(s_t) = \varphi^n(m_t, s_t)$, serves as the control being used next period.

Here is an outline of the algorithm:

1.   Start with initial guess $\varphi^0$
2.   Given current guess, find the current period's
    $\varphi^{n+1}(m_t,s_t)$ controls for any $(m_t,s_t)$ by solving (numerically)
    the arbitrage equation :
> $0 = E_t \left[ f\left(m_t, s_{t}, \varphi^{n+1}(m_t, s_t), g(s_t, \varphi^{n+1}(m_t, s_t)), \varphi^{n}(m_{t+1},g(s_t, \varphi^{n+1}(s_t))) \right) \right]$

3.  Compute successive approximation errors $\eta_n=|\varphi^{n+1}-\varphi^{n}|$.
    - if $\eta_n$ smaller thatn criterion $\epsilon_{\eta}$, return
    - otherwise return to step 2

```@docs
time_iteration
```

In some cases, the solution of the Euler equation, can be obtained
faster if a direct solution for optimal controls is known as a function expectation as in the following specification:


> $s_t = g\left(m_{t-1}, s_{t-1}, x_{t-1}, m_t \right)$
>
> $z_t = E_t \left[ h\left(m_{t+1}, s_{t+1}, x_{t+1} \right) \right]$
>
> $x_t = d(m_t, s_t, z_t)$

This information can be used by passing the `solver=Dict(:type=>:direct)` option to the time_iteration function,
or by using the devoted function:

```@docs
time_iteration_direct
```

## Improved Time iteration

```@docs
improved_time_iteration
```

## Perfect Foresight



Consider a series for the exogenous process $(m_t)_{0 \leq t \leq T}$.
The perfect foresight problem consists in finding the path of optimal
controls $(x_t)_{0 \leq t \leq T}$ and corresponding states
$(s_t)_{0 \leq t \leq T}$ such that:

> $s_t = g\left(m_{t-1}, s_{t-1}, x_{t-1}, m_t \right)$
>
> $0 = E_t \left( f\left(m_{t}, s_{t}, x_{t}, m_{t+1}, s_{t+1}, x_{t+1}\right) \right) \ \perp \ \underline{u} <= x_t <= \overline{u}$

Special conditions apply for the initial state and controls. Initial
state $s_0$ is given exogenously, or determined so that it corresponds for a steady-state corresponding to $m_0$. Final states and controls are
determined by assuming the exogenous process satisfies $m_t=m_T$ for all
$t\leq T$ and optimality conditions are satisfied in the last period:

$f(m_T, s_T, x_T, m_T,s_T, x_T) \perp \underline{u} \leq x_T \leq \overline{u}$.

We assume that $\underline{u}$ and $\overline{u}$ are constants. This is
not a big restriction since the model can always be reformulated in
order to meet that constraint, by adding more equations.

The stacked system of equations satisfied by the solution is:

>|      Transition         |   Arbitrage           |     
>| :-------------  |:-------------|
>| $s_0 = \overline{s_0}$   | $f(m_0, s_0, x_0, m_1, s_1, x_1) \perp \underline{u} <= x_0 <= \overline{u}$ |
>| $s_1 = g(m_0, s_0, x_0, m_1)$      |  $f(m_1, s_1, x_1, m_2, s_2, x_2) \perp \underline{u} <= x_1 <= \overline{u}$ |  
> | ... | ... |
>| $s_T = g(m_{T-1}, s_{T-1}, x_{T-1}, m_T)$ | $f(m_T, s_T, x_T, m_T, s_T, x_T) \perp \underline{u} <= x_T <= \overline{u}$  |

The system is solved using a nonlinear solver.

```@docs
perfect_foresight
```


## Local Analysis

```@docs
residuals
```

```@docs
find_deterministic_equilibrium
```

```@docs
perturbate
```
