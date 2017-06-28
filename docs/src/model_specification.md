Model Specification
===================

Variable types
--------------

The following types of variables can be used in models:

> -   `exogenous` (`m`)
> -   `states` (`s`)
> -   `controls` (`x`)
> -   `auxiliaries` (`y`)
> -   `rewards` (`r`)
> -   `values` (`v`)
> -   `expectations` (`z`)
> -   `parameters` (`p`)

Symbol types that are present in a model are always listed in that
order.

### State-space

Decisions are characterized by a vector $m$ of *exogenous* variables (exogenous states)
and by a $s$ of endogenous *states* (exogenous states).
The unknown vector of controls $x$ is a
function $\varphi$ of the states such that:

> $x = \varphi(m,s)$


The function $\varphi$ is typically approximated by the solution
algorithm. It can be either a Taylor expansion, or an interpolating
object (splines, smolyak). Once obtained, it can be evaluated at any point of the
the state space (represented by a couple of vectors), or at a list of points (of couple of matrices):

``` {.sourceCode .julia}
dr = time_iteration(model)
m0 = model.calibration[:states]
s0 = model.calibration[:states]
dr(m0,s0)                               # evaluates at the steady-state
dr([0.0;-0.01;-0.1], [2.5;2.5;2.5])     # evaluates at a list of point
```

A decision rule is defined on two discretized grids: one for the exogenous states
and one for the endogenous ones. If the exogenous grid points are numbered by $i$,
we can define non-ambiguously $x = \varphi(i,s)$ as  $\varphi(m_i,s)$, which corresponds to the following codes:
``` {.sourceCode .julia}
dr(2,s0)                 # evaluates at the second exogenous point
dr(2, [2.5;2.5;2.5])     # evaluates at a list of point
```




### Valid equations

The various equations understood by Dolo are descrbed below:

<!-- would be cool to know how to make comments ;-) -->

#### Transitions

    - name: `transition`
    - short name: `g`

Transitions are given by a function $g$ such that at all times:

$$s_t = g(m_{t-1}, s_{t-1}, x_{t-1}, m_t)$$

where $\m_t$ is a vector-valued exogenous process.

> **note**
>
> In the RBC model, the vector of endogenous states is $s_t=(a_t,k_t)$. The
> transitions are (note )
>
> > $a_t = \rho a_{t-1} + \epsilon_t$
> > $k_t = (1-\delta)k_{t-1} + i_{t-1}$
>
> If $\epsilon_t$ follow a normal law with variance $\sigma_{\epsilon}$, the yaml file is amended with:
>
> ``` {.sourceCode .yaml}
> symbols:
>     states: [a,k]
>     controls: [i]
>     exogenous: [epsilon]
>     ...
> equations:
>     transition:
>         a = rho*a(-1) + e
>         k = k(-1)*(1-delta) + i(-1)
> exogenous: !Normal
>     Sigma: [[sigma_epsilon]]
> ```
>
> Note that transition equations must list states in declaration order.
> Also, in this example productivity process $a_t$ is essentially an exogenous
> AR1 process. Declaring it as such (instead of $\epsilon$), provides additional
information to the solvers and can lead to faster solution time.


#### Auxiliary variables / Definitions

    - name: `auxiliary`
    - short name: `a`

In order to reduce the number of variables, it is useful to define
auxiliary variables $y_t$ using a function $a$ such that:

$$y_t = a(m_t, s_t, x_t)$$

These variables are defined in a special `definitions` block, outside of `equations`.
When auxiliary variables appear in an equation they are automatically substituted by
the corresponding expression in $m_t$,$s_t$ and $x_t$.


> **note**
>
> In the RBC model, three auxiliary variables are defined
> $y_t, c_t, r_{k,t}$ and $w_t$. They are a closed form function of
> $a_t, k_t, i_t, n_t$. Defining these variables speeds up computation
> since they are don't need to be solved for or interpolated.

#### Utility function and Bellman equation

    - name: `utility`
    - short name: `u`

The (separable) value equation defines the value $v_t$ of a given policy
as:

$$v_t = u(m_t,s_t,x_t) + \beta E_t \left[ v_{t+1} \right]$$

This gives rise to the Bellman eqution:

> $$v_t = \max_{x_t} \left( u(m_t, s_t,x_t) + \beta E_t \left[ v_{t+1} \right] \right)$$

These two equations are characterized by the reward function $u$ and the
discount rate $\beta$. Function $u$ defines the vector of symbols
`rewards`. Since the definition of $u$ alone is not sufficient, the
parameter used for the discount factor must be given to routines that
compute the value. Several values can be computed at once, if $U$ is a
vector function and $\beta$ a vector of discount factors, but in that
case in cannot be used to solve for the Bellman equation.

> **note**
>
> Our RBC example defines the value as
> $v_t = \frac{(c_t)^{1-\gamma}}{1-\gamma} + \beta E_t v_{t+1}$. This
> information is coded using: \#\# TODO add labour to utility
>
> ``` {.sourceCode .yaml}
> symbols:
>     ...
>     rewards: [r]
>
> equations:
>     ...
>     utility:
>         - r = c^(1-gamma)/(1-gamma)
>
> calibration:
>     ...
>     beta: 0.96   # beta is the default name of the discount
> ```

#### Value

    - name: `value`
    - short name: `w`

A more general updating equation can be useful to express non-separable
utilities or prices. the vector of (generalized) values $v^{*}$ are
defined by a function `w` such that:

$$v_t = w(m_t,s_t,x_t,v_t,m_{t+1},s_{t+1},x_{t+1},v_{t+1})$$

As in the separable case, this function can either be used to compute
the value of a given policy $x=\varphi()$ or in order solve the
generalized Bellman equation:

$$v_t = \max_{x_t} \left( w(m_t,s_t,x_t,v_t,m_{t+1},s_{t+1},x_{t+1},v_{t+1}) \right)$$

> **note**
>
> Instead of defining the rewards of the RBC example, one can instead
> define a value updating equation instead:
>
> ``` {.sourceCode .yaml}
> symbols:
>     ...
>     values: [v]
>
> equations:
>     ...
>     value:
>         - v = c^(1-gamma)/(1-gamma)*(1-n...) + beta*v(1)
> ```

#### Boundaries

    - name: `controls_lb` and `controls_ub`
    - short name: `lb` and `ub`

The optimal controls must also satisfy bounds that are function of
states. There are two functions $\underline{b}()$ and $\overline{b}()$
such that:

$$\underline{b}(m_t, s_t) \leq x_t \leq \overline{b}(m_t, s_t)$$

> **note**
>
> In our formulation of the RBC model we have excluded negative
> investment, implying $i_t \geq 0$. On the other hand, labour cannot be
> negative so that we add lower bounds to the model:
>
> ``` {.sourceCode .yaml}
> equations:
>     ...
>     controls_lb:
>         i = 0
>         n = 0
> ```
>
> Specifying the lower bound on labour actually has no effect since
> agents endogeneously choose to work a positive amount of time in order
> to produce some consumption goods. As for upper bounds, it is not
> necessary to impose some: the maximum amount of investment is limited
> by the Inada conditions on consumption. As for labour `n`, it can be
> arbitrarly large without creating any paradox. Thus the upper bounds
> are omitted (and internally treated as infinite values).

#### Euler equation

    - name: `arbitrage` (`equilibrium`)
    - short name: `f`

A general formulation of the Euler equation is:

$$0 = E_t \left[ f(m_t, s_t, x_t, m_{t+1}, s_{t+1}, x_{t+1}) \right]$$

Note that the Euler equation and the boundaries interact via
"complentarity equations". Evaluated at one given state, with the vector
of controls $x=(x_1, ..., x_i, ..., x_{n_x})$, the Euler equation gives
us the residuals $r=(f_1, ..., f_i, ...,
f_{n_x})$. Suppose that the $i$-th control $x_i$ is supposed to lie in
the interval $[ \underline{b}_i, \overline{b}_i ]$. Then one of the
following conditions must be true:

-   $f_i$ = 0
-   $f_i<0$ and $x_i=\overline{b}_i$
-   $f_i>0$ and $x_i=\underline{b}_i$

By definition, this set of conditions is denoted by:

-   $f_i = 0 \perp \underline{b}_i \leq x_i \leq \overline{b}_i$

These notations extend to a vector setting so that the Euler equations
can also be written:

$$0 = E_t \left[ f(m_t, s_t, x_t, m_{t+1}, s_{t+1}, x_{t+1}) \right] \perp \underline{b}(m_t, s_t) \leq x_t \leq \overline{b}(m_t, s_t)$$

Specifying the boundaries together with Euler equation, or providing
them separately is exactly equivalent. In any case, when the boundaries
are finite and occasionally binding, some attention should be devoted to
write the Euler equations in a consistent manner. In particular, note
that the Euler equations are order-sensitive.

The Euler conditions, together with the complementarity conditions
typically often come from Kuhn-Tucker conditions associated with the
Bellman problem, but that is not true in general.

> **note**
>
> The RBC model has two Euler equations associated with investment and
> labour supply respectively. They are added to the model as:
>
> ``` {.sourceCode .yaml}
> arbitrage:
>     - 1 - beta*(c/c(1))^(sigma)*(1-delta+rk(1))   | 0 <= i <= inf
>     - w - chi*n^eta*c^sigma                       | 0 <= n <= inf
> ```
>
> Putting the complementarity conditions close to the Euler equations,
> instead of entering them as separate equations, helps to check the
> sign of the Euler residuals when constraints are binding. Here, when
> investment is less desirable, the first expression becomes bigger.
> When the representative is prevented to invest less due to the
> constraint (i.e. $i_t=0$), the expression is then *positive*
> consistently with the complementarity conventions.

#### Expectations

    - name: `expectation`
    - short name: `h`

The vector of explicit expectations $z_t$ is defined by a function $h$
such that:

$$z_t = E_t \left[ h(m_{t+1}, s_{t+1},x_{t+1}) \right]$$

``` {.sourceCode .}
In the RBC example, one can define. the expected value tomorrow of one additional unit invested tomorrow:

.. math::

    m_t=\beta*(c_{t+1}^(-\sigma)*(1-\delta+r_{k,t+1})

 It is a pure expectational variable in the sense that it is solely determined by future states and decisions. In the model file, it would be defined as:

.. code: yaml

    symbols:
        ...
        expectations: [z]

    equations:
        ...
        - z = beta*(c(1))^(-sigma)*(1-delta+rk(1))
```

#### Generalized expectations

    - name: `expectation_2`
    - short name: `h_2`

The vector of generalized explicit expectations $z_t$ is defined by a
function $h^{\star}$ such that:

$$z_t = E_t \left[ h^{\star}(m_t,s_t,x_t,m_{t+1},s_{t+1},x_{t+1}) \right]$$

#### Euler equation with expectations

    - name: `arbitrage_2` (`equilibrium_2`)
    - short name: `f_2`

If expectations are defined using one of the two preceding definitions,
the Euler equation can be rewritten as:

$$0 = f(m_t, s_t, x_t, z_t) \perp \underline{b}(m_t, s_t) \leq x_t \leq \overline{b}(m_t, s_t)$$

> **note**
>
> Given the definition of the expectation variable $m_t$, today's
> consumption is given by: $c_t = z_t^\left(-\frac{1}{\sigma}\right)$ so the Euler
> equations are rewritten as:
>
> ``` {.sourceCode .yaml}
> arbitrage_2:
>     - 1 - beta*(c)^(sigma)/m   | 0 <= i <= inf
>     - w - chi*n^eta*c^sigma    | 0 <= n <= inf
> ```
>
> Note the type of the arbitrage equation (`arbitrage_2` instead of
> `arbitrage`).
>
> However $c_t$ is not a control itself,
>
> > but the controls $i_t, n_t$ can be easily deduced:
>
> ..math:
>
>     n_t = ((1-\alpha) z_t k_t^\alpha m_t/chi)^(1/(eta+\alpha))
>     i_t = z_t k_t^\alpha n_t^(1-\alpha) - (m_t)^(-1/\sigma)
>
> This translates into the following YAML code:
>
> ``` {.sourceCode .yaml}
> equations:
>     - n = ((1-alpha)*a*k^alpha*m/chi)^(1/(eta+alpha))
>     - i = z*k^alpha*n^(1-alpha) - m^(-1/sigma)
> ```

#### Direct response function

    - name: `direct_response`
    - short name: `d`

In some simple cases, there a function $d()$ giving an explicit
definition of the controls:

$$x_t = d(m_t, s_t, z_t)$$

Compared to the preceding Euler equation, this formulation saves
computational time by removing the need to solve a nonlinear system to
recover the controls implicitly defined by the Euler equation.
