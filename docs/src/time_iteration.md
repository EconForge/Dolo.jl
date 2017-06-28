Time iteration
==============

We consider a model with the form:

> $s_t = g\left(m_{t-1}, s_{t-1}, x_{t-1}, m_t \right)$
>
> $0 = E_t \left[ f\left(m_t, s_{t}, x_{t}, m_{t+1}, s_{t+1}, x_{t+1} \right) \right]$

where $g$ is the state transition function, and $f$ is the arbitrage
equation.

The time iteration algorithm consists in approximating the optimal
controls as a function of exogenous and endogenous controls
$x_t = \varphi(m_t,s_t)$.

At step $n$, the current guess for the control, $x(s_t) = \varphi^n(m_t, s_t)$, serves as the control being exercised next period :

:   -   Given current guess, find the current period's
        $\varphi^{n+1}(m_t,s_t)$ controls for any $(m_t,s_t)$ by solving
        the arbitrage equation :
        $0 = E_t \left[ f\left(m_t, s_{t}, \varphi^{n+1}(m_t, s_t), g(s_t, \varphi^{n+1}(m_t, s_t)), \varphi^{n}(m_{t+1},g(s_t, \varphi^{n+1}(s_t))) \right) \right]$


