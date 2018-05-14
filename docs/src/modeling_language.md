The dolo language
=================

The easiest way to code a model in dolo consists in using specialized
Yaml files also referred to as dolo model files.

YAML format
-----------

YAML stands for Yet Another Markup Language. It is a serialization
language that allows complex data structures in a human-readable way.
Atomic elements are floats, integers and strings. An ordered list can be
defined by separating elements with commas and enclosing them with
square brackets:

``` {.sourceCode .yaml}
[1,2,3]
```

Equivalently, it can be done on several lines, by prepending - to each
line

``` {.sourceCode .yaml}
- 'element'
- element         # quotes are optional there is no ambiguity
- third element   # this is interpreted as ``'third element'``
```

Associative arrays map keys(simple strings) to arbitrary values as
in the following example:

``` {.sourceCode .yaml}
{age: 18, name: peter}
```

Mappings can also be defined on several lines, and structures can be nested by using indentation (use spaces no tabs):

``` {.sourceCode .yaml}
age: 18
name: peter
occupations:
  - school
  - guitar
friends:
  paula: {age: 18}
```

The correspondance between the yaml definition and the resulting Julia
object is very transparent. YAML mappings and lists are converted to
Julia dictionaries and arrays respectively.


Special objects from the Dolo language can be created by adding a tag to yaml nodes
as in the following examples:
>
> ``` {.sourceCode .yaml}
> - !AR1:
>     rho: 0.9
>     Sigma: [[0.1]]
> ```

or

> ``` {.sourceCode .yaml}
> - !Product:
>      !AR1:
>         rho: 0.9
>         Sigma: [[0.1]]
>      !AgingProcess:
>           mu: 0.01     # death probability
>           K: 10        # number of ages
> ```


Any model file must be syntactically correct in the Yaml sense, before
the content is analysed further. More information about the YAML syntax
can be found on the [YAML website](http://www.yaml.org/), especially
from the [language specification](http://www.yaml.org/).



Example
-------

Here is an example model contained in the file
`examples\models\rbc.yaml`

This model can be loaded using the command:

``` {.sourceCode .python}
model = yaml_import("examples\models\rbc.yaml")
```

The function `yaml_import` will raise errors until the model
satisfies basic compliance tests. \[more of it below\]. In the following
subsections, we describe the various syntactic rules prevailing while
writing yaml files.


Sections
--------

A dolo model consists of the following:

- a `symbols:` section where all symbols used in the model must be
    defined
- an `equations:` section containing the list of equations
- a `definitions:` section contains variables that can be substituted directly into equations
- a `calibration:` section providing numeric values for the symbols
- a `domain:` section specifying the boundaries for the state variables
- an `exogenous:` section defining the kind of shocks which drive the model
- an `options:` section containing additional information

Sections `symbols:`, `equations:`, `calibration:` should all be present. Other sections can be omitted for specific applications.

These section have context dependent rules. We now review each of them
in detail:

### Declaration section

This section is introduced by the symbols keyword. All symbols appearing
in the model must be defined there.

Symbols must be valid Julia identifiers (alphanumeric not beginning
with a number) and are case sensitive. Greek letters (save for lambda
which is a keyword) are recognized. Subscripts and superscripts can be
denoted by \_ and \_\_ respectively. For instance beta\_i\_1\_\_d will
be printed nicely as $\beta_{i,1}^d$.

Symbols are sorted by type as in the following example:

``` {.sourceCode .yaml}
symbols:
  states: [a, b]
  controls: [c, d]
  exogenous: [e]
  parameters: [rho]
```

Note that each type of symbol is associated with a symbol list (as
\[a,b\]).

> **note**
>
> A common mistake consists in forgetting the commas, and using spaces
> only. This doesn't work since two symbols are recognized as one.

For some applications, other types of symbols can be declared, for instance, one can provide
a specific list of variable for expectations or values (see [Variable types](@ref)).

### Declaration of equations

The equations section contains blocks of equations sorted by type.

Epxressions follow (roughly) the Dynare conventions. Common arithmetic
operators (+,-,\*,/,\^) are allowed with conventional priorities as well
as usual functions
(sqrt, log, exp, sin, cos, tan, asin, acos, atan, sinh, cosh, tanh, asinh, acosh, atanh).
All symbols appearing in an expression must either be declared
in the symbols section or be one of the predefined functions. Any symbol
s that is not a parameter is assumed to be considered at date t. Values
at date t+1 and t-1 are denoted by s(1) and s(-1) respectively.

All equations are implicitly enclosed by the expectation operator
$E_t\left[\cdots \right]$. Consequently, the law of motion for the
capital

$$k_{t+1} = (1-\delta) k_{t} +  i_{t} + \epsilon_t$$

is written as:

``` {.sourceCode .yaml}
k = (1-delta)*k(-1) + i(-1)
```

while the Euler equation

$$E_t \left[ 1=\beta \left( \frac{c_{t+1}}{c_t} + (1-\delta)+r_{t+1} \right) \right]$$

is translated by:

``` {.sourceCode .yaml}
1 = beta*(c/c(1))^(sigma)*(1-delta+rk(1))
```

An equation can consist of one expression, or two expressions separated
by =. There are two types of equation blocks:

-   condition blocks

    > In these blocks, each equation `lhs = rhs` define the scalar value
    > `(rhs)-(lhs)`. A list of of such equations, i.e a block, defines a multivariate function of the appearing symbols.
    > Certain condition blocks, can be associated with complementarity conditions separated by |\`
    > as in `rhs-lhs | 0 <= x <= 1`. In this case it is advised to omit
    > the equal sign in order to make it easier to interpret
    > the complementarity. Also, when complementarity conditions are
    > used, the ordering of variables appearing in the complementarities
    > must match the declaration order (more in section Y).

- definition blocks

  > Definition blocks differ from condition blocks in that they define a
  > group of variables (for instance `states` or `expectations`) as a function of the
  > right hand side.

The types of variables appearing on the right hand side depend on the
block type. The variables enumerated on the left hand-side must appear
in the declaration order.

### Definitions

It is possible to define special variables that are substituted everywhere.

> **note**
> In the RBC model, output $y_t$, consumption $c_t$, return on investment $rk_t$ and wages $w_t$ are defined
recursively by:
> ``` {.sourceCode .yaml}
> definitions:
>    y: exp(z)*k^alpha*n^(1-alpha)
>    c: y - i
>    rk: alpha*y/k
>    w: (1-alpha)*y/n
> ```
> These variables, can then be used in equations, with any timing. For instance in the
> Euler equation `c(1)` which represents $c_{t+1}$ will implicitly be replaced by:
> `exp(z(1))*k(1)^alpha*n(1)^(1-alpha) - i(1)` (parameter `alpha` is not shifted).

### Calibration section

The role of the calibration section consists in providing values for the
parameters and the variables. The calibration of all parameters
appearing in the equation is of course strictly necessary while the
calibration of other types of variables is useful to define the
steady-state or an initial guess of the steady-state.

The calibrated values are also substituted in other sections, including
the `domain:`, `exogenous:` and `options:` sections. This is particularly useful to make the
covariance matrix depend on model parameters, or to adapt the
state-space to the model's calibration.

The calibration is given by an associative dictionary mapping symbols to
define with values. The values can be either a scalar or an expression.
All symbols are treated in the same way, and values can depend upon each
other as long as there is a way to resolve them recursively.

In particular, it is possible to define a parameter in order to target a
special value of an endogenous variable at the steady-state. This is
done in the RBC example where steady-state labor is targeted with
`n: 0.33` and the parameter `phi` calibrated so that the optimal labor
supply equation holds at the steady-state (`chi: w/c^sigma/n^eta`).

All symbols that are defined in the symbols section but do not appear in
the calibration section are initialized with the value nan without
issuing any warning.

> **note**
>
> No clear policy has been established yet about how to deal with
> undeclared symbols in the calibration section. Avoid them.

### Exogenous Shocks

Exogenous shock processes are specified in the section `exogenous` . Dolo accepts various exogenous processes such as normally distributed iid shocks, VAR1 processes, and Markov Chain processes. Dolo also allows for specific types of Markov Chains such as Poisson Processes, Aging Processes, and Death Processes.

Here are examples of how to define different processes. Note the use of yaml tags.
Normal Shock: It has mean zero and variance `Sigma`.
```{.sourceCode .yaml}
exogenous:!Normal
   Sigma: [[0.016^2]]
```
VAR1: `rho` is the persistence (only one allowed for now). `Sigma` is the covariance matrix. (Note: if a scalar `sigma` is given, it will be converted to `[[sigma]]` to indicate that it is the variance (not covariance) of the process).  
```{.sourceCode .yaml}
exogenous:!VAR1
    rho: 0.9
    sigma: [[0.01, 0.001],
             [0.001, 0.02]]
    N: 3
```
Markov Chain:
```{.sourceCode .yaml}
exogenous: !MarkovChain
  values: [[-0.01],[0.01]]
  transitions: [[0.9, 0.1], [0.1, 0.9]]
```
Poisson Process: `K` is the number of nodes and `mu` is the probability of a new arrival.
```{.sourceCode .yaml}
exogenous: !PoissonProcess
  mu: 0.05
  K: 10
```

Aging Process: `mu` is the probability of death and `K` is the maximum age. Note this also encompasses an indicator for death, so in the definition of exogenous variables you will need a variable for age and an indicator for death.

```{.sourceCode .yaml}
exogenous:!AgingProcess
  mu: 0.02
  K: 8
```

Death Process: `mu` is the probability of dying.
```{.sourceCode .yaml}
exogenous:!DeathProcess
  mu:0.02
```
We can also specify more than one process. For instance if we want to combine a VAR1 and an Aging Process we use the tag `Product` and write:
```{.sourceCode .yaml}
exogenous: !Product
    p1: !VAR1
         rho: 0.75
         Sigma: [[0.015^2]]

         N: 3

    p2: !AgingProcess
        mu: 0.02
        K: 8

```

### Domain
The domain section defines lower and upper bounds for the exogenous and endogenous states. For example, in the RBC model, we write:
``` {.sourceCode .yaml}
domain:
  z: [-2*sig_z/(1-rho^2)^0.5,  2*sig_z/(1-rho^2)^0.5]
  k: [ k*0.5, k*1.5]
```
The part for `z` sets the bounds for the productivity process to be two times its asymptotic standard deviation.

The boundaries for capital are a 50% bracket around its steady-state level.

### Options
The options sections contains extra information needed to solve the model.The
section follows the mini-language convention, with all calibrated values
replaced by scalars and all keywords allowed.

Here we can define the grids and specify their type, Cartesian. We can also specify how many grid points we want for `z` and `k`. Here we choose 5 points for `z` and 50 points for `k`. Note that the grid points are listed in accordance with the declaration order of the variables.
``` {.sourceCode .yaml}
options:
  grid: !Cartesian
        orders: [5, 50]
```
