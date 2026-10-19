# Dolo

![Build Status](https://github.com/EconForge/Dolo.jl/workflows/CI/badge.svg)


Documentation: [http://www.econforge.org/Dolo.jl/latest/](http://www.econforge.org/Dolo.jl/latest/)

Dolo is a tool to describe and solve economic models. It provides a simple classification scheme to describe many types of models, represents them as efficient Julia objects, and provides many reference solution algorithms to find the solution of these models under rational expectations.

Models are written directly in Julia (see [examples/ymodels](examples/ymodels)). To write them as simple YAML text files instead, use [DoloYAML.jl](https://github.com/EconForge/DoloYAML.jl), which imports them as models fully compatible with Dolo.jl.

Dolo understand several types of nonlinear models with occasionally binding constraints (with or without exogenous discrete shocks), as well as local pertubations models.

More docs coming soon...
