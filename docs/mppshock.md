<a id="top"></a>

# Mutation++ Post Shock Equilibrium (mppshock)
Given specific free stream conditions (P, V, T), mppshock computes the Rankine-Hugoniot
relations and the relaxed, equilibrium thermodynamic state and velocity. The
[Rankine-Hugoniot relations](https://en.wikipedia.org/wiki/Rankine–Hugoniot_conditions)
are an expression of the conservation of mass, momentum and energy and, since chemically
reacting mixtures are treated the chemical composition is assumed to be frozen across the
shock. In the post shock region the mixture relaxes and the local thermodynamic equilibrium
is reached (Equilibrated).  Two available numerical methods are available in mppshock for
calculating this post shock state, the bisection and the Newton-Raphson method. The first one
has increased robustness, but is significantly slower than the Newton-Raphson algorithm.

## Usage
```bash
Usage:   mppshock [OPTIONS] mixture
    -h, --help          prints this help message
    -P                  Free stream pressure range in Pa (default: 1 atm)
    -T                  Free stream temperature range in K (default: 300 K)
    -V                  Free stream velocity range in m/s (default: 10000 m/s)
    -m                  Solution algorithm; 0: Bisection - 1: Newton-Raphson (default)
Example: mppshock -P 1000 -V 10000 -T 300 -m 1 air_5
```

## Result for the example
```bash
                               Pressure (Pa)      Velocity (m/s)     Temperature (K)
    Free Stream:                      1000.0             10000.0               300.0
    Rankine-Hugoniot:               964161.6              1672.7             48382.6
    Equilibrated:                  1070709.9               751.5             12110.0
```

