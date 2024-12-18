# Sandlersteam

> Digitized steam tables from Sandler's 5th ed.

Sandlersteam implements a python interface to the steam tables found in Appendix III of _Chemical, Biochemical, and Engineering Thermodynamics_ (5th edition) by Stan Sandler (Wiley, USA). It should be used for educational purposes only.

The interface operates similarly to the IAPWS steam tables (which you should use instead of these).

## Installation 

Sandlersteam is available via `pip`:

```sh
pip install sandlersteam
```

## Usage example

Below we create a `State` object to define a thermodynamic state for steam at 100 deg. C and 0.1 MPa:

```python
>>> from sandlersteam.state import State
>>> state1 = State(T=100.0,P=0.1)
>>> state1.h  # enthalpy in kJ/kg
2676.2
>>> state1.u  # internal energy in kJ/kg
2506.7
>>> state1.v  # volume in m3/kg
1.6958
>>> state1.s  # entropy in kJ/kg-K
7.3614
```

Specifying a state requires values for two independent state variables.  The state variables recognized by sandlersteam are:

* `T` temperature in C
* `P` pressure in MPa
* `u` specific internal energy in kJ/kg
* `v` specific volume in m<sup>3</sup>/kg
* `h` specific enthalpy in kJ/kg
* `s` specific entropy in kJ/kg
* `x` quality; mass fraction of vapor in a saturated vapor/liquid system (between 0 and 1)

Initializing a `State` instance with any two of these values set by keyword parameters results in the other
properties receiving values by (bi)linear interpolation.

When specifying quality, a `State` objects acquires `Liquid` and `Vapor` attributes that each hold intensive, saturated single-phase property values.  The property value attributes owned directly by the `State` object reflect the quality-weighted sum of the respective single-phase values:

```python
>>> s = State(T=100,x=0.5)
>>> s.P
0.10135
>>> s.v
0.836972
>>> s.Vapor.v
1.6729
>>> s.Liquid.v
0.001044
>>> 0.5*(s.Vapor.v+s.Liquid.v)
0.836972
```

## Release History

* 0.1.2
    * Updated interpolators
* 0.1.1
    * Updated pyproject.toml and README.md
* 0.1.0
    * Initial version

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/sandlersteam/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
