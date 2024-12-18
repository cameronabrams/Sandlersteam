# Sandlersteam

> Digitized steam tables from Sandler's 5th ed.

Sandlersteam implements a python interface to the steam tables found in Appendix III of _Chemical, Biochemical, and Engineering Thermodynamics_ (5th edition) by Stan Sandler (Wiley, USA). It should be used for educational purposes only.

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

## Release History

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
