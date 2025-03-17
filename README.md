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

One can also import the `SteamTables` dictionary from the state `state` module and then generate LaTeX-compatible versions of either blocks in the superheated/subcooled steam tables or entire saturated steam stables, listed by temperature or pressure. For example:

```python
>>> from sandlersteam.state import SteamTables as st
>>> print(st['suph'].to_latex(P=1.0))  # generates latex for the 1 MPa block of the superheated steam table
```
```
\begin{minipage}{0.6\textwidth}
\footnotesize\vspace{5mm}
\begin{center}
$P$ = 1.0 MPa\\*[1ex]
\begin{tabular}{>{\raggedleft}p{8mm}@{}p{5mm}>{\raggedleft}p{4mm}@{}p{10mm}>{\raggedleft}p{10mm}@{}p{3mm}>{\raggedleft}p{10mm}@{}p{3mm}>{\raggedleft\arraybackslash}p{3mm}@{}p{8mm}}
\toprule
\multicolumn{2}{c}{$T$~($^\circ$C)} & \multicolumn{2}{c}{$\hat{V}$} & \multicolumn{2}{c}{$\hat{U}$} & \multicolumn{2}{c}{$\hat{H}$} & \multicolumn{2}{c}{$\hat{S}$}\\
\toprule
\midrule
179 & .91 & 0 & .19444 & 2583 & .6 & 2778 & .1 & 6 & .5865 \\
200 &  & 0 & .2060 & 2621 & .9 & 2827 & .9 & 6 & .6940 \\
250 &  & 0 & .2327 & 2709 & .9 & 2942 & .6 & 6 & .9247 \\
300 &  & 0 & .2579 & 2793 & .2 & 3051 & .2 & 7 & .1229 \\
350 &  & 0 & .2825 & 2875 & .2 & 3157 & .7 & 7 & .3011 \\
400 &  & 0 & .3066 & 2957 & .3 & 3263 & .9 & 7 & .4651 \\
500 &  & 0 & .3541 & 3124 & .4 & 3478 & .5 & 7 & .7622 \\
600 &  & 0 & .4011 & 3296 & .8 & 3697 & .9 & 8 & .0290 \\
700 &  & 0 & .4478 & 3475 & .3 & 3923 & .1 & 8 & .2731 \\
800 &  & 0 & .4943 & 3660 & .4 & 4154 & .7 & 8 & .4996 \\
900 &  & 0 & .5407 & 3852 & .2 & 4392 & .9 & 8 & .7118 \\
1000 &  & 0 & .5871 & 4050 & .5 & 4637 & .6 & 8 & .9119 \\
1100 &  & 0 & .6335 & 4255 & .1 & 4888 & .6 & 9 & .1017 \\
1200 &  & 0 & .6798 & 4465 & .6 & 5145 & .4 & 9 & .2822 \\
1300 &  & 0 & .7261 & 4681 & .3 & 5407 & .4 & 9 & .4543 \\
\bottomrule
\end{tabular}
\end{center}
\end{minipage}
```

This renders as

![a steam table block](https://github.com/cameronabrams/Sandlersteam/raw/main/stimage.png)
<!-- ![a steam table block](stimage.png "1 MPa superheated steam table block") -->

## Release History

* 0.4.2:
    * `RandomState` class introduced
* 0.3.3:
    * Included subcooled liquid data in the `Request` capability
* 0.2.1
    * bugfix:  set `left` and `right` parameters in all `numpy.interp()` to `numpy.nan` to override the pinning default for extrapolated values
* 0.2.0
    * Added the `Request` class for dynamically selecting and outputting (as LaTeX) steam table blocks requested by (for example) exam problems
* 0.1.9
    * bugfix: allows for specification of `x` as 0
* 0.1.8
    * Update readme
* 0.1.7
    * Updated pandas indexing for saturated table LaTeX printing
* 0.1.5
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
