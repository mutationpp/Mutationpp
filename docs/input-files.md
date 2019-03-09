<a id="top"></a>
# Input Files

**Contents** <br>

- [Introduction](#introduction)
  - [File Location](#file-location)
  - [Simplified XML Language](#xml-language)
- [Mixture](#mixtures)
  - [Mixture Options](#mixture-options)
  - [Species List Descriptor](#species-list-descriptor)
  - [Named Elemental Compositions](#elemental-compositions)
- [Elements](#elements)
- [Thermodynamic Databases](#thermodynamic-databases)
  - [NASA 7-Coefficient Polynomials](#nasa-7)
  - [NASA 9-Coefficient Polynomials](#nasa-9)
  - [Rigid-Rotor / Harmonic-Oscillator](#rrho)
- [Reaction Mechanisms](#reaction-mechanisms)
  - [Reaction Definitions](#reaction-definitions)
  - [Unit Specifiers](#unit-specifiers)
  - [Example Mechanism](#example-mechanism)
- [Collision Integral Database](#collision-integrals) (work in progress)


## Introduction
<a id="top"></a>

The data files distributed with Mutation++ are located in one of the subdirectories under the __data__ directory as follows:

- __data__
  + __mechanisms__ chemical [reaction mechanisms](#reaction-mechanisms)
  + __mixtures__ [mixture definitions](#mixtures)
  + __thermo__ [elemental](#elements) and species [thermodynamic databases](#thermodynamic-databases)
  + __transfer__ internal energy [transfer model databases](#transfer-databases)
  + __transport__ @subpage collisions "collision integral database"

### File location
<a id="file-location"></a>

Mutation++ is aware of its `data_directory`. It defaults to the `MPP_DATA_DIRECTORY` environment variable, which should have been setup during the [installation](@ref installation) of Mutation++ and usually points to the aforementioned data directory.

In addition, Mutation++ is aware of the `working_directory`, where an executable that utilizes Mutation++ is being run. [Mixture files](@ref mixtures) may also be in the working directory.

When a data file `name`.`ext` is requested, Mutation++ searches for it looking up locations in the following order:
1. `working_directory`/`name`.`ext`
2. `working_directory`/`dir`/`name`.`ext`
3. `data_directory`/`name`.`ext`
4. `data_directory`/`dir`/`name`.`ext`

where `dir` is any subdirectory. In practice, it allows a user to supersede a default data file with their own, by simply locating it in the working directory.

### Simplified XML Language
<a id="xml-language"></a>

Many of these files are written in a simplified version of the Extensible Markup
Language (XML).  XML provides a human readable, yet complex and extensible format for data
to be stored with only a few, limited rules.  An example XML fragment is shown below.

```xml
<!-- Comment string -->
<root_tag attribute="value">
    <child1_tag>
        text
    </child1_tag>
    <child2_tag attribute="another value" />
</root_tag>
```

An XML document begins with a _root element_.  Every _element_ (or _node_) must begin
with a _tag_ that identifies what type of element it is.  The root element depicted in
the example is of type `root_tag`.  Every element also ends with an _end-tag_
which signifies the end of the element.  Each element may have as many _attribute_ / _value
pairs_ as is desired following the element's tag.  Each pair must consist of an attribute
name followed by an equal sign and the value of the attribute in quotations.

Between the tag and end-tag of an element, an element may also contain one or more _child elements_ or _text_ (but not both).  From the figure, the root element contains two child elements named `child1_tag` and `child2_tag`.  Note that the first child element is an example of an element which contains text instead of more child elements.  The second child element is an example of a short-hand format for elements which only contain attributes.  For such elements, a full end-tag is not necessary.  Instead, simply putting `/>` after the attribute list is sufficient to end the element.  Finally, _comments_ can be inserted anywhere outside of element tags.  Comment strings begin with `<!--` and end with `-->` and can be spread over multiple lines.


## Mixtures 
<a id="mixtures"></a>

Mixture files are located in the `data/mixtures` directory.  They are the primary
input mode in Mutation++ because they provide the list of species to be loaded as
well as any options that can be used to control the behavior of Mutation++.

```xml
<!-- Example mixture file -->
<mixture thermo_db="RRHO">
    <!-- Species list -->
    <species>
        N2 N2+ N N+ e-
    </species>
</mixture>
```

### Mixture Options
<a id="mixture-options"></a>

The following options are available as attributes in the `mixture` element.  If an attribute is not given, then the __bold value__ is taken as the __default__.

Attribute              | Possible Values                                     | Description
-----------------------|-----------------------------------------------------|------------
`mechanism`            | __none__, name                                      | name of [reaction mechanism](#reaction_mechanisms)
`thermal_conductivity` | `CG`, __LDLT__, `Wilke`                             | choice of heavy particle translational thermal conductivity algorithm
`thermo_db`            | __RRHO__, `NASA-7`, `NASA-9`                        | choice of [thermodynamic database](#thermodynamic_databases)
`state_model`          | __ChemNonEq1T__, `ChemNonEqTTv`, `Equil`, `EquilTP` | choice of [state model](#statemodels)
`use_transport`        | `no`, __yes__                                       | whether or not to load transport data
`viscosity`            | `CG`, `Gupta-Yos`, __LDLT__, `Wilke`                | choice of viscosity algorithm

### Species List Descriptor
<a id="species-list-descriptor"></a>

A _species list descriptor_ tells Mutation++ which species to load from the chosen
thermodynamic database.  The simplest version of a descriptor is a list of species
names separated by white space as in the example mixture above.

A species list can also contain a more generic descriptor which allows the user to
choose species based on a given set of rules.  For now, the rule simply allows
to choose all species in a particular category which contian certain elements from a
list.  The format for this descriptor is

```
{ category with element_list }
```

where `category` can be one of the following strings:
- `gases` selects from all gases in the database
- `liquids` selects from all liquids in the database
- `solids` selects from all solids in the database
- `condensed` combines `liquids` and `solids`
- __all__ selects from all species in the database

The `element_list` should be a comma-separated list of [element names](#elements).
A rule can also be combined with a list of species names.  For example, if you
want a mixture with all gas species containing the elements C,H,O and graphite, the
species node would look like the following:

```xml
<species>
    {gases with C,H,O,e-} C(gr)
<species>
```

Note that the electron is also included in the above example.  This allows ions in
the mixture because they are treated as species with the fake "electron element"
(see the [Elements](#elements) section for more details).

**Species Names** <br>

Species names given in any data file, must obey the following rules:
- the characters `"`, `{`, `}`, `=`, `<`, `>`, and spaces are __not allowed__
- the first character __cannot be__ a numeric digit
- placing `(#)` (where `#` represents a whole number) at the end of a name indicates the
  `#`<sup>th</sup> electronic energy level of that species.
  + example: `N2(3)` is the 4<sup>th</sup> electronic level of N<sub>2</sub> because
    the index begins at 0
  + note that this only works if the thermodynamic database being used includes
    electronic energy levels for that species

**Species Order** <br>
<a id="species-order"></a>

A list of species __may not be__ represented internally in Mutation++ in the same order they
are specified in the species list descriptor.  Mutation++ places the species in an order which makes indexing the species the most convenient.  In general, the following steps are taken to order the species:
1. if present, the electron is placed at the beginning of the list
2. gas species are placed in front of condensed phase species
3. if electronic states are specified for any species, they are expanded in place
4. finally, species take the order given in the [list descriptor](#species-list-descriptor)

A good way to see how the species are actually ordered in Mutation++ is to check the
mixture with the [checkmix](checkmix.md) program.

### Named Elemental Compositions
<a id="elemental-compositions"></a>
Named elemental compositions can be included in the mixture file which can then
be retrieved inside Mutation++.  Many of the tools included with Mutation++ can also use this information to simplify input on the command line.  An example of a list of element compositions for the 11-species Air mixture are shown below.

```xml
<element_compositions default="air1">
    <composition name="air1"> e-:0.0, N:0.79, O: 0.21 </composition>
    <composition name="air2"> e-:0.0, N:0.80, O: 0.20 </composition>
</element_compositions>
```

Note that if the `default` attribute is used as above, the composition with the
name in the `default` value will be use as the default composition when computing
equilibrium calculations.  If no default is specified, then the first composition is used by default.


## Elements
<a id="elements"></a>

Element names and molecular weights are stored in the `elements.xml` file in the
`data/thermo` directory.  Each element is given as an XML node in the file.  An
example of the Nitrogen XML node is given as

```xml
<!-- Nitrogen -->
<element name="N">
    <mw units="g/mol">14.0067</mw>
</element>
```

__Note that the electron is treated like an element in Mutation++ and should not be
changed__ in the `elements.xml` file.  In general, most of the known elements are
already given in `elements.xml` and it is unlikely that it should need to be
modified.


## Thermodynamic Databases
<a id="thermodynamic-databases"></a>

### NASA 7-Coefficient Polynomials
<a id="nasa-7"></a>

Constants                                                        | Format       | Columns
-----------------------------------------------------------------|--------------|--------
<b>Line 1</b>                                                    | -            | -
Species name                                                     | `A18`        | `1-18`
Date code                                                        | `A6`         | `19-24`
Atomic symbols and formula                                       | `4(A2,I3)`   | `25-44`
Phase of species (S, L, G for solid, liquid, gas)                | `A1`         | `45`
Low temperature                                                  | `F10.0`      | `46-55`
High temperature                                                 | `F10.0`      | `56-65`
Common temperature (blank for default of 1000 K)                 | `F8.0`       | `66-73`
Atomic symbols and formula (blank if not needed)                 | `A2,I3`      | `74-78`
The integer '1'                                                  | `I1`         | `80`
<b>Line 2</b>                                                    | -            | -
Coefficients a_1-a_5 for upper temperature range           | `5(F15.8)`   | `1-75`
The integer '2'                                                  | `I1`         | `80`
<b>Line 3</b>                                                    | -            | -
Coefficients b_1 and b_2 for upper temperature range | `2(F15.8)`   | `1-30`
Coefficients  a_1-a_3 for lower temperature range          | `3(F15.8)`   | `31-75`
The integer '3'                                                  | `I1`         | `80`
<b>Line 4</b>                                                    | -            | -
Coefficients a_4-a_5 for lower temperature range           | `2(F15.8)`   | `1-30`
Coefficients b_1 and b_2 for lower temperature range | `2(F15.8)`   | `31-60`
The integer '4'                                                  | `I1`         | `80`


```
N                 L 6/88N   1    0    0    0G   200.000  6000.000 1000.        1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226244E-10-0.20360983E-14    2
 0.56133775E+05 0.46496095E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104638E+05 0.41939088E+01 0.56850013E+05    4
```

### NASA 9-Coefficient Polynomials
<a id="nasa-9"></a>

[McBride and Sanford (1996)](bibliography.md#McBride1996)

Constants                                               | Format       | Columns
--------------------------------------------------------|--------------|--------
<b>Line 1</b>                                           | -            | -
Species name                                            | `A24`        | `1-24`
Comments (data source)                                  | `A56`        | `25-80`
<b>Line 2</b>                                           | -            | -
Number of T intervals                             | `I2`         | `2`
Optional identification code                            | `A6`         | `4-9`
Chemical formulas, symbols, and numbers                 | `5(A2,F6.2)` | `11-50`
Zero for gas, nonzero for condensed phases              | `I1`         | `52`
Molecular weight                                        | `F13.5`      | `53-65`
Heat of formation at 298.15 K, J/mol                    | `F13.5`      | `66-80`
<b>Line 3</b>                                           | -            | -
Temperature range                                       | `2F10.3`     | `2-21`
Number of coefficients for C_p^\circ/R_u          | `I1`         | `23`
T exponents in polynomial for C_p^\circ/R_u | `8F5.1`      | `24-63`
H^\circ (298.15)-H^\circ (0), J/mol               | `F15.3`      | `66-80`
<b>Line 4</b>                                           | -            | -
First five coefficients for C_p^\circ/R_u         | `5F16.8`     | `1-80`
<b>Line 5</b>                                           | -            | -
Last three coefficients for C_p^\circ/R_u         | `3F16.8`     | `1-48`
Integration constants  b_1  and  b_2        | `2F16.8`     | `49-80`
<i>Repeat 3, 4, and 5 for each interval...</i>          | -            | -


```
N2                Ref-Elm. Gurvich,1978 pt1 p280 pt2 p207.
 3 tpis78 N   2.00    0.00    0.00    0.00    0.00 0   28.0134000          0.000
    200.000   1000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8670.104
 2.210371497D+04-3.818461820D+02 6.082738360D+00-8.530914410D-03 1.384646189D-05
-9.625793620D-09 2.519705809D-12                 7.108460860D+02-1.076003744D+01
   1000.000   6000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8670.104
 5.877124060D+05-2.239249073D+03 6.066949220D+00-6.139685500D-04 1.491806679D-07
-1.923105485D-11 1.061954386D-15                 1.283210415D+04-1.586640027D+01
   6000.000  20000.0007 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0         8670.104
 8.310139160D+08-6.420733540D+05 2.020264635D+02-3.065092046D-02 2.486903333D-06
-9.705954110D-11 1.437538881D-15                 4.938707040D+06-1.672099740D+03
```

### Rigid-Rotator Harmonic-Oscillators
<a id="rrho"></a>

```xml
<!-- Nitrogen diatomic -->
<species name="N2">
	<stoichiometry>
		N: 2
	</stoichiometry>
	<thermodynamics type="RRHO">
		<linear>
			yes
		</linear>
		<rotational_temperature units="K">
			2.886
		</rotational_temperature>
		<steric_factor>
			2
		</steric_factor>
		<vibrational_temperatures units="K">
			3408.464
		</vibrational_temperatures>
		<electronic_levels units="1/cm">
			<level degeneracy="1" energy="0.0" />
			<level degeneracy="3" energy="50203.5" />
			<level degeneracy="6" energy="59619.2" />
			<level degeneracy="6" energy="59808.0" />
			<level degeneracy="3" energy="66272.5" />
			<level degeneracy="1" energy="68152.7" />
			<level degeneracy="2" energy="69283.0" />
			<level degeneracy="2" energy="72097.6" />
			<level degeneracy="5" energy="77600.0" />
			<level degeneracy="1" energy="85200.0" />
			<level degeneracy="6" energy="86800.0" />
			<level degeneracy="6" energy="89136.7" />
			<level degeneracy="10" energy="93000.0" />
			<level degeneracy="6" energy="97603.0" />
			<level degeneracy="6" energy="104600.0" />
		</electronic_levels>
		<formation_enthalpy T="298 K" P="1 atm" units="J/mol">
			0.0
		</formation_enthalpy>
	</thermodynamics>
</species>
```


## Reaction Mechanisms
<a id="reaction-mechanisms"></a>

A mechanism name referenced in a mixture file will correspond to a file with an .xml
extension in the mechanisms directory.  A mechanism file describes a reaction mechanism
which is a set of elementary reactions.

Mechanism files begin with a root element called `mechanism`.  A `name` attribute may be
given but it is ignored.  The name of the mechanism is determined by the filename (without
the `.xml` extension).

The root `mechanism` element may have any number of child elements which correspond to
__reaction elements__ or __unit specifiers__.

### Reaction Definitions
<a id="reaction-definitions"></a>

Below is an example `reaction` element
```xml
<reaction formula="N2+M=2N+M">
    <arrhenius A="7.0E+21" n="-1.6" T="113200." />
    <M>N2:3,N:4</M>
</reaction>
```

`reaction` elements must contain a `formula` attribute which specifies the reaction
formula for the given reaction.  The `formula` string must obey the following rules:
- species names are separated by the `+` symbol
- an `M` in the reactant and product lists indicate a thirdbody reaction (note that this
  character is reserved and cannot be a species name)
- stoichiometric coefficients must directly precede a species name (if the coefficient is
  a `1`, then it can be left out)
- reactants and products must be separated by one of the following
  + `=` indicates a reversible reaction
  + `=>` indicates an irreversible reaction
- __a maximum of 3 reactants and products are allowed__

If the reaction is a thirdbody reaction, thirdbody efficiency factors can be specified in
a `M` child node of the `reaction` element.  The `M` text should be a list of of
colon-separated pairs of species names and efficiency factor, each separated by commas
such as in the example above.  __Any efficiency factors which are not specified default to__
__a value of 1.  Free electrons do not participate in thirdbody reactions.__

**Rate Laws** <br>

A `reaction` element __must contain a child node which describes which rate law to use__
when evaluating the forward rate coefficient.  The rate law node is specified by the tag
and its corresponding attributes correspond to the parameters for that given rate law.  The
following rate rate laws are currently supported.

- `arrhenius`  k(T) = A T<sup>n</sup> exp[-Ea / (Ru T)]

Att. | Value
-----|--------
`A`  | pre-exponential factor
`n`  | temperature exponent
`Ea` | activation energy
`T`  | characteristic temperature Ea / Ru

_note: only one of `Ea` or `T` may be used, not both_

### Unit Specifiers
<a id="unit-specifiers"></a>

A unit specifier element is used to specify the units of the parameters in each rate law
if needed.  A unit specifier element is denoted by a tag which has `_units` appended to
the rate law tag (ie: `arrhenius_units`).  Unit specifiers can be placed anywhere in a
reaction mechanism and will apply to all rate laws of that type for reactions listed after
the unit specifier.  The attributes of the unit specifier element correspond to the
parameters of the rate law.  Some parameters need units specified for several physical
quantities (such as mass, length, time, etc.).  These will be separated by commas.

- `arrhenius_units`

Att. | Value                               | Description
-----|-------------------------------------|-------------
`A`  | quantity, length, time, temperature | units of pre-exponential factor
`E`  | energy, quantity, temperature       | units of activation energy and characteristic temperature

### Example Mechanism
<a id="example-mechanism"></a>

As an example of how to create a mechanism file, consider the following example reaction
mechanism for a 5-species Nitrogen mixture of N2, N2+, N, N+ and e- with each reaction
controlled by an Arrhenius rate law.

\# | Formula                                  | A [mol,cm,s,K] | n     | Ea [K]
---|------------------------------------------|----------------|-------|--------
1  | N2 + N_2 <-> 2N + N2    | 1.0e21         | -1.6  | 113200
2  | N2 + N <-> 2N + N       | 3.0e21         | -1.6  | 113200
3  | N2 + N^+ <-> 2N + N+    | 1.0e21         | -1.6  | 113200
4  | N2 + e- <-> 2N + e-     | 7.0e22         | -1.6  | 113200
5  | N + e- <-> N+ + e- + e- | 2.5e30         | -3.82 | 168200
6  | N + N <-> N2+ + e-      | 4.4e7          | 1.5   | 67500

Reactions 1-4 in the above example have the same rate constants except for the
pre-exponential factor differ only by the thirdbody species, making them good candidates
to combine into a generic thirdbody reaction with efficiency factors.  However, free
electrons cannot be thirdbodies in Mutation++, so reaction 4 must be treated separately.
Note also that N2+ is not a possible thirdbody in reactions 1-4, thus its efficiency
factor should be zero.  The remaining reactions, 5-6, should be specified normally.  The
resulting mechanism file `example.xml` is given below.

```xml
<!-- Example mechanism for 5-species N2 mix -->
<mechanism name="example">

    <arrhenius_units A="mol,cm,s,K" E="kcal,mol,K" />

    <!-- # 1-3 -->
    <reaction formula="N2+M=2N+M">
        <arrhenius A="1.0E+21" n="-1.6" T="113200." />
        <M> N:3, N2+:0 </M>
    </reaction>

    <!-- # 4 -->
    <reaction formula="N2+e-=2N+e-">
        <arrhenius A="7.0E+22" n="-1.6" T="113200." />
    </reaction>

    <!-- # 5 -->
    <reaction formula="N+e-=N++e-+e-">
        <arrhenius A="2.5E+30" n="-3.82" T="168200." />
    </reaction>

    <!-- # 6 -->
    <reaction formula="N+N=N2++e-">
        <arrhenius A="4.4E+07" n="+1.50" T="67500." />
    </reaction>

</mechanism>
```

## Collision Integral Database
<a id="collision-integrals"></a>

_This is a workin in progress._

@section collisions_intro Introduction
The collision integral database is provided in the `data/transport/collisions.xml` file.
This page documents how collision integral data is stored and configured.  In particular,
- how to define [collision pairs](\ref collisions_pairs),
- how to use different [collision integral models](\ref collisions_types),
- how to tune the [default behavior](\ref collisions_defaults),
- how to specify [global options](\ref collisions_global_options),
- and finally, how to provide [species-specific data](\ref collisions_species_data) needed
  for some collision integral models.

All collision integral databases must be contained within a `collisions` XML element as
follows.
\code{xml}
<collisions>
    <!-- All data corresponding to collision pairs and integrals go here -->
</collisions>
\endcode

----
@section collisions_pairs Defining collision pairs
The most basic node type which can be given as a child of `collisions` is a `pair`
node, which specifies an individual collision pair.
\code{xml}
<pair s1="first species name" s2="second species name">
    <!-- Specify models for each collision integral kind here -->
</pair>
\endcode
The `pair` node has two attributes, `s1` and `s2`, which are the names of the two
interacting species in the pair.  Child nodes can then be given which provide the
data necessary to compute any of the collision integrals for that pair.  The format
of the collision integral nodes is described next.

----
@section collisions_types Collision integral types
Specific collision integral data are provided as XML nodes with the node tag as the 
name of the collision integral kind.  For example, to specify the 
\f$\overline{Q}^{(1,1)}\f$ integral for the N-O interaction pair, use
\code{xml}
<pair s1="N" s2="O">
    <Q11 type="model" />
</pair>
\endcode
where the `type` attribute is the collision integral model to use.  The attributes and
possible child nodes or text needed in the collision integral XML node depends on the
type of model being used.
The following table summarizes the various collision integral types which are
available.  Additional information, including their respective XML format, can
be found in the following subsections.

type              | Description
------------------|------------
`Bruno-Eq(11)`    | Implements the curve-fit found in Eq. (11) of \cite Bruno2010
`Bruno-Eq(17)`    | Implements the curve-fit found in Eq. (17) of \cite Bruno2010
`Bruno-Eq(19)`    | Implements the curve-fit found in Eq. (19) of \cite Bruno2010
`constant`        | Constant value
`Debye-Huckel`    | Based on Debye-Huckel potential (Coulomb potential, screened by Debye length)
`exp-poly`        | Exponential polynomial
`from A*`         | Uses A* relationship to compute missing integral
`from B*`         | Uses B* relationship to compute missing integral
`from C*`         | Uses C* relationship to compute missing integral
`Langevin`        | Based on the Langevin potential
`Murphy`          | Implements Eq. (5) of \cite Murphy1995 for combining elastic and inelastic interactions
`Pirani`          | Implements the Pirani potential, using the approximate curve-fits of Laricchiuta et al. \cite Laricchiuta2007
`ratio`           | Ratio of another integral
`table`           | Tabulated data
`warning`         | Prints a warning if used

Regardless of the `type`, all collision integral nodes can specify the following 
attributes.

Attribute  | Description
-----------|------------
`accuracy` | The estimated accuracy of the integral in percent off true value
`multpi`   | (`yes` or `no`) Should the integral be multipled by \f$\pi\f$?
`ref`      | A reference for the data
`units`    | The units of the integral provided by the data

Note that while theoretically any name can be used as a tag for specifying collision
integrals, the following names are reserved for use within the various transport 
algorithms: `Q11`, `Q12`, `Q13`, `Q14`, `Q15`, `Q22`, `Q23`, `Q24`, `Ast`, `Bst`, and
`Cst`.

The following code snippet taken from the `collisions.xml` file provided with the library
shows the specification of the collision integrals needed for the N2-N2 interaction.
\code{xml}
<pair s1="N2" s2="N2">
    <Q11 type="table" units="K,Å-Å" multpi="yes" ref="Wright2005" accuracy="10">
            300   600  1000  2000  4000  6000  8000 10000,
          12.23 10.60  9.79  8.60  7.49  6.87  6.43  6.06 </Q11>
    <Q22 type="table" units="K,Å-Å" multpi="yes" ref="Wright2005" accuracy="10">
            300   600  1000  2000  4000  6000  8000 10000,
          13.72 11.80 10.94  9.82  8.70  8.08  7.58  7.32 </Q22>
</pair>
\endcode

Note that the values of B* and C* which are also need for the neutral-neutral interaction
are not shown because they are specified through [default behaviour](\ref collisions_defaults).

@subsection collision_types_bruno_11 Bruno-Eq(11)
Implements the curve-fit found in Eq. (11) of Bruno et al. \cite Bruno2010,
\f[
\sigma^2\Omega^{(l,s)*} = 
    [a_1 + a_2 x] \frac{\exp[(x-a_3)/a_4]}{\exp[(x-a_3)/a_4] + \exp[(a_3-x)/a_4]} +
    a_5 \frac{\exp[(x-a_6)/a_7]}{\exp[(x-a_6)/a_7] + \exp[(a_6-x)/a_7]},
\f]
where \f$a_i\f$ are given coefficients and \f$x = \ln(T)\f$.  The following code
snippet details how the 7 coefficients are provided in XML format.
\code{xml}
<Q11 type="Bruno-Eq(11)">
    15.09506044  -1.25710008   9.57839369  -3.80371463   
     0.98646613   9.25705877  -0.93611707 
</Q11>
\endcode
        
@subsection collision_types_bruno_17 Bruno-Eq(17)
Implements the curve-fit found in Eq. (17) of Bruno et al. \cite Bruno2010 for
modeling charge exchange interactions,
\f[
\sigma^2\Omega^{(l,s)*} = d_1 + d_2 x + d_3 x^2,
\f]
where \f$d_i\f$ are given coefficients and \f$x = \ln(T)\f$.  The following code
snippet details how the 3 coefficients are provided in XML format.
\code{xml}
<Q11 type="Bruno-Eq(17)"> 6.3544e+01 -5.0093e+00  9.8797e-02 </Q11>
\endcode

@subsection collision_types_bruno_19 Bruno-Eq(19)
Implements the curve-fit found in Eq. (19) of Bruno et al. \cite Bruno2010 for
modeling electron-neutral interactions,
\f[
\sigma^2\Omega^{(l,s)*} = 
    g_3 x^{g_5} \frac{\exp[(x-g_1)/a_2]}{\exp[(x-g_1)/g_2] + \exp[(g_1-x)/g_2]} +
    g_6 \exp\big[-\big(\frac{x-g_7}{g_8}\big)^2\big] + g_4,
\f]
where \f$g_i\f$ are given coefficients and \f$x = \ln(T)\f$.  The following code
snippet details how the 8 coefficients are provided in XML format.
\code{xml}
<Q11 type="Bruno-Eq(19)">
    1.035291134e+1 -1.583011620e+0 1.245844036e+1 -2.328519000e-1 
    5.366285730e-2 -5.343729290e+0 9.355617520e+0 -2.154634270e+0 
</Q11>
\endcode

@subsection collision_types_constant constant
Uses a constant value for the collision integral.
\code{xml}
<Q11 type="constant" value="1.0e-20"/>
\endcode

@subsection collision_types_debye_huckel Debye-Huckel
Models interactions between charged particles using the screened Coulomb potential
shielded by the Debye length, known as the Debye-Huckle potential.  The tables of
precomputed collision integrals provided by Mason, Munn, and Smith \cite Mason1967
and Devoto \cite Devoto1973 are used to implement the integrals.  The following
code snippet shows how to use this model.
\code{xml}
<Q11 type="Debye-Huckel"/>
\endcode
Note the actual value of the integral will depend on the kind of collision integral
(or ratio) such as `Q11` or `Ast` and whether or not the interaction is repulsive
or attractive.

@subsection collision_types_exp_poly exp-poly
Implements an exponential polynomial curve-fit expression of the form
\f[
\overline{Q}^{(l,s)}_{i,j}(T) = \exp\big(\sum_{i=0}^{n-1} a_i x^i \big),
\f]
where \f$a_i\f$ are given coefficients and \f$x = \ln(T)\f$.  The following XML
snippet shows how to format the coefficients in the database.
\code{xml}
<Q11 type="exp-poly">
    -4.942712e-03 1.113337e-01 -9.844008e-01 6.435295e+00
</Q11>
\endcode
Note that the **coefficients are given in reverse order** and any number of 
coefficients can be used.

@subsection collision_types_from_A from A*
Uses relationship \f$ A^{*}_{ij} = \overline{Q}^{(2,2)}_{i,j} / 
\overline{Q}^{(1,1)}_{i,j}\f$ to compute either \f$A^{*}_{ij}\f$, 
\f$\overline{Q}^{(1,1)}_{i,j}\f$, or \f$\overline{Q}^{(2,2)}_{i,j}\f$, given
the other two.  Simply use the `from A*` type for the missing integral, as in
\code{xml}
<Q11 type="from A*"/>
\endcode
and provide concrete types for the other two.

@subsection collision_types_from_B from B*
Uses relationship \f$ B^{*}_{ij} = (5\overline{Q}^{(1,2)}_{i,j} - 4
\overline{Q}^{(1,3)}_{i,j})/\overline{Q}^{(1,1)}_{i,j}\f$ to compute either 
\f$B^{*}_{ij}\f$, \f$\overline{Q}^{(1,1)}_{i,j}\f$, \f$\overline{Q}^{(1,2)}_{i,j}\f$,
or \f$\overline{Q}^{(1,3)}_{i,j}\f$, given the other three.  Simply use the `from B*` 
type for the missing integral, as in
\code{xml}
<Q11 type="from B*"/>
\endcode
and provide concrete types for the other three.

@subsection collision_types_from_C from C*
Uses relationship \f$ C^{*}_{ij} = \overline{Q}^{(1,2)}_{i,j} / 
\overline{Q}^{(1,1)}_{i,j}\f$ to compute either \f$C^{*}_{ij}\f$, 
\f$\overline{Q}^{(1,1)}_{i,j}\f$, or \f$\overline{Q}^{(1,2)}_{i,j}\f$, given
the other two.  Simply use the `from C*` type for the missing integral, as in
\code{xml}
<Q11 type="from C*"/>
\endcode
and provide concrete types for the other two.

@subsection collision_types_langevin Langevin
The Langevin (polarization) potential has been used extensively for modeling 
ion-neutral interactions and represents a special case of the inverse power
potential.  A closed form analytical solution of the collision integral using
this potential is given as
\f[
\overline{Q}^{(l,s)}_{i,j} = C^{(l,s)} z \pi \sqrt{\frac{\alpha}{T}},
\f]
where \f$C^{(l,s)}\f$ is a constant depending on \f$l\f$ and \f$s\f$, \f$z\f$ is
the elementary charge of the ion, and \f$\alpha\f$ is the dipole polarizability
of the neutral.  An example is given in the following code snippet.
\code{xml}
<Q11 type="Langevin">
\endcode
Note that the values of \f$l\f$ and \f$s\f$ are automatically determined from the
kind of integral (in this case `Q11`) and \f$z\f$ is determined from the interaction
pair that the integral is assigned to.  The dipole polarizability of the netural
must be provided in a [polarizability table](@ref collisions_polarizabilities) 
somewhere in the database.

@subsection collision_types_murphy Murphy
Implements Eq. (5) of \cite Murphy1995 for combining elastic and inelastic
interactions, where
\f[
\overline{Q}_{ij}^{(l,s)}(T) = \sqrt{Q_1(T)^2 + Q_2(T)^2}
\f]
where \f$Q_1\f$ and \f$Q_2\f$ are the elastic and
inelastic parts.  This combination rule is often used for ion-parent
interactions where charge exchange interactions become important.  The XML format is
as follows.
\code{xml}
<Q11 type="Murphy">
    <Q1 type="Bruno-Eq(11)">
        46.68783791  -0.33303803   4.25686770  -2.03851201  
        14.98170958   8.59618369  -1.65616736 </Q1>
    <Q2 type="Bruno-Eq(17)">
        6.3544e+01 -5.0093e+00  9.8797e-02 </Q2>
</Q11>
\endcode
Note that the types of `Q1` and `Q2` can be any valid collision integral type.

@subsection collision_types_pirani Pirani
The phenomenological Pirani potential has been introduced in \cite Pirani2004, 
taking the form
\f[
\phi = \epsilon_0 \big[ \frac{m}{n(x)-m} \big(\frac{1}{x}\big)^{n(x)} -
\frac{n(x)}{n(x)-m} \big(\frac{1}{x}\big)^{m} \big],
\f]
where \f$ x = r/r_e\f$ and \f$n(x) = \beta + 4x^2\f$.  For netural-neutral and
ion-neutral interactions, m has the value of 6 and 4, respectively.  The value
of \f$\beta\f$ ranges from 6 to 10 depending on the "hardness" of the the 
interacting electronic distribution densities and can be estimated as 
\f[
\beta = 6 + \frac{5}{s_1 + s_2},
\f]
where \f$s_1\f$ and \f$s_2\f$ are the softness values of the colliding partners
1 and 2.  This is defined as the cubic root of the dipole polarizabilities of 
the two species.  The remaining parameters, \f$r_e\f$ and \f$\epsilon_0\f$ must
be either fit to match experimental observations or estimated based on correlations.
Both options are available:
\code{xml}
<!-- Option 1: provide necessary parameters -->
<Q11 type="Pirani" beta="7.2644" eps0="0.00798" re="3.621" />

<!-- Option 2: estimate parameters
<Q11 type="Pirani" />
\endcode
Laricchiuta et al. \cite Laricchiuta2007 provide the necessary correlation
formulas for determining the missing parameters if option 2 is chosen.  The
correlations rely on needing the [dipole polarizabilities](@ref collisions_polarizabilities) 
of the interacting species and, for neutral-neutral interactions, the 
[effective electrons](\ref collisions_effective_electrons) contributing to the 
polarization of the neutral species.

Curve-fits of the collision integrals computed using the Pirani potential have
been provided in \cite Laricchiuta2007.  These are used to evaluate collision
integrals using this potential, automatically taking into account the type of
interaction and integral.

@subsection collision_types_ratio ratio
Use a ratio of another defined collision integral.
\code{xml}
<Q22 type="ratio" ratio="1.1" integral="Q11" />
\endcode
The XML syntax takes an integral name and the ratio to multiply that integral by. In this
case `Q22` will be computed as 1.1 times `Q11`.

@subsection collision_types_table table
Interpolates tabulated collision integral data.  The table is provided in XML as two
space delimited lists, separated by a comma.
\code{xml}
<Q11 type="table" units="K,Å-Å" interpolator="l=Linear" clip="yes">
    1000  2000  4000  5000  6000  8000 10000 15000 20000,
    0.12  0.19  0.43  0.55  0.66  0.89  1.12  1.68  2.23 
</Q11>
\endcode
The first list indicates the temperature points and the second the collision integral
values, with the units specified by the units attribute.  The `interpolator`
attribute specifies which interpolator is used to interpolate on the table.  By 
default, `Linear` is used if no interpolator is specified.  `Chebyshev` and
`MonotoneCubic` may also be used.  The `clip` attribute specifies whether or not to
clip the collision integral at the bounds of the table.  By default this is `yes`.

@subsection collision_types_warning warning
Uses a constant value for the integral and prints a warning message.  This can be
used when no data is available for a given species pair and kind of integral, without
crashing.
\code{xml}
<Q11 type="warning" value="1.0e-20"/>  
\endcode

----
@section collisions_defaults Specifying default behavior
Specifying collision integral data for every collision pair in your mixture may
sound like a daunting task.  Chances are that most of the integrals you
need are already provided in the database.  For the ones that aren't available
already, the database provides a facility to specify default behavior by
specifying the `defaults` node.
\code{xml}
<!-- Default collision integral data -->
<defaults>
    <!-- Neutral-Neutral interactions -->
    <neutral-neutral>
        <!-- Collision integral definitions -->
    </neutral-neutral>
    
    <!-- Ion-Neutral interactions -->
    <ion-neutral>
        <!-- Collision integral definitions -->
    </ion-neutral>
    
    <!-- Electron-Neutral interactions -->
    <electron-neutral>
        <!-- Collision integral definitions -->
    </electron-neutral>
    
    <!-- Charged interactions -->
    <charged>
        <!-- Collision integral definitions -->
    </charged>
</defaults>
\endcode

The `defaults` node allows the user to configure how to provide default
collision integrals for missing data, based on the four general types of 
interactions: neutral-neutral, ion-neutral, electron-neutral, and charged.
Whenever the user requests to load a specific integral, the library will first
search the database for the explicit collision pair/integral kind requested. If
the pair or integral type is not found explicitely, the library will then resort
to the default integral provided in `defaults`, according to the type of 
interaction specified by the pair.  This approach has several benefits:
- The default behavior is self-documenting and easily tunable without writing
new code or recompiling the library
- Complex expressions can be created to find the most accurate collision model
when data is missing
- No extra work is required to add new species pairs if the default behavior is
good enough.

This last point is particularly valid in the case of charged interactions, for 
which the screened Coulomb potential is nearly always a valid approximation.
Therefore, all charged interaction pairs can be supported easily by setting
\code{xml}
<defaults>
    <!-- Charged interactions -->
    <charged>
        <Q11 type="Debye-Huckel"/>
        <Q12 type="Debye-Huckel"/>
        <Q13 type="Debye-Huckel"/>
        <Q13 type="Debye-Huckel"/>
        <Q14 type="Debye-Huckel"/>
        <Q15 type="Debye-Huckel"/>
        <Q22 type="Debye-Huckel"/>
        <Q23 type="Debye-Huckel"/>
        <Q24 type="Debye-Huckel"/>
        <Ast type="Debye-Huckel"/>
        <Bst type="Debye-Huckel"/>
        <Cst type="Debye-Huckel"/>
    </charged>

    <!-- Other interactions -->
</defaults>
\endcode
Other cases can be envisaged as well.  For example, using the 
[Langevin](@ref collision_types_langevin) model as default for all ion-neutral 
interactions means that only the [dipole polarizability](\ref collisions_polarizabilities) 
of the neutral species needs to be specified to add a new ion-neutral pair.

----
@section collisions_global_options Global options
The `global-options` node specifies options that act on the whole collision
integral database.  For the moment, there are only two options that can be set.
\code{xml}
<global-options>
    <tabulate Tmin="100" Tmax="50000" dT="100" />
    <integral type="table" interpolator="Linear" clip="true"/>
</global-options>
\endcode

The `tabulate` option instructs the database manager to tabulate all collision
integrals after loading them.  The temperature grid used in the tabulation is
specified through the `Tmin`, `Tmax`, and `dT` attributes which are the minimum
and maximum temperatures, and the constant temperature spacing, respectively.

The second option shown, `integral`, is used to specify global options for a 
specific collision integral type.  In the example, the default `interoplator`
and `clip` attributes for all [table](\ref collision_types_table) collision 
integrals are set to `Linear` and `true` respectively.  Note that for the moment, 
only attributes of [table](\ref collision_types_table) can be specified in this
way.

----
@section collisions_species_data Species data
Some of the collision integral models detailed above require additional information
about individual species, such as their dipole polarizabilities.  These parameters are
stored in generic XML nodes as children of the root `collisions` element, with the 
following general structure.
\code{xml}
<species-property units="some units">
    <species name="some name" value="some value" ref="some reference" units="other units" />
    <species name="some name" value="some value" ref="some reference" />
    <species name="some name" value="some value" ref="some reference" />
</species-property>
\endcode
Here, `species-property` denotes the name of the spcies property stored in the node.
Default units can be given as shown in the `units` attribute of the `species-property`
node.  The property node then contains child nodes for each individual species for which
data is provided.  Those nodes provide the species `name`, the `value` of the property
for that species, an optional `ref` attribute for attributing a reference for the data,
and the possibility to override the default `units`.  

__A note for developers:__ Species properties stored in this
way can be easily found in any child class of Mutation::Transport::CollisionIntegral by 
calling the Mutation::Transport::CollisionIntegral::loadSpeciesParameter() function.

@subsection collisions_polarizabilities Dipole polarizabilities
Species dipole polarizabilities are used as inputs to some collision integral models,
including [Langevin](@ref collision_types_langevin) and 
[Pirani](@ref collision_types_pirani).  An example listing of dipole polarizabilities
is given below.
\code{xml}
<!-- Dipole Polarizabilities -->
<dipole-polarizabilities units="Å-Å-Å">
    <species name="N"    value="1.10"   ref="Wright2007"/>
    <species name="N2"   value="1.74"   ref="Wright2007"/>
    <species name="NO"   value="1.70"   ref="Wright2007"/>
    <species name="O"    value="0.80"   ref="Wright2007"/>
    <species name="O2"   value="1.58"   ref="Wright2007"/>
</dipole-polarizabilities>
\endcode

@subsection collisions_effective_electrons Number of effective electrons in polarization
The number of effective electrons in the polarization of a given species is used
in the [Pirani](@ref collision_types_pirani) collision integral model.  These
values are stored in the species database `effective-electrons`, following the
example below.
\code{xml}
<effective-electrons>
    <species name="C"    value="4.0"   />
    <species name="C2"   value="6.00"  />
    <species name="C2H"  value="8.11"  />
    <species name="C2H2" value="10.00" />
    <species name="C2H4" value="12.00" />
</effective-electrons>
\endcode
The value of the effective electrons for new species can be estimated using the
approach of Cambi et al. \cite Cambi1991.

