/**
@page input_files Input Files

@tableofcontents

@section input_intro Introduction

The data files distributed with Mutation++ are located in one of the subdirectories under the __data__
directory as follows:

- __data__
  + __mechanisms__ chemical [reaction mechanisms](@ref reaction_mechanisms)
  + __mixtures__ [mixture definitions](@ref mixtures)
  + __thermo__ [elemental](@ref elements) and species [thermodynamic databases](@ref thermodynamic_databases)
  + __transfer__ internal energy [transfer model databases](@ref transfer_databases)
  + __transport__ @subpage collisions "collision integral database"

@section input_location File location

Mutation++ is aware of its `data_directory`. It defaults to the `MPP_DATA_DIRECTORY` environment variable,
which should have been setup during the [installation](@ref installation) of Mutation++ and usually points
to the aforementioned data directory.

In addition, Mutation++ is aware of the `working_directory`, where an executable that utilizes Mutation++ is
being run. [Mixture files](@ref mixtures) may also be in the working directory.

When a data file `name`.`ext` is requested, Mutation++ searches for it looking up locations in the
following order:

1. `working_directory`/`name`.`ext`
2. `working_directory`/`dir`/`name`.`ext`
3. `data_directory`/`name`.`ext`
4. `data_directory`/`dir`/`name`.`ext`

where `dir` is any subdirectory. In practice, it allows a user to supersede a default data file with their
own, by simply locating it in the working directory.

@subsection simple_xml Simplified XML Language
Many of these files are written in a simplified version of the Extensible Markup
Language (XML).  XML provides a human readable, yet complex and extensible format for data
to be stored with only a few, limited rules.  An example XML fragment is shown below.

\code{.xml}
<!-- Comment string -->
<root_tag attribute="value">
    <child1_tag>
        text
    </child1_tag>
    <child2_tag attribute="another value" />
</root_tag>
\endcode

An XML document begins with a __root element__.  Every __element__ (or __node__) must begin
with a __tag__ that identifies what type of element it is.  The root element depicted in
the example is of type `root_tag`.  Every element also ends with an __end-tag__
which signifies the end of the element.  Each element may have as many __attribute/value
pairs__ as is desired following the element's tag.  Each pair must consist of an attribute
name followed by an equal sign and the value of the attribute in quotations.

Between the tag and end-tag of an element, an element may also contain one or more __child elements__
or __text__ (but not both).  From the figure, the root element contains two child elements
named `child_tag` and `child2_tag`.  Note that the first child element is an example of an
element which contains text instead of more child elements.  The second child element is
an example of a short-hand format for elements which only contain attributes.  For such
elements, a full end-tag is not necessary.  Instead, simply putting `/>` after the attribute
list is sufficient to end the element.  Finally, __comments__ can be inserted anywhere outside
of element tags.  Comment strings begin with `<!-``-` and end with `-``->` and can be spread
over multiple lines.

----
@section mixtures Mixtures

Mixture files are located in the `data/mixtures` directory.  They are the primary
input mode in Mutation++ because they provide the list of species to be loaded as
well as any options that can be used to control the behavior of Mutation++.

@code{.xml}
<!-- Example mixture file -->
<mixture thermo_db="RRHO">
    <!-- Species list -->
    <species>
        N2 N2+ N N+ e-
    </species>
</mixture>
@endcode

@subsection miture_opts Mixture Options

The following options are available as attributes in the `mixture` element.  If an
attribute is not given, then the __bold value__ is taken as the __default__.

Attribute              | Possible Values                                     | Description
-----------------------|-----------------------------------------------------|------------
`mechanism`            | __none__, name                                      | name of [reaction mechanism](@ref reaction_mechanisms)
`thermal_conductivity` | `CG`, __LDLT__, `Wilke`                             | choice of heavy particle translational thermal conductivity algorithm
`thermo_db`            | __RRHO__, `NASA-7`, `NASA-9`                        | choice of [thermodynamic database](@ref thermodynamic_databases)
`state_model`          | __ChemNonEq1T__, `ChemNonEqTTv`, `Equil`, `EquilTP` | choice of [state model](@ref statemodels)
`use_transport`        | `no`, __yes__                                       | whether or not to load transport data
`viscosity`            | `CG`, `Gupta-Yos`, __LDLT__, `Wilke`                | choice of viscosity algorithm

@subsection species_list Species List Descriptor

A _species list descriptor_ tells Mutation++ which species to load from the chosen
thermodynamic database.  The simplest version of a descriptor is a list of species
names separated by white space as in the example mixture above.

A species list can also contain a more generic descriptor which allows the user to
choose species based on a given set of rules.  For now, the rule simply allows
to choose all species in a particular category which contian certain elements from a
list.  The format for this descriptor is

     { category with element_list }

where `category` can be one of the following strings:
- __gases__ selects from all gases in the database
- __liquids__ selects from all liquids in the database
- __solids__ selects from all solids in the database
- __condensed__ combines __liquids__ and __solids__
- __all__ selects from all species in the database

The `element_list` should be a comma-separated list of [element names](@ref elements).
A rule can also be combined with a list of species names.  For example, if you
want a mixture with all gas species containing the elements C,H,O and graphite, the
species node would look like the following:

@code{.xml}
<species>
    {gases with C,H,O,e-} C(gr)
<species>
@endcode

Note that the electron is also included in the above example.  This allows ions in
the mixture because they are treated as species with the fake "electron element"
(see the [Elements](@ref elements) section for more details).

@subsubsection species_names Species Names

Species names given in any data file, must obey the following rules:
- the characters `"`, `{`, `}`, `=`, `<`, `>`, and spaces are __not allowed__
- the first character __cannot be__ a numeric digit
- placing `(#)` (where `#` represents a whole number) at the end of a name indicates the
  `#`<sup>th</sup> electronic energy level of that species.
  + example: `N2(3)` is the 4<sup>th</sup> electronic level of N<sub>2</sub> because
    the index begins at 0
  + note that this only works if the thermodynamic database being used includes
    electronic energy levels for that species

@subsubsection species_order Species Order

A list of species may not be represented internally in Mutation++ in the order they
are specified in the species list descriptor.  Mutation++ places
the species in an order which makes indexing the species the most convenient.  In
general, the following steps are taken to order the species:
-# if present, the electron is placed at the beginning of the list
-# gas species are placed in front of condensed phase species
-# if electronic states are specified for any species, they are expanded in place
-# finally, species take the order given in the [list descriptor](@ref species_list)

A good way to see how the species are actually ordered in Mutation++ is to check the
mixture with the [checkmix](@ref checkmix) program.

@subsection compositions Named Elemental Compositions

Named elemental compositions can be included in the mixture file which can then
be retrieved inside Mutation++ (see Mutation::Mixture::getComposition()).  Many of the
tools included with Mutation++ can also use this information to simplify input
on the command line.  An example of a list of element compositions for the
11-species Air mixture are shown below.

@code{.xml}
<element_compositions default="air1">
    <composition name="air1"> e-:0.0, N:0.79, O: 0.21 </composition>
    <composition name="air2"> e-:0.0, N:0.80, O: 0.20 </composition>
</element_compositions>
@endcode

Note that if the `default` attribute is used as above, the composition with the
name in the `default` value will be use as the default composition when computing
equilibrium calculations.

----
@section elements Elements

Element names and molecular weights are stored in the `elements.xml` file in the
`data/thermo` directory.  Each element is given as an XML node in the file.  An
example of the Nitrogen XML node is given as

@code{.xml}
<!-- Nitrogen -->
<element name="N">
    <mw units="g/mol">14.0067</mw>
</element>
@endcode

__Note that the electron is treated like an element in Mutation++ and should not be
changed__ in the `elements.xml` file.  In general, most of the known elements are
already given in `elements.xml` and it is unlikely that it should need to be
modified.

----
@section thermodynamic_databases Thermodynamic Databases

@subsection nasa_7 NASA 7-Coefficient Polynomials

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
Coefficients \f$a_1-a_5\f$ for upper temperature range           | `5(F15.8)`   | `1-75`
The integer '2'                                                  | `I1`         | `80`
<b>Line 3</b>                                                    | -            | -
Coefficients \f$b_1\f$ and \f$b_2\f$ for upper temperature range | `2(F15.8)`   | `1-30`
Coefficients \f$ a_1-a_3\f$ for lower temperature range          | `3(F15.8)`   | `31-75`
The integer '3'                                                  | `I1`         | `80`
<b>Line 4</b>                                                    | -            | -
Coefficients \f$a_4-a_5\f$ for lower temperature range           | `2(F15.8)`   | `1-30`
Coefficients \f$b_1\f$ and \f$b_2\f$ for lower temperature range | `2(F15.8)`   | `31-60`
The integer '4'                                                  | `I1`         | `80`


@code{.txt}
N                 L 6/88N   1    0    0    0G   200.000  6000.000 1000.        1
 0.24159429E+01 0.17489065E-03-0.11902369E-06 0.30226244E-10-0.20360983E-14    2
 0.56133775E+05 0.46496095E+01 0.25000000E+01 0.00000000E+00 0.00000000E+00    3
 0.00000000E+00 0.00000000E+00 0.56104638E+05 0.41939088E+01 0.56850013E+05    4
@endcode

@subsection nasa_9 NASA 9-Coefficient Polynomials

@cite McBride1996
Constants                                               | Format       | Columns
--------------------------------------------------------|--------------|--------
<b>Line 1</b>                                           | -            | -
Species name                                            | `A24`        | `1-24`
Comments (data source)                                  | `A56`        | `25-80`
<b>Line 2</b>                                           | -            | -
Number of \f$T\f$ intervals                             | `I2`         | `2`
Optional identification code                            | `A6`         | `4-9`
Chemical formulas, symbols, and numbers                 | `5(A2,F6.2)` | `11-50`
Zero for gas, nonzero for condensed phases              | `I1`         | `52`
Molecular weight                                        | `F13.5`      | `53-65`
Heat of formation at 298.15 K, J/mol                    | `F13.5`      | `66-80`
<b>Line 3</b>                                           | -            | -
Temperature range                                       | `2F10.3`     | `2-21`
Number of coefficients for \f$C_p^\circ/R_u\f$          | `I1`         | `23`
\f$T\f$ exponents in polynomial for \f$C_p^\circ/R_u\f$ | `8F5.1`      | `24-63`
\f$H^\circ (298.15)-H^\circ (0)\f$, J/mol               | `F15.3`      | `66-80`
<b>Line 4</b>                                           | -            | -
First five coefficients for \f$C_p^\circ/R_u\f$         | `5F16.8`     | `1-80`
<b>Line 5</b>                                           | -            | -
Last three coefficients for \f$C_p^\circ/R_u\f$         | `3F16.8`     | `1-48`
Integration constants \f$ b_1 \f$ and \f$ b_2 \f$       | `2F16.8`     | `49-80`
<i>Repeat 3, 4, and 5 for each interval...</i>          | -            | -


@code{.txt}
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
@endcode

@subsection rrho Rigid-Rotator Harmonic-Oscillators

@code{.xml}
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
@endcode


----
@section reaction_mechanisms Reaction Mechanisms

A __mechanism name__ _referenced in a mixture file will correspond to a file with an .xml
extension in the mechanisms directory_.  A mechanism file describes a reaction mechanism
which is a set of __elementary reactions__ of the form

\f[ \sum_i \nu'_{ij} \mathcal{A}_i \rightleftharpoons \sum_i \nu''_{ij} \mathcal{A}_i \f]

where \f$ \mathcal{A}_i \f$ represents the i'th species and \f$\nu'_{ij}\f$ and
\f$\nu''_{ij}\f$ are the forward and reverse stoichiometric coefficients of the i'th
species in the j'th reaction.

Mechanism files begin with a root element called `mechanism`.  A `name` attribute may be
given but it is ignored.  The name of the mechanism is determined by the filename (without
the `.xml` extension).

The root `mechanism` element may have any number of child elements which correspond to
__reaction elements__ or __unit specifiers__.

@subsection reaction_definitions Reaction Definitions

Below is an example `reaction` element
\code{.xml}
<reaction formula="N2+M=2N+M">
	<arrhenius A="7.0E+21" n="-1.6" T="113200." />
    <M>N2:3,N:4</M>
</reaction>
\endcode

`reaction` elements must contain a `formula` attribute which specifies the reaction
formula for the given reaction.  The `formula` string must obey the following rules:
- species names are separated the `+` symbol
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

@subsubsection rate_laws Rate Laws

A `reaction` element __must contain a child node which describes which rate law__ to use
when evaluating the forward rate coefficient.  The rate law node is specified by the tag
and its corresponding attributes correspond to the parameters for that given rate law.  The
following rate rate laws are currently supported.

- `arrhenius`  \f$ k_f(T_f) = A T_f^n \exp(\frac{E_a}{R_u T_f}) \f$
Att. | Value
-----|--------
`A`  | pre-exponential factor
`n`  | temperature exponent
`Ea` | activation energy
`T`  | characteristic temperature (\f$E_a/R_u\f$)
_note: only one of `Ea` or `T` may be used, not both_


@subsection unit_specifiers Unit Specifiers

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


@subsection example_mechanism Example Mechanism

As an example of how to create a mechanism file, consider the following example reaction
mechanism for a 5-species Nitrogen mixture of N2, N2+, N, N+ and e- with each reaction
controlled by an Arrhenius rate law.

\# | Formula                                      | A [mol,cm,s,K] | n     | Ea [K]
---|----------------------------------------------|----------------|-------|--------
1  | \f$N_2 + N_2 \rightleftharpoons 2N + N_2 \f$ | 1.0e21         | -1.6  | 113200
2  | \f$N_2 + N \rightleftharpoons 2N + N \f$     | 3.0e21         | -1.6  | 113200
3  | \f$N_2 + N^+ \rightleftharpoons 2N + N^+ \f$ | 1.0e21         | -1.6  | 113200
4  | \f$N_2 + e^- \rightleftharpoons 2N + e^-\f$  | 7.0e22         | -1.6  | 113200
5  | \f$N+e^-\rightleftharpoons N^++e^-+e^-\f$    | 2.5e30         | -3.82 | 168200
6  | \f$N+N\rightleftharpoons N_2^++e^-\f$        | 4.4e7          | 1.5   | 67500

Reactions 1-4 in the above example have the same rate constants except for the
pre-exponential factor differ only by the thirdbody species, making them good candidates
to combine into a generic thirdbody reaction with efficiency factors.  However, free
electrons cannot be thirdbodies in Mutation++, so reaction 4 must be treated separately.
Note also that there N2+ is not a possible thirdbody in reactions 1-4, thus its efficiency
factor should be zero.  The remaining reactions, 5-6, should be specified normally.  The
resulting mechanism file `example.xml` is given below.

\code{.xml}
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
\endcode

----
@section transfer_databases Energy Transfer Model Data

@section gsi Gas-Surface Interaction

The gas-surface interaction module of Mutation++ is responsinble for
treating chemically reacting surfaces in thermochemical non-equilibrium.
With the aim to obtain the appropriate surface boundary conditions
for a chemically reacting gas, the conservation equations around the surface
interface have to be written. In general, this idea can be applied for several
different categories of interfaces, with different phases and materials on each
side. Since this Mutation++ module emphases on gas-surface interaction it is limited
to having gas on one side (underscore g) and a catalytic or ablative solid or liquid
material on the other (underscore b). The surface is not simulated, but approximately modeled.
In most of the cases it can be considered
impermeable, or simplistically porous for the pyrolysis gases. The interface
can promote chemical reactions or radiate following the Planck's law of radiation
to the gas phase at the surface temperature. The set of balance equations obtained
are solved with respect to the conserved quantities of the gas, density,
momentum and energy, which is imposed as boundary values for the full Navier-Stokes equations,
or the porous material.

The general procedure to obtain the surface balances is the following.
Assuming steady state on the surface, the time derivative term, \f$\partial \mathcal / \partial t\f$,
is equal to zero. By taking the limit that one dimension of the volume goes to zero, the
three dimensional fluxes reduce to the fluxes normal to the interface while the volume
source terms such as chemical reactions rates go to zero; only the surface sources remain.
Generally, the set of balance equations can be written as:
\f[
[ \mathbf{F}_{\text{g}} - \mathbf{F}_{\text{g}} ] \cdot \mathbf{n} = \dot{\Omega}_{\text{s}},
\f]
where \f$\mathbf{F}\f$ are the fluxes from the gas (g) and the bulk (b) phases, while
\f$\dot{\Omega}_{\text{s}}\f$ are the source terms associated to the surface (s) processes.
Only the normal to surface flux component should be considered, denoted by
the inner product of the flux with the surface unit vector, (\f$\mathbf{n}\f$).
The specific form of the balance equations for mass and energy will be seen in the following
sections.

@subsection surf_chem Surface Chemical Production Terms

Below the input file when a single catalytic reactions is considered
modeled with the gamma (\f$ \gamma \f$) model with reaction probability equal to 1.

\code{.xml}
<gsi gsi_mechanism = "gamma">

    <surface_properties type = "none" >
    </surface_properties>

    <surface_chemistry>
        <!-- 1 -->
        <reaction type= "catalysis" formula="O + O => O2">
            <gamma_const> O:1. </gamma_const>
        </reaction>
    </surface_chemistry>

</gsi>
\endcode

The \f$ \gamma \f$ model, introduced by Goulard in the late 50s, is arguably the most popular
way to treat catalysis in the aero-thermodynamics community. It describes catalytic
reactions as macroscopic, non-elementary processes of the form:
\f[
A + A \rightarrow A_2
\f]
In order to determine the chemical production term for this kind of catalytic reactions a
probability for recombination \f$ \gamma \f$ is defined for each recombining species \f$ A \f$ as:

\f[
\gamma_{A} = \frac{ F_{A}, \text{rec} }{F_{A}^\downarrow},
\f]

where \f$ F_{A}^\downarrow \f$ is the flux of species \f$ A \f$ impinging the surface and
\f$ F_{{A}, \text{rec}}\f$ is the flux of species recombining at the surface. This probability is
the input parameter for the model. A fully catalytic wall has $ \gamma_{A} $ equal to 1, which
means that all the particles of species \f$ A \f$ impinging the surface recombine at the wall.
\f$ \gamma_{A} \f$ equal to 0 means that no reaction takes place and corresponds
to a non catalytic, or chemically inert wall. Anything between these two extreme cases is a
partially catalytic wall, which is the case for most of the surfaces. When the porbability
\f$ \gamma \f$ is defined, the surface chemical source term is determined as:
\f[
\dot{\omega}_{A,\text{cat}} = \gamma_{A} m_{A} F_{A}^\downarrow \quad \text{and} \quad
\dot{\omega}_{A_2,\text{cat}} = - \gamma_{A} m_{A} F_{A}^\downarrow
\f]
with \f$ m_{A} \f$ being the mass of species \f$ A \f$.

The only parameter that still needs to be defined is the impinging flux on the surface.
When the distribution function of species at the wall is well approximated by a Maxwellian
and there is no temperature slip, the impinging flux, \f$ F_{A}^\downarrow \f$, is equal to
\f[
F_{A}^\downarrow = n_{A} \sqrt{\frac{k_B T_{\text{s}}}{2 \pi m_{A}}}
\f]
according to the kinetic theory of gases.

Initially, the gamma model was developed to describe homonuclear reactions, such as the one presented in above.
Soon, though, it was observed that heteronuclear reactions were also probable
to occur at the wall, which can take the form of
\f[
A + B \rightarrow AB,
\f]
with the chemical rate for the produced molecule \f$ AB \f$ equal to:
\f[
\dot{\omega}_{AB,\text{cat}}  = - \gamma_{AB} m_{A} F_{A}^\downarrow - \gamma_{BA} m_{B} F_{B}^\downarrow.
\f]
One atom of \f$ A \f$ recombines with one atom of \f$ B \f$ on the surface to produce a molecule \f$ AB \f$.
In other words the number of atoms of \f$ A \f$ that recombine into \f$ AB \f$ should be equal to the number of \f$ B \f$
atoms that recombine into \f$ AB \f$. This restriction should be explicitly imposed in order to conserve mass:
\f[
\gamma_{AB} F_{A}^\downarrow = \gamma_{BA} F_{B}^\downarrow.
\f]
As a result, in a reaction of this type two gamma recombination coefficients should be defined not
necessarily equal for the two processes, the one activated when the catalytic reaction is limited by the
flux of \f$ A \f$ atoms and the opposite. In practice the two recombination number fluxes,
\f$ \gamma_{AB} F_{A}^\downarrow \f$ and \f$ \gamma_{BA} F_{B}^\downarrow \f$, are compared and the limiting one
determines which of the two gammas is chosen.

These gamma coefficients cannot take arbitrary values, they should be limited between 0 and 1,
just like in the homonuclear case. When extra catalytic reactions are added, such as
\f[
A + A \rightarrow A_2 \quad \text{and} \quad  B + B \rightarrow B_2
\f]
the gammas should be further constrained, as:
\f[
0 \le \gamma_{A} + \gamma_{AB} \le 1 \quad \text{and} \quad 0 \le \gamma_{B} + \gamma_{BA} \le 1,
\f]
in order for mass to be conserved. This approach is consistent with similar approaches
considering the catalytic recombination occurring in Martian atmospheres, where \f$ O \f$ can recombine
into both  \f$ O_2 \f$ and  \f$ CO_2 \f$ due to catalytic reactions.

An example input file of catalytic reactions in air including the formation of \f$ NO \f$ are below.
\code{.xml}
<gsi gsi_mechanism = "gamma">

    <surface_properties type = "none" >
    </surface_properties>

    <surface_chemistry>
        <!-- 1 -->
        <reaction type= "catalysis" formula="N + N => N2">
            <gamma_const> N:.001 </gamma_const>
        </reaction>
        <!-- 2 -->
        <reaction type= "catalysis" formula="O + O => O2">
            <gamma_const> O:.001 </gamma_const>
        </reaction>
        <!-- 3 -->
        <reaction type= "catalysis" formula="N + O => NO">
            <gamma_const> N:.002 O:.003 </gamma_const>
        </reaction>
    </surface_chemistry>

</gsi>
\endcode

It is still unclear which is the proper boundary conditions for the
ions and electrons on the surface. One of the most common model is to
assume full ion recombination on the surface can be expressed as
catalytic reaction with probability 1.
The following example shows how to impose full ion recombination in
Mutation.
\code{.xml}
<gsi gsi_mechanism = "gamma">

    <surface_properties type = "none" >
    </surface_properties>

    <surface_chemistry>
        <!-- 1 -->
        <reaction type= "catalysis" formula="N + N => N2">
            <gamma_const> N:.001 </gamma_const>
        </reaction>
        <!-- 2 -->
        <reaction type= "catalysis" formula="O + O => O2">
            <gamma_const> O:.001 </gamma_const>
        </reaction>
        <!-- 3 -->
        <reaction type= "catalysis" formula="N + O => NO">
            <gamma_const> N:.002 O:.003 </gamma_const>
        </reaction>
        <!-- 4 -->
        <reaction type= "catalysis" formula="N+ + e- => N">
            <gamma_const> N+:1. e-:1. </gamma_const>
        </reaction>
        <!-- 5 -->
        <reaction type= "catalysis" formula="O+ + e- => O">
            <gamma_const> O+:1. e-:1. </gamma_const>
        </reaction>
        <!-- 6 -->
        <reaction type= "catalysis" formula="NO+ + e- => NO">
            <gamma_const> NO+:1. e-:1. </gamma_const>
        </reaction>
        <!-- 7 -->
        <reaction type= "catalysis" formula="N2+ + e- => N2">
            <gamma_const> N2+:1. e-:1. </gamma_const>
        </reaction>
        <!-- 8 -->
        <reaction type= "catalysis" formula="O2+ + e- => O2">
            <gamma_const> O2+:1. e-:1. </gamma_const>
        </reaction>
    </surface_chemistry>

</gsi>
\endcode

An example of ablation model can be seen below.
\code{.xml}
<gsi gsi_mechanism = "gamma">

    <surface_properties type = "ablation" >
        <surface label="b" species="C" />
    </surface_properties>

    <surface_chemistry>
        <!-- 1 -->
        <reaction type= "ablation" formula="C-b + O => CO">
            <gamma_T pre_exp="0.63" T="1160.0" />
        </reaction>
        <!-- 2 -->
        <reaction type= "ablation" formula="C-b + N => CN">
            <gamma_const> N:0.003 </gamma_const>
        </reaction>
        <!-- 3 -->
        <reaction type= "ablation" formula="3C-b => C3">
            <sublimation vap_coef="0.1" pre_exp="5.19E15" T="90845." />
        </reaction>
        <!-- 4 -->
        <reaction type= "catalysis" formula="N+N=>N2">
            <gamma_const> N:0.001 </gamma_const>
        </reaction>
    </surface_chemistry>

</gsi>
\endcode

@subsubsection coxid Carbon Oxidation Model

The first ablation reaction presented above is the oxidation
of the solid carbon by atomic oxygen. The reaction reads as:

\f[
C_{\text{b}} + O \rightarrow CO
\f]

and is exothermic, releasing \f$ 3.74 \f$ eV per molecule produced. Its reaction
rate coefficient is given by defining a recombination probability
\f$ \gamma_{CO} \f$. This probability is an Arrhenius type function of temperature
and is given by the formula:
\f[
\gamma_{CO} = 0.63 \exp(-1160/T_{\text{s}}).
\f]
Carbon oxidation with molecular oxygen is also possible
(\f$ C_{\text{b}} + O_2 \rightarrow CO + O \f$), but it is considered a less significant
process.

@subsubsection cnitr Carbon Nitridation Model

The second shown in the example is carbon nitridation,
\f[
C_{\text{b}} + N \rightarrow CN,
\f]
an exothermic reaction with \f$ 0.35 \f$ eV of energy released per reacting atom.
The reaction rate is given using a constant recombination probability
\f$ \gamma_{CN} \f$ like in the catalytic case.

@subsubsection csubl Carbon Sublimation Model

At high temperatures carbon removal from the surface is dominated by phase
change processes like sublimation. The production of \f$ C_3 \f$ is considered here.
It should be noted that this type of reactions are invertible, with formula:
\f[
3C_\text{b} \rightleftharpoons C_3.
\f]
The chemical production rate of this reaction is equal to:
\f[
\dot{\omega}_{\text{subl},C_3} = \beta_{C_3} ( \rho_{C_3, \text{equil}} - \rho_{C_3})
 \sqrt{ \frac{ k_B T_\text{s}}{2 \pi m_{C_3}}}.
\f]
The equilibrium partial density of the \f$ C_3 \f$ species is obtained from the
saturated vapor pressure of carbon, which is equal to
\f[
p_{C_3, \text{sat}} = c \exp (-T_{\text{act}}/T_{\text{w}})
\f]
with \f$ \beta_{C_3} \f$ being the evaporation coefficient, \f$ c \f$ the pre-exponential coefficient,
and \f$  T_{\text{act}}  \f$ the activation temperature.
Even though here only sublimation is presented in the example, evaporation processes
can be considered with the same model.

@subsection smb Surface Mass Balance

@subsubsection smbc Surface Mass Balance for Catalytic Surfaces

Heterogeneous catalysis is an important gas-surface interaction phenomenon occurring
during re-entry of vehicles equipped with re-usable thermal protection system. It
describes the recombination of the dissociated atoms in the flow using the thermal
protection system as a catalyst. It is called heterogeneous, because the recombining
species and the catalyst are in a different phase, here gas and solid. The catalyst,
without being consumed, increases the rate of chemical reactions by offering an alternative,
energetically favored path. It is important to note that catalysis does not change the chemical
equilibrium of reactions, since it favors equally both the forward and the backward reaction rates.
This is a constrain that should be respected by the catalytic model chosen.

In hypersonic the recombination reactions that occur on the surface are in general exothermic.
The energy released to the wall is a substantial percentage of the total heat flux experienced
by space vehicles. Not necessarily all of the recombination energy is directly deposited to the
surface. A part of it is used to excite the internal energy of the produced molecules. Another
reason why calculating
the actual heat released on the surface is a complicated task, is that the gas phase chemistry
and diffusion play an important role in determining the overall catalytic rates.
If all of the phenomena above are modeled with accuracy, the re-entry heat load can be predicted
and the size of the thermal protection system can be determined.
The mass balance on a catalytic surface reads:
\f[
\mathbf{j}_i \cdot \mathbf{n} = \dot{\omega}_{\text{cat},i},
\f]
where one mass balance equation should be solved for each distinct species in the flow.
The equation above states, that the catalytic activity of every species is equal to the diffusion
flux of these species to the surface. This leads as to two cases. In the first one, the rate
with which the chemical species are produced or destroyed at the wall is higher than the rate
they diffuse to the surface while in the second one the opposite happens. The first case is called
diffusion limited, since diffusion is the mechanism controlling the chemical process, while the
opposite is called reaction limited. When a species $i$ is inert at the surface, then its chemical
rate is exactly equal to zero, which at steady state imposes that its net diffusion flux is also zero.
Even though, in principle, these equations could be omitted, since they impose that the mole fractions
of the species in question do not change with respect to the ones in the gas phase, the full system was
chosen to be solved.

@subsubsection smbac Surface Mass Balance for Ablative and Catalytic Surfaces

Ablation is the chemical gas-surface interaction phenomenon which occurs on non-reusable thermal
protection systems of re-entry vehicle. The term ablation describes the category of chemical
reactions during which the dissociated atoms in the flow field recombine directly with the
thermal protection system, which by burning protects the vehicle. Contrary to catalysis, this
burning destroys the material itself making it unable for reuse. Not only chemical reactions can
cause the degradation of ablative thermal protection systems. Mechanical removal processes, such
as spillation, can also occur. The particle injected in the flow field due to these phenomena are
not necessarily in a gaseous form and their modeling requires approaches beyond the scope of this
work. The type of chemical reactions studied here are assumed to produce only gaseous species
and occur only the surface of the material. In cases where the material is porous, ablation
processes also in the bulk, such as in the case of pyrolysis.

Taking these ideas in mind, the surface mass balance accounting only for surface reactions becomes:
\f[
[ \rho_i (\mathbf{u}_{\text{g}} - \mathbf{u}_{\text{r}} ) + \mathbf{j}_i - \mathbf{F}_{\text{b},i} ]
\cdot \mathbf{n} = \dot{\omega}_{\text{s},i},
\f]
with
\f$ \dot{\omega}_{\text{s},i} = \dot{\omega}_{\text{cat},i}  + \dot{\omega}_{\text{abl},i} \f$ and
term \f$ \mathbf{F}_{\text{b},i} \f$ is the flux of species \f$ i \f$ entering the interface due
to solid process, like pyrolysis and solid-solid chemical reactions.
Just like before, one mass balance equation should be solved for each distinct species in the flow.
It is most of the time reasonable to consider that recession velocity is
orders of magnitude lower than the gas velocity \f$ \mathbf{u}_{\text{r}}  \f$.

@subsection seb Surface Energy Balance Solver

In order to determine the surface temperature a surface energy balance should be solved along
with the mass balances. It takes the form:
\f[
[ \rho ( \mathbf{u}_{\text{g}} - \mathbf{u}_{\text{r}} ) H + \mathbf{q}_{\text{g}}
    - \mathbf{F}_{\text{b},{e}} ] \cdot \mathbf{n} =
    \dot{\Omega}_{\text{s},e}
\f]
where \f$ \mathbf{q}_{\text{g}} \f$ is the heat flux to the gas phase equal to:
\f[
    \mathbf{q}_{\text{g}} = -\lambda \nabla T + \sum_{i = 1}^{n_\text{s}} \mathbf{j}_i h_i.
\f]
The radiative heat flux can be also included and will be discussed along with the surface
radiation.

Term \f$ \mathbf{F}_{\text{s},{e}} \f$ describes the energy exchanged between the interface
and the bulk of the solid and is composed of three contributions: the first one is the
thermal conduction exiting the surface \f$ \mathbf{q}_{\text{b,cond}} \f$, the second one is the
enthalpy entering the interface due to the movement of the surface with the recession velocity,
\f$ \mathbf{u}_{\text{r}} \rho_{\text{s}} h_{\text{s}} \f$ and the third one appears only in cases of
porous material, describing the enthalpy of the solid pyrolysis gases convected in the interface,
denoted as \f$ \mathbf{u}_{\text{p}} \rho_{\text{p}} h_{\text{p}} \f$. The subscript p symbolizes
the pyrolysis gas properties, with the \f$ \rho_{\text{p}} h_{\text{p}} \f$ being actually the
sum \f$ \sum_{i=i}^{n_\text{s}} = \rho_i h_i \f$ for the pyrolysis gas densities. The surface
enthalpy \f$ h_s \f$ is an input to the code, as an attribute to the surface_properties element.

An example input file for solving both mass and energy balance can be seen below.
\code{.xml}
<gsi gsi_mechanism = "gamma_energy">

    <surface_properties type = "ablation" >
        <surface label="b" species="C"
                 enthalpy_surface="0." />
    </surface_properties>

    <solid_properties
        virgin_to_surf_density_ratio = "1."
        enthalpy_virgin = "0."
    />

    <surface_features
        solid_conduction = "steady_state"
        surface_in_thermal_equil = "true"
        gas_radiation = "false"
    />

    <surface_chemistry>
        <!-- 1 -->
        <reaction type= "ablation" formula="C-b + O => CO">
            <gamma_T pre_exp="0.63" T="1160.0" />
        </reaction>
        <!-- 2 -->
        <reaction type= "ablation" formula="C-b + N => CN">
            <gamma_const> N:0.003 </gamma_const>
        </reaction>
        <!-- 3 -->
        <reaction type= "ablation" formula="3C-b => C3">
            <sublimation vap_coef="0.1" pre_exp="5.19E15" T="90845." />
        </reaction>
        <!-- 4 -->
        <reaction type= "catalysis" formula="N + N => N2">
            <gamma_const> N:0.001 </gamma_const>
        </reaction>
    </surface_chemistry>

    <surface_radiation emissivity="0.86" T_env="0."/>
</gsi>
\endcode

The only source term taken into account for the surface energy balance in the
example above is radiation. It is considered by adding a new element with tag
<surface_radiation/>. The surface is assumed to be in thermodynamic equilibrium
at a temperature \f$ T_{\text{s}} \f$  emitting energy following the Stefan-Boltzmann law,
\f[
\dot{\Omega}_{\text{s},e} = \dot{\Omega}_{\text{rad}} = \epsilon \left( \sigma T_{\text{s}}^4
- q_{\text{rad},g} \right)
\f]
with $\sigma$ being the Stefan-Boltzmann constant and \f$ \epsilon \f$ being the emissivity
coefficient.  and \f$ T_{\text{env}} \f$ is the surrounding environment temperature.

\f$ \sigma T_{\text{env}}^4 \f$.

In order to obtain the value for the conductive heat flux on the solid,
either a material code should be used or one should at least solve the
energy equation in the solid. Even though these approaches can be very accurate,
in cases where the material has low thermal conductivity or the recession rates
are high, approximate methods can be used, without compromising the accuracy of
the simulations. Such an approximation was adopted here for the modeling of the
conductive heat flux inside the material, the steady-state ablation approach for a
semi-infinite surface. By writing the steady-state energy equation for the solid phase
and integrating over the semi-infinite material, where on one side is the solid
properties at exactly the interface (subscript s), while on the other, at infinity,
is the virgin material (subscript v), the steady state heat flux is given by the formula:
\f[
q_{\text{cond}}^{\text{SS}} = -\mathbf{u}_{\text{r}}(\rho_\text{v} h_{\text{v}} - \rho_\text{s} h_{\text{s}}).
\f]
By replacing the formula above in the surface energy balance and after simplifications one gets:
\f[
[
\rho ( \mathbf{u}_{\text{g}} - \mathbf{u}_{\text{r}} ) H + \mathbf{q}_{\text{g}}
+ \mathbf{u}_{\text{r}}\rho_\text{v} h_{\text{v}}
] \cdot \mathbf{n} =
\dot{\Omega}_{\text{s},e},
\f]
which is equally valid for both porous and non-porous materials. In order to use the steady state ablation
approximation, the attribute of surface_feature element solid_conduction should be set to steady state.
Instead of inputing the
virgin material density, the ratio between the virgin and surface density minus one is often used,
refered to as
\f[ \phi = ( \frac{\rho_\text{v}}{\rho_\text{s}} - 1) \f] and is the attribute
virgin_to_surf_density_ratio in the solid_properties element defaulting to 1.
The enthalpy_virgin can also be an input with default
value equal to 0. Note that these two last options are only necessary when the steady
state assumption for the solid conduction is considered.

When a material solver is available the conductive heat flux should be an input to the
library for increased accuracy. This can be achieved by setting the solid_conduction
surface feature to "input" and the setSolidCondHeatFlux function can be invoked.
The option enthalpy_surface should be set in this case, otherwise it is automatically set
to zero.

Note that only one surface balance equation is solved regardless of the thermodynamic
state model. When multitemperature models are considered for the thermodynamics additional
temperatures should be imposed on the surface. It is often a reasonable assumption to
impose thermal equilibrium at the wall. This is achieved by setting the surface_feature
option surface_in_thermal_equil "true". If it is "false", then any additional temperature
beyond translations will be left unchanged, an assumption which can be used to impose an
adiabatic boundary condition for the internal energy modes.

@subsection gsi_example Using the GSI module

In order for Mutation++ to take into account the Gas Surface Interaction
features the gsi_mechanism attribute should be assigned to the gsi
mechanism file.
@code{.xml}
<!-- Example mixture file for gsi -->
<mixture gsi_mechanims="name_of_the_gsi_file">
    <!-- Species list -->
    <species>
        ...
    </species>
</mixture>
@endcode

Below an example code of how to use the gas surface interaction features
of Mutation is presented.
@code{.cpp}
    const int set_state_with_rhoi_T = 1;

    MixtureOptions opts("mixture_name");
    Mixture mix(opts);

    const int iter = 5;
    mix.setIterationsSurfaceBalance(iter);

    // Setting the state and setting up the library
    mix.setState(rhoi_surf.data(), Tsurf.data(), set_state_with_rhoi_T);
    mix.setSurfaceState(rhoiw.data(), Tsurf.data(), set_state_with_rhoi_T);
    mix.setDiffusionModel(xi_edge.data(), dx);
    mix.setGasFourierHeatFluxModel(Tedge.data(), dx); // Only works with energy balance
    // The .data() function returns the pointer at the first element of the data container.

    // Additional Options
    double gas_rad = 0.;
    setGasRadHeatFlux(*gas_rad); // Only called if feature gas_radiation is true
    double solid_cond = 0.;
    setSolidCondHeatFlux(*solid_cond); // Only called if feature solid_cond is set to input

    // Solving the surface Mass Balance and requesting the solution
    mix.solveSurfaceBalance();
    mix.getSurfaceState(rhoiw.data(), Tsurf.data(), set_state_with_rhoi_T);

    // Getting mass blowing rate.
    double mblow;
    mix.getMassBlowingRate(mblow);

    // Getting surface reaction rates.
    mix.getSurfaceReactionRates(wdot.data());

    // Getting number of reactions and surface reaction rates per reaction.
    int m_nr = mix.nSurfaceReactions();
    mix.getSurfaceReactionRatesPerReaction(wdot_reac.data());
@endcode

Note that this example is not supposed to compile or run, but is there to indicate
the most important features of the library.

*/
