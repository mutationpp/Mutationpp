/**
@page input_files Input Files

@tableofcontents

@section input_intro Introduction

Data files should be provided in one of the subdirectories under the `data` 
directory as follows:

- __data__
  + __mechanisms__ chemical [reaction mechanisms](@ref reaction_mechanisms)
  + __mixtures__ [mixture definitions](@ref mixtures)
  + __thermo__ [elemental](@ref elements) and species [thermodynamic databases](@ref thermodynamic_databases)
  + __transfer__ internal energy [transfer model databases](@ref transfer_databases)
  + __transport__ [collision integral database](@ref collision_integrals)

The path to the __data__ directory is given by the `MPP_DATA_DIRECTORY` environment
variable which should have been setup during the [installation](@ref installation) of Mutation++.
Note that [mixture files](@ref mixtures) may also be in the local directory where
an excutable that utilizes Mutation++ is being run.

@subsection simple_xml Simplified XML Language
Many of these files are written in a simplified version of the Extensible Markup 
Language (XML).  XML provides a human readable, yet complex and extensible format for data
to be stored with only a few, limited rules.

\code{.xml}
<!-- Comment string -->
<root_tag attribute="value">
    <child1_tag>
        text
    </child1_tag>
    <child2_tag attribute="another value" />
</root_tag>
\endcode

The above code fragment shows a small example of how the simplified XML format works for 
Mutation++.  An XML document begins with a root XML element.  Every element must begin 
with a tag that identifies what type of element it is.  The root element depicted in 
the example is of type `root_tag`.  Every element also ends with an end-tag 
which signifies the end of the element.  Each element may have as many attribute/value 
pairs as is desired following the element's tag.  Each pair must consist of an attribute 
name followed by an equal sign and the value of the attribute in quotations.  Between the 
begin- and end-tags of an element, an element may also contain one or more child elements 
or text (but not both).  From the figure, the root element contains two child elements 
named `child_tag` and `child2_tag`.  Note that the first child element is an example of an 
element which contains text instead of more child elements.  The second child element is 
an example of a short-hand format for elements which only contain attributes.  For such 
elements, a full end-tag is not necessary.  Instead, simply putting `/>` after the attribute 
list is sufficient to end the element.  Finally, comments can be inserted anywhere outside
of element tags.  Comment strings begin with `<!-``-` and end with `-``->` and can be spread 
over multiple lines.

@section mixtures Mixtures

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

@section thermodynamic_databases Thermodynamic Databases

@subsection nasa_7 NASA 7-Coefficient Polynomials

@subsection nasa_9 NASA 9-Coefficient Polynomials

@subsection rrho Rigid-Rotator Harmonic-Oscillators


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
`T`  |characteristic temperature (\f$E_a/R_u\f$)
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
6  | \f$N+N\rightleftharpoons N2^++e^-\f$         | 4.4e7          | 1.5   | 67500

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


@section collision_integrals Collision Integrals

Collision integral data is currently the only input data which does not use the
[XML format](@ref simple_xml).  This will be updated in the future to allow a more
expressive data format for a wide range of integrals.

For now however, all collision integral data should be given in the
`data/transport/heavy.dat` file.  For each binary interaction, the \f$\Omega_{ij}^{11}\f$,
\f$\Omega_{ij}^{22}\f$, and \f$B^*_{ij}\f$ functions are expressed in terms of 
an exponential polynomial curve-fit versus temperature with units of Angstroms squared.

\f[\Omega_{ij} = \exp[c_1(\ln T)^3 + c_2(\ln T)^2 + c_3(\ln T) + c_4] \f]

An example entry in the `heavy.dat` file is given below for the binary collision
of two Nitrogen molecules.

     .N2-N2.
       -0.0065954     0.13921     -1.1559      6.9352
       -0.0087373     0.19482     -1.6023      8.1845
                0           0           0      0.1398


@section transfer_databases Energy Transfer Model Data

*/
