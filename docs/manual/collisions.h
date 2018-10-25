/**
@page collisions Collision Integral Database

@tableofcontents

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
*/