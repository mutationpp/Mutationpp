<a id="top"></a>

# checkmix
This is a utility program which loads a mixture and prints information about the various elements, species, and reactions in the [mixture](input-files.md#mixtures).  A mixture is loaded just as it would be in any other application, so this tool is useful to check for any syntax errors or missing data in a given mixture before using it elsewhere.  It can also be used to see [exactly which order species](input-files.md#species-order) and reactions are stored internally in Mutation++ for a given mixture.

## Usage
**Method 1:** From a Mixture file
```bash
checkmix mixture
```

**Method 2:** From a Species Descriptor
```bash
checkmix [NASA-7 | NASA-9 | RRHO] "species-descriptor"
```
`species-descriptor` should follow the rules given [here](@ref species_list).

## Example
Assuming you have a [mixture file](input-files.md#mixtures) called `air5.xml` and a
[mechanism file](input-files.md#reaction-mechanisms) called `air5-chem.xml` (inther
your local directory or the `data/mixtures` and `data/mechanisms`
directories) which are as follows
```xml
<mixture name="air5" mechanism="air5-chem">
    <species>
        N O NO N2 O2
    </species>
</mixture>
```
and
```xml
<mechanism name="air5">
    <arrhenius_units A="mol,cm,s,K" E="kcal,mol,K" />
    <!-- 1 -->
    <reaction formula="N2+M=2N+M">
        <arrhenius A="3.0E+22" n="-1.6" T="113200.0" />
       <M>N2:0.2333, NO:0.2333, O2:0.2333</M>
    </reaction>
    <!-- 2 -->
    <reaction formula="O2+M=2O+M">
        <arrhenius A="1.0E+22" n="-1.5" T="59360.0" />
        <M>N2:0.5, NO:0.5, O2:0.5</M>
    </reaction>
    <!-- 3 -->
    <reaction formula="NO+M=N+O+M">
        <arrhenius A="5.0E15" n="+0.0" T="75500.0" />
        <M>NO:20.0, N:20.0, O:20.0</M>
    </reaction>
    <!-- 4 -->
    <reaction formula="N2+O=NO+N">
        <arrhenius A="5.69E+12" n="+0.42" T="42938.0" />
    </reaction>
    <!-- 5 -->
    <reaction formula="O2+N=NO+O">
        <arrhenius A="2.49E+09" n="+1.18" T="4005.5" />
    </reaction>
</mechanism>
```
then the command `checkmix air5` will produce the following output:

```
5 species containing 2 elements
5 reactions
Species info:
-------------
    N   O  Mw (g/mol)    Charge       Phase
Gas Species (5):
N     1   0     14.0067         0         gas
O     0   1     15.9994         0         gas
NO    1   1     30.0061         0         gas
N2    2   0     28.0134         0         gas
O2    0   2     31.9988         0         gas
Default elemental composition:
------------------------------
   N  :   0.5
   O  :   0.5
Reaction info:
--------------
Type ID Key
   3: heavy particle impact dissociation
   7: exchange
Reactions
   #  Formula             Type  Rate Law     A (m,s,mol)      n    Ta (K)
   1: N2+M=2N+M           3     Arrhenius:     3.000e+16  -1.60  113200.0
      N2: 0.23, NO: 0.23, O2: 0.23
   2: O2+M=2O+M           3     Arrhenius:     1.000e+16  -1.50   59360.0
      N2: 0.50, NO: 0.50, O2: 0.50
   3: NO+M=N+O+M          3     Arrhenius:     5.000e+09   0.00   75500.0
      NO: 20.00, N: 20.00, O: 20.00
   4: N2+O=NO+N           7     Arrhenius:     5.690e+06   0.42   42938.0
   5: O2+N=NO+O           7     Arrhenius:     2.490e+03   1.18    4005.5
```
