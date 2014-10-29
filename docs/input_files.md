Input Files
===========

Introduction
------------
This section details how to create element, species, mixture, and reaction mechanism data 
that will be loadable by Mutation++.  In general, data files are listed in 
one of the subdirectories under the _data_ directory.  

- data
  + __mechanisms__ _chemical reaction mechanisms_
  + __mixtures__ _mixture definitions_
  + __thermo__ _elemental and species thermodynamic databases_
  + __transfer__ _internal energy transfer model data_ 
  + __transport__ _collision integral data_

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

Reaction Mechanisms
-------------------

\code{.xml}
<!-- Park Air-11 Reaction Mechanism from Park 2001:               -->
<!--     Journal of Thermophysics and Heat Transfer Vol. 15 No. 1 -->
<mechanism name="air11_park01">
    
    <arrhenius_units A="mol,cm,s,K" E="kcal,mol,K" />
    
    <!-- 1a Park 2001-->
    <reaction formula="N2+M=2N+M">
        <arrhenius A="7.0E+21" n="-1.6" T="113200." />
        <M>N:4.28571428571, O:4.28571428571, O2+:0.0</M>
    </reaction>
    
    <!-- 1b Park 2001-->
    <reaction formula="N2+e-=2N+e-">
        <arrhenius A="3.0E+24" n="-1.6" T="113200." />
    </reaction>

    <!-- 2 Park 2001-->
    <reaction formula="O2+M=2O+M">
        <arrhenius A="2.0E+21" n="-1.5" T="59360."/>
        <M>N:5.0, O:5.0, O2+:0.0</M>
    </reaction>

    <!-- 3 Park 2001-->
    <reaction formula="O+e-=O++e-+e-">
        <arrhenius A="3.9E+33" n="-3.78" T="158500." />
    </reaction>

    <!-- 4 Park 2001-->
    <reaction formula="N+e-=N++e-+e-">
        <arrhenius A="2.5E+34" n="-3.82" T="168200." />
    </reaction>

    <!-- 5 Park 2001-->
    <reaction formula="N+O=NO++e-">
        <arrhenius A="5.3E+12" n="+0.00" T="31900." />
    </reaction>

    <!-- 6 Park 2001-->
    <reaction formula="N+N=N2++e-">
        <arrhenius A="4.4E+07" n="+1.50" T="67500." />
    </reaction>

    <!-- 7 Park 2001-->
    <reaction formula="N2+O=NO+N">
        <arrhenius A="5.70E+12" n="+0.42" T="42938." />
    </reaction>

    <!-- 8 Park 2001-->
    <reaction formula="O+NO=O2+N">
        <arrhenius A="8.40E+12" n="+0.00" T="19400." />
    </reaction>

</mechanism>
\endcode

