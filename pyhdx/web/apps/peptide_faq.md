FAQ
===

### What does this app do?

This app calculates the D-uptake of a peptide, given the parameters of the Linderstrøm-lang
model, $$k_{open}$$ and $$k_{close}$$, as well as the experimental parameters pH and temperature.
This allows HDX-MS users to find which $$\Delta G_{EX}$$ (Hereafter $$\Delta G$$) values can be probed, 
given time range, pH and temperature.

For example, with pH 6.5 and a temperature of 5&deg;C, if $$\Delta G$$ for an amino 
acid is $$\gt 35$$ kJ/mol, no exchange is expected before $$10^5$$ seconds (~ 1 day) of exchange.

### What is this 'dependent variable' in peptide controls?

The parameters for the model are $$k_{open}$$ and $$k_{close}$$, which relate to $$\Delta G$$
via $$\tfrac{k_{close}}{k_{open}} = e^{\frac{\Delta G}{RT}} $$. Therefore, if two out of three
quantities are known, the third can be calculated. The 'dependent variable' settings lets you
choose which one is calculated, while the other two can be set by the user.

### How do you know what values for $$k_{open}$$ are?

We don't, typically only the ratio of $$k_{open}$$ and $$k_{close}$$ can be determined from 
HDX experiments. Literature information on opening and closing rates is scarce, so for I'm 
only aware of a paper from 1993 ([Pederson](https://doi.org/10.1006/jmbi.1993.1176)) and a
more recent work from 2021 ([Peacock](https://10.1021/acs.biochem.1c00277)). 

The former finds $$k_{open}}$$ rates ranging from ~$$10^{-3}$$ to $$10^{-8} s^{-1}$$. However,
with such a slow opening rate reaction, to satisfy typical $$\Delta G$$ values of 10-50 kJ/mol,
the corresponding $$k_{close}$$ rates would be 0.1-100 $$s^{-1}$$, and given that typical
values of $$k_{int}$$ (Intrinsic/chemical rate of exchange) are similar or even lower, this
would place most of these HD kinetics outside the EX2 regime (where $$k_{close} \gg k_{int}$$), 
contrary to observation.

The latter work by Peacock at al finds $$k_{open} \approx 10^2 s^{-1}$$, and for reasonably
flexible proteins/residues with $$\Delta G$$ = 15 kJ/mol, this would mean $$k_{close}$$ of
~$$10^{5} s^{-1}$$. The protein reported is likely more flexible than $$\Delta G$$ = 15 kJ/mol, 
and assuming flexibility is reduced by a decrease in $$k_{open}$$, rates for $$k_{close}$$ 
are will be slower than $$10^{5} s^{-1}$$.


### What are the timescales involved?

If we do an experiment at 20 &deg;C, pH 7.5 with a middle-of-the road protein flexibility of
$$\Delta G$$ = 20 kJ/mol D-uptake is expected to start at about 100s of D-exposure, where 
~95% D is reached at roughly 1000s. If we assume $$k_{open}$$ is $$10^{2} s^{-1}$$, $$k_{close}$$ 
is about $$10^6 s^{-1}$$. This means that the lifetime of the open state is on the order of 
microseconds, and the opening event will occur about 10000 times before we start to see a
measurable amount of D-uptake, after 100 seconds of D-exposure. 


### What about EX1 vs EX2?

In this simulation we know both $$k_{open}$$ and $$k_{close}}$$, so there is no need to make
mathematical approximations depending on the kinetics regime. D-uptake is calculated from:

$$
k_{ex} = \frac{1}{2} (k_{open} + k_{close} + k_{int} - \sqrt{((k_{open} + k_{close} + k_{int})^2 - 4 k_{open} k_{int}})
$$

As described in the Pederson paper, and is accurate for at least the $$ \Delta G \gt 10 $$ kJ/mol regime.

Since in EX2 the exchange only depends on the ratio of the  closing and opening rates, and 
not their absolute values, we can probe in which regime we are by playing with the sliders. 
Typically, with $$\Delta G$$ at 30 kJ/mol, and 'Dependent variable' set to $$k_{close}$$, 
we see no effect in the uptake curve (lower panel for per-residue uptake curve) when changing
the rate of $$k_{open}$$ over 6 order of magnitude from $$10^{-2}$$ to $$10^4 s^{-1}$$. 
However, at about 15 kJ/mol, changes in the D-uptake curve begin to show as function of 
$$k_{open}$$, with about a 2-fold difference in the half-life of D-uptake. These differences
become more pronounced at 10 kJ/mol (up to 10-fold) and relate to accuracy of $$\Delta G$$ values
found by PyHDX in these low $$\Delta G$$ / high flexibility regimes given that PyHDX assumes
EX2 kinetics. 

### The first amino acid is missing!

For simplicity, we are only considering 'middle' amino acids, so when calculating intrinsic 
rates of exchange no N-terminal or C-terminal effects are taken into account. However, since 
these rates do also depend on the chemical nature of the amino acid on the N-terminal side of
the amide in question, we need to know which amino acid this is. So to model exchange for
the peptide 'LGPLTAGHH', find in your protein which amino acid is N-terminal to this 
peptide, and then for example enter 'KLGPLTAGHH'.

### What about back exchange?

There is no back exchange modelled here. 

### Why do the sliders look funny/colored or move janky?

There are some bugs with the sliders and possibly the way they are coupled to eachother is
not ideal, we're looking at improving this at some point.

### Can I import my data and compare to simulation?

No sorry, at the moment this is not possible. We always welcome feature requests through the 
[Issues](https://github.com/Jhsmit/PyHDX/issues) page on GitHub if you any suggestions to improve PyHDX. 

### Can I change the limits of the time axis or slider values?

Same as the previous question, although this one is a little higher on the priority list. 
