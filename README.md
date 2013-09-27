# xeno-fix
This is a clone of the [`XenoCat`](https://code.google.com/p/r-xenocat/) R package, updated to work with the latest versions of the [`lme4`](https://github.com/lme4/lme4) package. `XenoCat` is a great R package for analyzing tumor growth experiments -- stats described [here](http://clincancerres.aacrjournals.org/content/early/2012/07/19/1078-0432.CCR-11-3215.full.pdf). It is much better than running a t-test at some bogus timepoint you chose before your xenograft experiment started. 

## What's going on?
I really like the `XenoCat` package and use it regularly in my daily work. Unfortunately, it is no longer compatible with the latest versions of lme4. 

#### Problem 1
In older versions of `lme4`, fitted models were returned as objects of class `mer`; the latest versions now return fitted models as objects of class `{ngl}merMod`. There are about 1,000 places in the `XenoCat` code where functions error out if the object they are operating on is not of class `mer`. My main goal is to fix that issue so I can `XenoCat` my data again. 

#### Problem 2
All of `XenoCat`'s p-value function are broken as well. These were all based on `lme4`'s `mcmcsamp` procedure that is no longer available. Here is the relevant verbiage from the `lme4` [vignette]():

> the new version of lme4 does not provide an mcmcsamp (post-hoc MCMC sampling) method, because this was deemed to be unreliable. Alternatives for computing p-values include parametric bootstrapping (bootMer) or methods implemented in the pbkrtest package and leveraged by the lmerTest package and the Anova function in the car package (see pvalues for more details).
  
I will be implementing new [p-value calculations](http://www.phdcomics.com/comics/archive.php?comicid=905) using the [lmerTest](http://cran.r-project.org/web/packages/lmerTest/lmerTest.pdf) package. 

## Other bonuses
I will also be extending some of the graphical capabilities of the package to allow for more flexibility. Maybe I will add some other things too. I will document them all here

* Definable axis limits! 

## How can I use it?
You can install my patched `XenoCat` package straight from GitHub with the [`devtools`](https://github.com/hadley/devtools) package. Just run:
```
library(devtools) 
install_github("xeno-fix",username = "borgmaan") # install it
library(XenoCat) # loads the patched version now
```


