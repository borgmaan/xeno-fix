# xeno-fix
This is a clone of the [`XenoCat`](https://code.google.com/p/r-xenocat/) R package, updated to work with the latest versions of the [`lme4`](https://github.com/lme4/lme4) package. `XenoCat` is a great R package for analyzing tumor growth experiments -- stats described [here](http://clincancerres.aacrjournals.org/content/early/2012/07/19/1078-0432.CCR-11-3215.full.pdf). It is much better than running a t-test at some bogus timepoint you chose before your xenograft experiment started. 

## What's going on?
I really like the `XenoCat` package and use it regularly in my daily work. Unfortunately, it is no longer compatible with the latest versions of lme4. In older versions of `lme4`, fitted models were returned as objects of class `mer`; the latest versions now return fitted models as objects of class `{ngl}merMod`. There are about 1,000 places in the `XenoCat` code where functions error out if the object they are operating on is not of class `mer`. My main goal is to fix that issue so I can `XenoCat` my data again. 

## Other bonuses
I will also be extending some of the graphical capabilities of the package to allow for more flexibility. Maybe I will add some other things too. I will document them all here

* Great feature 1...
* Bad feature 2...
* Useless feature 3...

## How can I use it?
You can install my patched `XenoCat` package straight from GitHub with the [`devtools`](https://github.com/hadley/devtools) package. Just run:
```
library(devtools)
install_github("xeno-fix",username = "borgmaan")
```

