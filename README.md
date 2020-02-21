JAMM Peak Finder
======

**JAMM is a peak finder for NGS datasets (ChIP-Seq, ATAC-Seq, DNase-Seq..etc.) that can integrate replicates and assign peak boundaries accurately. JAMM is applicable to both broad and narrow datasets. Read more [in JAMM's publication](https://academic.oup.com/bioinformatics/article/31/1/48/2365061).** 

JAMM was developed at the [Ohler lab](http://ohlerlab.mdc-berlin.net/) at [BIMSB-MDC](https://www.mdc-berlin.de/bimsb) in Berlin.

If you have questions or need help running JAMM please email us at [this email](http://scr.im/jammpro), we will be happy to help you.



Download JAMM
------

**The latest JAMM version is _JAMMv1.0.7.6_. [Click here to download it](https://github.com/mahmoudibrahim/JAMM/archive/JAMMv1.0.7.6.zip).**

If you are interested in older versions, check the [Downloads Archive](https://github.com/mahmoudibrahim/JAMM/wiki/JAMM-Downloads-Archive) or the [Github release page](https://github.com/mahmoudibrahim/JAMM/releases). 



Latest News and Updates
------

* **Feb. 21 2020:** *JAMM will be updated incrementally over the next few months under version 1.0.8.x* to improve performance and speed. These updates are stable, synced directly to the master repository and released as pre-releases. However, they may change output somewhat significantly, so use at your risk. The most recent stable and tested version is [v1.0.7.6](https://github.com/mahmoudibrahim/JAMM/archive/JAMMv1.0.7.6.zip) which you can get from the [releases page](https://github.com/mahmoudibrahim/JAMM/releases).

* **May 17 2019:** *New version of JAMM released (v1.0.7rev6)*. This version features some small changes to input checks.

_Visit the [the Wiki Homepage](https://github.com/mahmoudibrahim/JAMM/wiki) for older updates and news._



JAMM Documentation
------

 * [How to Install JAMM](https://github.com/mahmoudibrahim/JAMM/wiki/Installing-JAMM)
 * [How to Use JAMM](https://github.com/mahmoudibrahim/JAMM/wiki/Usage)
 * [JAMM's Output](https://github.com/mahmoudibrahim/JAMM/wiki/JAMM%27s-Output)
 * [Frequently Asked Questions](https://github.com/mahmoudibrahim/JAMM/wiki/Frequently-Asked-Questions)
 * [Parameter Recommendations](https://github.com/mahmoudibrahim/JAMM/wiki/JAMM-Parameter-Recommendation)
 * Of course it is always good to read [the actual paper](https://academic.oup.com/bioinformatics/article/31/1/48/2365061) :)

Please check [the Wiki Homepage](https://github.com/mahmoudibrahim/JAMM/wiki) for more information on JAMM. 

If you still have questions or need help running JAMM please email us at [this email](http://www.google.com/recaptcha/mailhide/d?k=01vPijd2GG0LEbZV2NyE_rSA==&c=49GEiFp47dZQV_120clczwxYcP3tQ98VWBJtNl6_dBw=).

Note: JAMM produces a large list of peaks on purpose so that you can choose how to threshold it. If you want JAMM to try and automatically threshold it, use "-e auto" option. Read more about this in the [FAQs on the wiki](https://github.com/mahmoudibrahim/JAMM/wiki/Frequently-Asked-Questions).


Example Dataset
------

*You can use [this small dataset](https://drive.google.com/file/d/0B8nxBVNVchN9cFFzQnQxMnNQUjQ/edit?usp=sharing) to test if JAMM is installed correctly. Installation instructions above!*

This dataset is chromosome 21, in K562 CTCF ChIP-Seq from [ENCODE](https://genome.ucsc.edu/ENCODE/).



Quick Tutorial
------

On [this link](https://github.com/mahmoudibrahim/JAMM/wiki/Basic-Tutorial), you will find a quick tutorial to show how you can run JAMM on a real dataset.


Galaxy Wrapper
------
JAMM has a Galaxy wrapper which you can get from here: [https://github.com/mahmoudibrahim/JAMMalaxy](https://github.com/mahmoudibrahim/JAMMalaxy)

The Galaxy wrapper was developed by [Clemens Messerschmidt](https://github.com/messersc)


How to Cite JAMM?
------

Please cite:

**Ibrahim MM, Lacadie SA, Ohler U (2015). JAMM: A Peak Finder for Joint Analysis of NGS Replicates. _Bioinformatics_, 31(1): 48-55. doi: [10.1093/bioinformatics/btu568] (https://academic.oup.com/bioinformatics/article/31/1/48/23650618)**


---

*Hint: While you run JAMM, you can also go on jammin...*
[![IMAGE ALT TEXT HERE](http://img.youtube.com/vi/HSs1HgM0Wos/0.jpg)](http://www.youtube.com/watch?v=HSs1HgM0Wos)
