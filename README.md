mete-spatial
===========

Code for working with the spatial predictions of John Harte and colleagues' Maximum Entropy Theory of Ecology.

The code in this project is for two projects, one that investigates the species-area relationship (SAR) and another that investigates the distance decay relationship
(DDR). 

Species-area Project
--------------------

The SAR code provides the supplemental material to: 

McGlinn, D.J., X. Xiao, E.P. White. 2013. An empirical evaluation of four
variants of a universal species-area relationship. PeerJ. 1: e212. (http://dx.doi.org/10.7717/peerj.212)


### Setup

Requirements: R ≥ 2.7.0 with the R package RCurl installed. Python 2.x and the following Python modules: numpy, scipy, matplotlib, mpmath. You will also need two of our custom Python modules: METE (https://github.com/weecology/METE) and macroecotools (https://github.com/weecology/macroecotools).
These modules can be installed by running the following commands from the command
line (with sufficient permissions):

```sh
git clone https://github.com/weecology/METE.git
cd METE
python setup.py install
git clone https://github.com/weecology/macroecotools.git
cd macroecotools
python setup.py install
```

### Replicate analyses

The analyses can be replicated by running the following commands from the
command line.

Run all analyses and generate figures:
`Rscript sar_run_all.R`

Please note that these analyses involve both a large amount of data and a lot of
computational work and therefore take a while to run. Expect the empirical
analysis to take a few hours. 

Distance-decay Project
----------------------

The DDR code provides the supplemental material to:  

McGlinn, D.J., X. Xiao, J. Kitizes E.P. White. in prep. Exploring the spatially explicit predictions of the Maximum Entropy Theory of Ecology

### Setup

Requirements: R ≥ 2.7.0 with the R packages vegan, RCurl, and bigmemory installed. Python 2.x and the following Python modules: numpy, scipy, matplotlib, mpmath. You will also need two of our custom Python modules: METE (https://github.com/weecology/METE) and macroecotools (https://github.com/weecology/macroecotools).
These modules can be installed by running the following commands from the command
line (with sufficient permissions)

### Replicate analyses

The analyses can be replicated by running the following commands from the
command line.

Run all analyses and generate figures:
`Rscript ddr_run_all.R`

Please note that these analyses involve both a large amount of data and a lot of
computational work and therefore take a while to run. Expect the empirical
analysis to take a few hours. 


License
-------
mete-spatial is licensed under the open source [MIT License](http://opensource.org/licenses/MIT)

Copyright (c) 2012 Daniel McGlinn

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Acknowledgments
---------------
Development of this software was funded by the [National Science Foundation](http://nsf.gov/) as part of a [CAREER award to Ethan White](http://www.nsf.gov/awardsearch/showAward?AWD_ID=0953694)
