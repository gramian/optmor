# optmor

OPTimization-based Model Order Reduction (Version 2.5) ![DOI badge](https://zenodo.org/badge/doi/10.5281/zenodo.46683.png)


## About

optimization-based model order reduction,
for the computation of combined state and parameter
reduced order models of state-space input-output systems


## License

All source code is licensed under the open source 
[BSD 2-clause license](http://opensource.org/licenses/BSD-2-Clause):

#### Copyright (c) 2013-2016, Christian Himpe

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.</p>


## Usage

General Usage: `{X,P} = optmor(f,g,s,t,w,pr,nf,ut,us,xs,um,xm);`

Minimal Usage: `{X,P} = optmor(f,g,s,t,w);`

About Info Usage: `v = optmor('version');`


## Arguments

`f` - system function handle; signature: xdot = f(x,u,p)

`g` - output function handle; signature: y = g(x,u,p)

`s` - system dimensions [inputs,states,outputs]

`t` - time discretization [step,stop]

`r` - reduced order(>1) or error threshold(<1)

`q` - nominal parameter


### Optional Arguments

`nf` - options, 6 components

`ut` - input; default: delta impulse

`x0` - initial state; default: zero

`co` - covariance matrix; default: unit 

`yd` - experimental data; default: empty


## Flags

`nf(1)` - Optimization Algorithm: fminunc(0), fminsearch(1), custom(-1)

`nf(2)` - Lasso Regularization Weight: default(0)

`nf(3)` - Tikhonov Regularization Weight: default(0.1)

`nf(4)` - Data-Driven Regularization Weight: default(0)

`nf(5)` - Number of Maximum Optimizer Iterations: default(4)

`nf(6)` - Initial Parameter: last(0), random(1)


## References

* C. Himpe, M. Ohlberger; "[Data-driven combined state and parameter reduction for inverse problems](http://dx.doi.org/10.1007/s10444-015-9420-5)"; Advances in Computational Mathematics 41(5): 1343--1364, 2015.
    + One-Sentence Abstract: Extensions for combined state and parameter reduction for inverse problems.


## Cite

* Cite as:

    C. Himpe (2016). optmor - Optimization-Based Model Order Reduction (Version 2.5) [Software]. http://gramian.github.io/optmor

* BibTeX:

        @MISC{optmor,
         author = {C. Himpe},
         title = {{optmor} - Optimization-Based Model Order Reduction (Version 2.5)},
         howpublished = {\url{http://gramian.github.io/optmor}},
         year = {2016}
        }
* DOI: [10.5281/zenodo.46683](http://dx.doi.org/10.5281/zenodo.46683) (Version 2.5)
* Last Change: 2016-02-28


## Links

* http://gramian.de emgr - Empirical Gramian Framework
* http://modelreduction.org Model Order Reduction Wiki
* http://morepas.org Model Reduction for Parametrized Systems

