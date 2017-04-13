wrightOmegaq
========
##### Complex double-precision evaluation of the Wright omega function, a solution of ```W+LOG(W) = Z```.
###### Version 1.0, 3-12-13
##### Download Repository: [ZIP Archive](https://github.com/horchler/wrightOmegaq/archive/master.zip)

--------

[```wrightOmegaq(Z)```](https://github.com/horchler/wrightOmegaq/blob/master/wrightOmegaq.m) performs floating point evaluation of the Wright omega function. ```Z``` is an array and may be complex. If ```Z``` is an array of symbolic values, it is converted to double-precision for computation and then recast as symbolic.

```wrightOmegaq``` is up three to four orders of magnitude faster than ```wrightOmega``` for double-precision arrays. Additionally, it has much less numeric error, properly evaluates values greater than ```2^28```, supports single-precision evaluation, and handles ```NaN``` inputs.
&nbsp;  

--------

References  

 1. P.W. Lawrence, R.M. Corless, and D.J. Jeffrey, &#8220;Algorithm 917: Complex Double-Precision Evaluation of the Wright omega Function,&#8221; *ACM Transactions on Mathematical Software*, Vol. 38, No. 3, Article 20, pp. 1&ndash;17, Apr. 2012. [[http://dx.doi.org/10.1145/2168773.2168779](http://dx.doi.org/10.1145/2168773.2168779)]
 2. R.M. Corless and D.J. Jeffrey, &#8220;The Wright omega Function,&#8221; *In: Artificial Intelligence, Automated Reasoning, and Symbolic Computation*, Joint International Conferences, AISC 2002 and Calculemus 2002, Marseille, France, July 2002, (Jacques Calmet, Belaid Benhamou, Olga Caprotti, Laurent Henocque, and Volker Sorge, Eds.), Berlin: Springer-Verlag, pp. 76&ndash;89, 2002. [[http://orcca.on.ca/TechReports/2000/TR-00-12.html](http://orcca.on.ca/TechReports/2000/TR-00-12.html)]
&nbsp;  

--------

Andrew D. Horchler, *horchler @ gmail . com*, [biorobots.case.edu](http://biorobots.case.edu/)  
Created: 7-12-12, Revision: 1.0, 3-12-13  

This version tested with Matlab 9.0.0.341360 (R2016a)  
Mac OS X 10.11.4 (Build: 15E65), Java 1.7.0_75-b13  
Compatibility maintained back through Matlab 7.4 (R2007a)  
&nbsp;  

--------

Acknowledgment of support: This material is based upon work supported by the [National Science Foundation](http://www.nsf.gov/) under [Grant No.&nbsp;1065489](http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=1065489). Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.  
&nbsp;  

Copyright &copy; 2012&ndash;2017, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.