ChaosBook.org/extras/KSEproject/00ReadMe.txt
$Author: predrag $ $Date: 2007-04-23 21:43:07 -0400 (Mon, 23 Apr 2007)
----------------------------------------------------------------------

							24 Apr 2007

Examples of what can be done to explore KS solutions using ksfmstp solver,
download file 
	KSEproject.zip       [should always be the latest updated version]

index.m : 
	run from Matlab command line produces the figures:  

KS: 	   equilibria, 
	   travelling waves, 
	   L=22 - unstable manifold of E2 equilibirum.  

Exercises for students:
	= plot unstable manifold of another equilibrium.

----------------------------------------------------------------------

/KSEproject/ : m-files to produce figures in the /html/ directory.  

To generate html output, make sure you are in KSEproject directory, open 
	KSEproject/index.m 
in the Matlab editor, then 
	Menu -> File -> Publish To HTML.

ksfmstp.m :
	essentially Kassam and Trefethen's kursiv.m, converted into a
	function, with Jacobian computation (optional) added.

	The solution of variational equation (Jacobian) can be used 
	for computing Lyapunov exponents.

ksfm2real.m :  
	converts solution from FM to real space.

ks22equilibria.mat :
ks22statespace.mat :
	L=22 data files

gsorth.m : ?

  Ruslan L. Davidchack, Senior Lecturer in Applied Mathematics
  Address:  Department of Mathematics, University of Leicester
            University Road, Leicester LE1 7RH, United Kingdom
  Office:   210 Michael Atiyah Building
  Tel:      +44 (0)116 252-3819   Fax:  +44 (0)116 252-3915
  E-mail:   rld8@mcs.le.ac.uk     URL:
  http://www.math.le.ac.uk/people/rld8/
----------------------------------------------------------------------
							 5 may 2007

KSbiffs.tex is the start of a homepage about the class fishing trip.
On linux, converts to HTML using:

> latex KSEbiffs
> bibtex KSEbiffs
> latex KSEbiffs
> latex KSEbiffs
> tth -a -i -e2 -LKSEbiffs -v < KSEbiffs.tex > KSEbiffs.html

----------------------------------------------------------------------------

00ReadMe.txt : this file
KSEbiffs.tex : used to generate list of bifurcation figures KSEbiffs.html
srcltx.sty : used by Kile in linux (can comment out)
  
----------------------------------------------------------------------------
