These are the deprojected, azimuthally-averaged radial profiles of the 
continuum surface brightnesses for the DSHARP targets.  These profiles have 
been computed as described in Section 3.1 of DSHARP II: Characteristics of 
Annular Substructures by Huang et al. (2018, ApJL, XXX, XXX).  

From left to right, the columns are:

1) Radius in au.  This assumes the distances in Table 1 of DSHARP I.

2) Radius in arcseconds.

3) Surface brightness (Inu) in units of mJy / beam, (Jansky / beam).  Beam dimensions are in 
   Table 4 of DSHARP I (more precisely, in the BMAJ and BMIN header variables 
   in the fiducial continuum FITS files provided in this Data Release).

4) 1-sigma uncertainties on 3), again in units of mJy / beam.

5) Brightness temperatures in units of K, calculated assuming the 
   Rayleigh-Jeans approximation is valid; i.e., Tb where

        Inu = 2 nu**2 k Tb / c**2

   with nu (the frequency) specified in Table 4 of DSHARP I (more precisely, in 
   the CRVAL3 header variables in the fiducial continuum FITS files provided in 
   this Data Release).  Note that this requires the conversion of the Inu units 
   in 3) into something more physically meaningful.  It is useful to know that 
   the area of a Gaussian beam with bmaj and bmin the FWHM of the major and 
   minor beam axes is 

        pi * bmaj * bmin / (4 ln 2)

6) 1-sigma uncertainties on 5), again in units of K.

7) Brightness temperatures in units of K, in this case using the full Planck 
   equation; i.e., Tb where

        Inu = (2 h nu**3 / c**2) / (exp(h*nu / k*Tb) - 1)

   and unit conversions as discussed in the description of 5).  This is the 
   more appropriate Tb to use in comparisons to dust temperatures (and thereby 
   optical depths) if that is of interest, since the mm part of the spectrum 
   for typical disk temperatures is decidedly not in the R-J regime.  Note that 
   the Tb calculation in this case has a singularity when Inu <= 0 (i.e., in 
   the noisy regions at large radius).  In those cases, Tb = NAN.

8) A crude estimate of the 1-sigma uncertainties on 7), again in units of K.  
   The general approach here is to take the average of the deviations in Tb 
   from the mean in 7) induced by the 1-sigma uncertainties on Inu in 4).  [The 
   average because in general the nonlinearity of the Tb definition means that 
   the proper pdf for Tb is asymmetric if the pdf for Inu is Gaussian and 
   symmetric; that asymmetry is generally small here, so we ignore it.]  But 
   again, in cases where the lower uncertainty bound on Inu is <= 0 (i.e., in 
   the noisy outer disk, or even in deep gaps) the corresponding singularity 
   poses a problem.  There we simply adopt a symmetric uncertainty based on the 
   upper bound uncertainty instead.  There are certainly different (better) 
   approaches for how to deal with these issues; users should feel free to 
   adopt a different set of assumptions by re-calculating the Tb curves from 
   the Inu curves (or the images directly).


*** Special notes on AS 205 and HT Lup:
    These profiles include the contributions from companion disk(s).  They 
    were constructed using the procedure described in DSHARP II, assuming the 
    geometries derived in the spiral models in DSHARP IV (Kurtovic et al. 
    2018).  They are centered on the primary disk peaks.
