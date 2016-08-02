
Implementation of ZOGY method v1
--------------------------------

** Note: the input sci (new) and ref images will already have been gain-matched,
   background-matched, registered astrometrically, and with the ref-image
   resampled (using a Lanczos kernel) onto the sci-image frame. 

--
For this first version, all we need is a difference image (D), a 
matched-filtered score image for D (Scorr) and its PSF (P_D). IPAC
will handle the detection step using Scorr and the flux estimation
step using D and P_D.

--
Contents of this directory:
- testsci_newscibmtch.fits: example prepped new (sci) image. 
- testsci_resampref.fits: accompanying prepped ref image.
- sci_psf.fits: PSF template for new (sci) image.
- ref_psf.fits: PSF template for ref image.
- testsci_pmtchscimref.fits: diff-image from PTFIDE for comparison.

--
Some input parameters pertaining to the above example:
- effective electronic gain in both new and ref images: 1.539 e-/DN 
- background level in new image (mode): 303.922 DN 
- pixel noise-sigma of background in new image (mad): 52.764 DN 
- background level in ref image (mode): 309.590 DN
- pixel noise-sigma of background in ref image (mad): 2.259 DN 

--
Steps for first version of code that implements ZOGY method:

Equation numbers below refer to those in:
http://arxiv.org/pdf/1601.02655v2.pdf
Definitions of all inputs are summarized on page 10.

(1) compute FFT of the difference image D from FFTs of the inputs
    using eq. 13. The flux-based gain factors Fr, Fn can be set to 1
    since the input images are already gain-mached.

(2) compute inverse FFT of output from (1) to obtain the real-space
    difference image D.

(3) compute FFT of the PSF image for D from FFTs of the inputs using eq. 14. 

(4) compute inverse FFT of output from (3) to obtain the real-space PSF 
    image P_D.

(5) compute the point-source match-filtered image S using eq. 16.

(6) compute the score image (or effective signal-to-noise ratio image) Scorr
    from the ratio "S/sigma_S" where S is the output from (5) and sigma_S is
    an image of the 1-sigma pixel uncertainties in S. I.e., sigma_S effectively 
    includes background noise (and implicitly readnoise) and any excess
    photon-noise associated with signals above the background. The sigma_S
    image will also include the effects of smoothing from P_D that created S. 
    I suggest for now only accounting for photometric noise when computing
    sigma_S, i.e., sigma_S = sqrt[V(Sn) + V(Sr)] where V(Sn) and V(Sr)  
    can be computed using equations 26, 27, 28, 29.

-- END --
