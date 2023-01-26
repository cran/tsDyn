## tsDyn 11-0.4
 * By  Matthieu Stigler (2023-01-26).
 * Fixing test output as numerical differences on exotic platforms implied 11-0.3 didn't pass. 
 * Fix doi formatting

## tsDyn 11-0.3 (not submitted to CRAN)
 * By  Matthieu Stigler (2023-01-21).
 * Minor fix, integrate pull request  #47 "Replace rgl.* functions with *3d" by @dmurdoch

## tsDyn 11-0.2
 * By  Matthieu Stigler (2022-03-09).
 * Add argument `seed` in `setarTest`
 * Minor fix in using class()== ... instead of inherits (x, ...)
 * Minor fix of URL: remove URL markup as created errors. 
 * Now use only NEWS, no changelog


## tsDyn 11-0.1 (not submitted to CRAN)
 * By  Matthieu Stigler (2022-03-05).
 * Minor fix in internal test tests/setar.sim_TEST.R that showed different results under ATLAS
 * Fixing bug where TVECM.SeoTest() had not been exported


## tsDyn 11-0.0
 * By  Matthieu Stigler (2022-02-21).
 * Minor fixes for CRAN.
 * Fixing many Github issue #13, #12, #22, #23 https://github.com/MatthieuStigler/tsDyn/issues/23.
 * Improve page for coefB, (@philturk, #21).

## Version 11-0.0: Matthieu Stigler (2020-12-06)
 * Minor fixes for CRAN
 * Bunch of small fixes
 * Fixing Github issue #13, #12, #22, #23 https://github.com/MatthieuStigler/tsDyn/issues/23
 * Improve page for coefB, Github issue #21 https://github.com/MatthieuStigler/tsDyn/issues/21

## Version 10-1.2: Ho Tsung-wu (2020-02-04)
 * Minor fixes for CRAN

## Version 10-1.1: Ho Tsung-wu (2020-01-10)
 * Minor fixes for CRAN

## Version 0.9-50: Matthieu Stigler (20/5/2018)
 * new function: ECT_plot(), showing the (T)VECM response to ECT deviations
 * newfunction: ar_mean(), to compute the long-term mean of an AR process
 * MAJOR: tsDyn now depends on R 2.13, not 2.10 (i.e. at least 2011 version, not 2010)

## Version 0.9-48: Matthieu Stigler (20/5/2018)
 * minor: fix two issues in R test, reduce tolerance
 * minor: create internal coefVec
 * minor: new nlVar methods tets file
 * minor: adjust slightly the example from TVAR, TVECM

## Version 0.9-47: Matthieu Stigler (15/5/2018)   never published on CRAN
 * fix small issue in lstar code.

## Version 0.9-46: Matthieu Stigler (21/1/2018)
 * new: predict.TVAR
 * new: rewrite TVAR.gen, cleaner TVAR.boot and TVAR.sim
 * fix bug: nameB in TVAR had issue for some cases.
 * many minor updates in doc

## Version 0.9-45: Matthieu Stigler (29/8/2017)  never published on CRAN
 * minor: change registration of c code
 * minor: add dataset m.unrate, previously in (now archived) FinTs
 * minor: small documentation updates, use canonical urls, etc

## Version 0.9-44: Matthieu Stigler (22/5/2016)
 * new: vcov.TVAR
 * minor fix: TVAR_LRtest now uses robust QR instead of manual solve()...
 * fix bug: vcov.VAR corrected, was wrong
 * minor fixes:  VAR.sim did not return all args, internal insertCol did not handle well just one row, typos in man file VECM
 * minor new: summary.VAR now outputs component coefMat, more consistent output
 * minor new: TVAR outputs df.residual
 * vcov: improve vcov for VAR/VECM, using faster chol2inv saving previous qr
 * minor fixes in NAMESPACE, urls citations, etc 

## Version 0.9-43: Matthieu Stigler (23/4/2015)
 * minor fixes: add checking in selectLSTAR in case of convergence issue, rename lstar.Rout.save as created issues on other platforms

## Version 0.9-42: Matthieu Stigler (15/4/2015) never published on CRAN
 * full rewrite of VAR.sim, VAR.gen, now cleaner, allows for exogeneous predictions, etc.
 * remove predict.VECM
 * fix small issue in DESC file (pkg description in title case, only one maintainer).

## Version 0.9-41: Matthieu Stigler (26/8/2014)
 * change description for new url to github wiki, as well as new field mailinglist
 * fix suggests issues (sm::foo()) newly detected by R CMD check devel

## Version 0.9-40: Matthieu Stigler (12/8/2014), never published on cran
 * new functions coefB(), coefA(), coefPI() returning the A, B and PI coefs for VECMs
 * new function VECM.sim()
 * new function TVECM.boot()
 * new method getTh.lstar(), getTh.default()
 * new method vcov.setar()
 * new method coef.setar(), coef.lstar()
 * new method logLik.cajo.test()
 * extend argument in lineVar/VECM: allow lag=0 
 * TVECM.sim() extend for r>1, new argument seed
 * VECM(): extend ML estimation of restricted model
 * lstar(): change coefficient naming, consistency with setar()
 * minor: robustify OLS computations using lm.fit() in SETAR(), in TVAR()


## Version 0.9-33: Matthieu Stigler (28/3/2014)
 * minor: correct use of partial naming in all.equal(), which resulted in package being taken out of CRAN
 * bug correction: correct VARrep in presence of external regressors. Reported by Paul Schreck
 * bug correction: correct resVar
 * minor: removed various calls of :::
 * minor: change package call, most notably remove the enhance forecast, as forecast already in Imports

## Version 0.9-32: Matthieu Stigler (13/7/2013)
 * change some tests to avoid difference in outputs depending on platform/architecture
 * remove the data(x) used at bottom of many scripts

## Version 0.9-31: Matthieu Stigler (12/7/2013)
 * adopt half roxygen documentation
 * minor correction in rank.test: add print.summary()
 * minor correction in predict_rolling for forecast objects. 
 * changes in tests/lstar and /compareVecm to lower printing digits, to reduce changes in architecture/platform

## Version 0.9-3: Matthieu Stigler (5/5/2013)
 * change: package loading: put most packages from depends to import
 * change: TVECM.SeoTest, check with try() if inversion ok, otherwise use return NA
 * change: VAR, check with try() if inversion ok, otherwise use ginv()
 * vignette: change GIRF section
 * New argument: label in some toLatex functions
 * New argument: exogen in VECM, VAR

## Version 0.9-2: Matthieu Stigler (20/12/2012)
 * GUI: remove due to problems in assign(yy, globalEnv) not alowed anymore.
 * New function: accuracy_stat: give forecast accuracy either comparing two datasets, or giving a pred_roll object
 * Extend function: predict_rolling: allows to give rolling predictions for many models nlar, VAR/VECM, forecast

## Version 0.9-12: Matthieu Stigler (11/12/2012)
 * Correct bug: TVAR_LRtest now handles extern ThVar
 * Correct bug: TVAR had bug in naming arguments
 * Minor: many very minor changes, check class in predict, add message for ADF in VARrep, exclude lag <1 in lineVar

## Version 0.9-1: Matthieu Stigler (05/11/2012)
 * Extend method: vcov.lstar
 * Extend method: confint.lstar
 * Extend method: getTh.lstar
 * Extend method: regime.lstar
 * New argument: thInt as control for starting values of lstar()
 * New argument: trace in selectLSTAR and selectNNET
 * New argument: strategy for selectLSTAR: allows recovering previous optimal values for starting values in all LSTAR searches. 
 * Minor: changed procedure in lstar(), recover only parameters with optimHessian, instead of runnign second optim round. 
 * Correct bug: corrected bug in exporting fake ca.jo object, reported by Uwe Ligges. 


## Version 0.9-0: Matthieu Stigler (26/10/2012)
 * New function: rank.select() and derived lags.select(): compute the rank and lags according to AC/BIC criteria. 
 * New function: rank.test() for Johansen LR rank test, using gamma approximated p-values. 
 * New function: fevd() and irf() from package vars are now available. 
 * New function: predict() now implemented for linear multivariate models. 
 * New function: predict_rolling(): allow rolling forecasts, using sub-samples. 
 * New function: VARrep(): VAR representation of a VECM, or of a VAR in differences. 
 * New function: VECM.sim(): wrapper of TVECM.sim()
 * New function: VAR.sim(): wrapper of renamed TVAR.gen()
 * New function: VAR.boot(): wrapper of renamed TVAR.gen()
 * New argument: fitted() can be retunred either as in the model (so sometimes in diff) or on the original series with arg "level". 
 * New argument: VAR.sim() can now fix seed when generating/resamping innovations, with arg "seed".  
 * New argument: predict() now can predict using bootstrap or MC, with arg "type". 
 * New argument: lineVar() can compute also model with type="ADF" as in univariate case. 
 * New argument: AIC, BIC are corrected and have new argument fitMeasure. 
 * New help page: document lokLik()
 * New test page: nlar-methods. 
 * New test page: VAR. 
 * Correct bug: Correct bug in resVar()
 * minor: TVECM.sim(): print more information on input coefficient matrix when demanded. 
 * minor: plot.setar(): change lty
 * minor: adapt citation to new person() function. 
 * minor: adapt documentation pages for cross-links to other packages using \pgk{}
 * minor: corrections in VECM() for restricted constant and intercepts in MLE. 
 * minor: corrections in VECM() for restricted constant and intercepts in OLS. 
 * minor: getTh() now alos works for linear multivariate models. 

## Version 0.8-1: Matthieu Stigler (13/02/2012)
 * minor: avoid partial matching in call to legend() in TVECM.HStest
 * minor: re-run tests with newest R 2.14 in english
 * minor: move vignettes from /inst/doc to new folder vignettes/

## Version 0.8: Matthieu Stigler (12/02/2012)
 * New feature: Add argument include in lstar: possibility to estimate lstar with constant, none, trend or both.
 * New feature: Add argument starting.control in lstar: possibility to specify grid search parameters
 * New feature: Add argument hpc in selectSETAR to run on parallel cores. 
 * Minor: add explicit warning in lstar when grid search select bound. 
 * Minor: add explicit warning in lstar when optim code is 1, referencing argument maxit. 
 * Minor: Cleaning of help files, correcting many typos. 

## Version 0.7-62: Matthieu Stigler (19/12/2011)
 * Add standard errors for lstar()
 * Minor: change sd(matrix) to apply(matrix,2,sd) (chage in R 2.14)
 * Minor: change partial argument matching in a few calls
 * Minor: correct vignette ThCointegration

## Version 0.7-61: Matthieu Stigler (02/08/2011)
 * TVECM: allow for K>2 variables when the cointegrating vector is pre-specified
 * TVECM summary: change output, adding newlines after coint vector and percentage of observations

## Version 0.7-6: Matthieu Stigler (14/05/2011)
 * correct bug in star (bugreport by Petri H) by Christoph B
 * change return value of star() in case of no star regime detected: returns linear object instead of series itself

## Version 0.7-52: Matthieu Stigler (08/05/2011)
 * add tests for lstar
 * correct lstar so that optim function re-evaluate each time linear parameters
 * remove VECM.Rout as useless

## Version 0.7-51: Matthieu Stigler 
 * update NAMESPACE as previous (0.7-50) did not work
 * add VECM.Rout
 * add econometrica .bst and amsplain.bst

## Version 0.7-50: Matthieu Stigler (30/03/2011)
 * remove own BIC function and use that of pkg nlme, as recommend by Prof Ripley
 * fix bug in Rd files, and description (missing suggest packages)
 * fix bug in 'autopairs' and 'autotriples': bogus 'showvalue' option was passed to the 'scale' widget
 (bugreport by Budi Hari Priyanto)

## Version 0.7-40: Matthieu Stigler (11/08/2010)
 * Add paralell option for TVAR.LRtest(), setarTest(), TVECM.Seotest()
 * regime() works also for nlVar, also add arguments initVal and timeAttr
 * RENAME function TVECM.HSTest in TVECM.HStest
 * correct bug in TVECM.HSTest

## Version 0.7-30: Matthieu Stigler (30/06/2010)
 * Add function TVECM.HStest (before was HanSeo_TVECM). Has new argument hpc for parallel computing.
 * Add function VECM(), simple wrapper for lineVar(..., model="VECM")
 * Add Johansen estimator for VECM()
 * Add methods AIC, BIC, logLik for class VECM
 * Correct numerous bugs

## Version 0.7-23: Matthieu Stigler (01/04/2010)
 * correct bug TVAR.sim
 * export function VAR.sim, simpler wrapper for TVAR.sim

## Version 0.7-22: Matthieu Stigler (22/02/2010)
 * correct bug in HanSeo_TVECM, which can now be used (but still not exported)

## Version 0.7-2: Matthieu Stigler (20/02/2010)
 * Correct degrees of freedom for VAR covariance matrix
 * Standardize output of lineVar with lm
 * update slightly vignette
 * add specific vignette on threshold cointegration

## Version 0.7-1: Matthieu Stigler (04/09/2009)
 * Change args of TVECM from nn to ngridBeta, from ngridG to ngridTh, model to common
 * Add functions regime to extract variabl indicating state regime
 * Add function getTh to extract threshold values
 * Adapt toLatex.setar to handle MTAR models. Still needs improvement
 * Fix bugs in TVECM, in HansSeo_TVECM and in setar

## Version 0.7-0: Matthieu Stigler (02/07/2009)
 * added possibilty to have two thresolds and hence three regimes in setar and selectSETAR (arg nthresh)
 * created new functions for unit roots tests: KapShinTest() and BBCTest()
 * created new functions for estimating VAR and VECM: lineVar
 * created new function for estimating TVECM: function TVECM()
 * created new function for estimating TVAR: function TVAR()
 * created new function to test for setar: function setarTest()
 * created new function to test for TVAR: function TVAR.LRtest()
 * created new function to test for TVECM: functions TVECM.SeoTest() and HanSeo_TVECM()
 * created new function to simulate/bootstrap a TVAR: function TVAR.sim()
 * created new function to simulate/bootstrap a TVECM: function TVECM.sim()
 * created new function to simulate/bootstrap a setar: function setar.sim()
 * created new function to estimate regime-specific variance in setar: function resVar()
 * created new function to extend a bootstrap replication in setarTest: function extendBoot()
 * added 'zeroyld' dataset
 * added 'unemp' dataset
 * added 'IIPUs' dataset
 * added in selectSETAR() and setar() following args: include, common, model, trim, MM, ML, MH, model, restriction
 * added MakeThSpec() function for specification of the threshold in the search
 * added in selectSETAR(): criterion "SSR" (sum of squares residual) and argument max.iter 
 * extended arg th in selectSETAR to search inside an intervall or around a point or on the whole grid
 * added k parameter to the AIC method and fixed a bug
 * standardized license naming in description
 * fixed typo in print.star

## Version 0.6-1: Antonio Fabio di Narzo
 * fixed CITATION file to the new R format
 * fixed bug in nnetTs: control args were not passed to 'nnet' fitter
	reported by prof. Joachim Zietz
 * fixed bug in llar.predict (spotted by Markus Gross)

## Version 0.6-0: Antonio Fabio di Narzo
 * fixed Matrix dependency in DESCRIPTION
 * added minimal docs for newly created internal objects
 * added STAR model option to the tcltk GUI
 * added a new STAR-fitting function with automatic selection of the number of regimes
 * fixed bugs in lstar
 * fixed bugs in some internal functions
 * added google codebase project URL to DESCRIPTION

## Version 0.5-6: Antonio Fabio di Narzo
 * fixed bug in setar model when mL and/or mH=0

## Version 0.5-5 released
 * fixed bug in lstar model out-of-sample forecasting

## Version 0.5-4 released
 * fixed some warnings issued in R-2.4.0 check script

## Version 0.5-3 released
 * fixed some warnings issued in R-devel check script
 * added URL field to DESCRIPTION

## Version 0.5-2 released
 * fixed gui design error in autopairs and autotriples

## Version 0.5-1 released
 * added draft of internal design vignette
 * minor fixes in vignette, added AAR model fitting example
 * fixed bug in setar and lstar (introduced in version 0.5)

## Version 0.5 released
 * added minor fixes for passing checks under R-2.3.0
 * added tsDyn package overview man page, with getting started indications
 * upgraded vignette to the new models usage
 * fixed bug in plot.llar
 * changed llar interface to nlar models schema
 * *eliminated* nlar function: now each model has its own function. Using more explicitely classes inheritance and polymorphism. Adding more nlar models now should be much simpler.

## Version 0.4-1 released
 * updated all package dialogs to the new GUI system
 * written a whole (basic) object oriented GUI system, with pass-by-reference semantics
 * added input checking to first dialog of nlarDialog()
	
## Version 0.4 released
 * added hidden 'series' name argument to 'nlar'
 * added early version of nlarDialog, a GUI for 'nlar' model fitting
 * fixed latex apparence of some manual pages
 * fixed notational errors in vignette: Z_t=Z_{t-d}, z_t changed in Z_t in LSTAR model formulation
 * added nlar.fit class documentation
 * changed text entry with slide bar in autopairs and autotriples	
 * added availableModels function, upgraded docs	
 * added and documented aar model	
 * added and documented tarch model	
 * added and documented arch model
	
## Version 0.3-5 released
 * largely improved apparence of regime-switching plot, changed palette for lstar diagnostic plots
 * added latex printing of fitted setar models
	
## Version 0.3-4 released
 * changed lag-plots for setar and lstar models when more than 300 points are to be plotted
 * fixed notation error in vignette (equation (4) on neural networks)
 * deleted 0-lag in acf, pacf and ami summary plots for nlar.fit objects
 * extended threshold variable printing informations for setar and lstar models
 * added alternative selection criteria for selectSETAR
 * added regime proportions info in print.setar and summary.setar
 * added regime switching plot to plot.setar and plot.lstar
 * fixed graphical pars. resetting in plot.nlar.fit
 * added convenience wrappers 'setar' and 'lstar', improved relative documentation
	
## Version 0.3-3 released
 * fixed minor bug in selectSETAR output presentation
	
## Version 0.3-2 released
 * changed general AIC formula
 * added dummy (non-working), undocumented version of FAR (Functional AutoRegressive) model
 * added early, undocumented version of AAR (Additive AutoRegressive) model
 * modified selectSETAR, selectLSTAR, selectNNET functions: now more coherent results
	
## Version 0.3-1 released
 * added autotriples.rgl (now also the rgl package is suggested)
 * added the possibility to stop the llar procedure
 * fixed y-scale in llar diagnostic plot
 * improved autotriples, upgraded vignette 
	
## Version 0.3 released
 * added llar.predict and llar.fitted functions, improved llar output object

## Version 0.2-4 released
 * fixed vignette
 * fixed external variable management in predict.setar and predict.lstar

## Version 0.2-3 released
 * added vignette
 * improved SETAR and LSTAR models: now 'predict' interface is sensible to the type of threshold variable
 * fixed and exported delta, delta.test, delta.lin and delta.lin.test

## Version 0.2-2 released
 * added optional legend to plot.setar and plot.lstar

## Version 0.2-1 released
 * added autopairs and autotriples functions, with a minimal (optional) TclTk GUI

## Version 0.2 released
 * changed builtin nonlinear models documentation structure: now one page per model
 * added possibility of external threshold variables in setar models
 * added unexported delta tests of indipendence and linearity
 * added hidden tests in nlar-methods and llar man pages examples
 * improved summary method for linear and setar models
 * added llar routine, i.e., Casdagli test of nonlinearity

## Version 0.1-2 released
 * removed unuseful sum-of-squares functions
 * added selectSETAR, selectLSTAR and selectNNET functions, for automatic selection of model hyper-parameters
 * minor bug fixes
 * added 'predict' method to nlar.fit objects, subclasses linear, setar, lstar and nnet
 * fixed nlar.fit objects time series storage redundancy

## Version 0.1-1 released
 * added neural network model documentation
 * added dummy predict method to nlar.fit
 * minor bug fix in lstar model
 * added lstar and setar plot method
 * fixed lstar model documentation
 * bug fix in summary.lstar when embedding dimension is 1
 * bug fix in nlar and plot.setar when embedding dimension is 1
 * added dummy 'summary', 'print' and 'plot' function to nnet subclass
 * added utility 'sum-of-squares' functions for implemented models
 * added plot method to nlar.fit class

## Version  0.1: initial release
