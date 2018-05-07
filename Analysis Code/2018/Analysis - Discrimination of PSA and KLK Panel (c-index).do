/*
PROGRAM: Analysis - Discrimination of PSA and KLK Panel (c-index)
PROGRAMMER: Daniel
PURPOSE: 
The following code copmutes the c-index for various models: PSA alone, free to total ratio alone, 
total+PSA+free PSA+age, and the KLK panel+age.
The discrimination index is calculted in variables PSA subgroups (e.g. PSA>=1, PSA<1, PSA>50th centile).  
The dataset contains imputed data (10-times).  
1. The first section of code calculates the discimrination
2. The second combines the results accross the imputations.
3. And the third prepares the results to be displayed.
*/

cd "O:\Outcomes\Andrew\Analytic Projects\1 Active\Lilja MDC KLK predicts PCa Outcomes"
use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear

****************************************
**  Section 1: CALCULATE THE C-INDEX  **
****************************************
* this porgmam calculated Harrell's c-index on the imputed data	
capture program drop cindexmi
program cindexmi
	syntax , 	cohort(string) /// specifes the subset of patients to perform caluclations on (typically age subsets)
				cuttype(string) /// either pselevel of psacentile
				cutsymbol(string) /// >= or < (aids in deifning PSA subset we're intrested in, e.g. tpsa >= 2)
				cutlist(numlist) /// list of values to use to create PSA subsets, e.g. tpsa >= 2, tpsa >= 3
				outcome(string)  /// varname of the outcome
				covar(string) /// predictor varname
				[saving(string)] 

	*setting up empty dataset to append estimates
	tempname memhold
	tempfile results
	postfile `memhold' str50(outcome covar cuttype cutsymbol cohort) cutraw psacut psactlcut m n casen cindex str15(error) using `results'
		
	noi display `"`outcome' `cohort' $S_DATE $S_TIME"'
	qui foreach cut of numlist `cutlist' {
		noi display "`cuttype'`cutsymbol'`cut' $S_DATE $S_TIME"
		
		*looping over each of the imputed datasets
		qui foreach m of numlist 1/10 {
			use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear
			* using mth imputed set of data
			mi extract `m', clear 
			noi display "    m=`m'"
			keep if `cohort'

			*creating local to store range of PSA values to keep
				* PSA LEVEL is a known cutpoint for PSA
				if "`cuttype'"=="psalevel" {
					local psacut=`cut'
					*getting the centile
					count if tpsa<`psacut'
						local nbelow=`r(N)'
					count if tpsa<=`psacut'
						local nabove=`r(N)'
					local ctlcut=((`nabove' + `nbelow')/2)/`=_N'	
				}
				* cutpoint for a PSA percentile
				else if "`cuttype'"=="psacentile" {
					local ctlcut=`cut'/100
					centile tpsa, c(`cut')
						local psacut=`r(c_1)'
				}
				* print error if misspecified
				else {
					disp as err "cuttype must be psalevel psacentile"
					exit 100
				}

			*keeping subset for analyses
				keep if tpsa `cutsymbol' `psacut'
			
			*prepping the data for survival analyses
				stset tt`outcome', f(`outcome')
				
			*getting counts for output table
				count 
					local N=`r(N)'
				count if `outcome'==1 
					local caseN=`r(N)'

			*as there are so many imputed datasets, ensuring that there is enough data to calculate the c-index
				if `caseN'>=7 & `caseN'<`N' {
					* calculating the c-index and posting results
					stcox `covar' 
					estat concord
					post `memhold' ("`outcome'") ("`covar'") ("`cuttype'") ("`cutsymbol'") ///
									 ("`cohort'") (`cut') (`psacut') (`ctlcut') (`m') (`N') (`caseN') (`r(C)') ("")	
				}
			*if not enough data, posting a missing with note why c-index coulud not be calculated
				else if !(`caseN'>=7) post `memhold' ("`outcome'") ("`covar'") ("`cuttype'") ("`cutsymbol'") ///
									 ("`cohort'") (`cut') (`psacut') (`ctlcut') (`m') (`N') (`caseN') (.) ("Too Few Events")	
				else if !(`caseN'<`N') post `memhold' ("`outcome'") ("`covar'") ("`cuttype'") ("`cutsymbol'") ///
									 ("`cohort'") (`cut') (`psacut') (`ctlcut') (`m') (`N') (`caseN') (.) ("No Controls in Cohort")	
		}
	}
	
	* loading results into memory
	postclose `memhold'
	use `results', clear
	
	*saving results
		if "`saving'"!="" save "`saving'", replace
end

						
*  HERE we loop through age groups (defined in te variable cohort), and the various models and calculate the c-index
*  The c-=index is calculated for both cancer diagnosis (pca) and death from cancer (dod)
	x
	use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear
	levelsof cohort, local(cohortlevels)
	set more off
	foreach outcome in /*pca*/ dod {
		foreach cohortn in `cohortlevels' {
			foreach covar in tpsa protectriskhg ftpsa ftpsariskhg {
					* above cuts (e.g. PSA>=2)
					cindexmi, cohort(cohort==`cohortn') cuttype(psalevel) cutsymbol(>=) cutlist(0 1 1.5 2 3) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' above psalevel.dta")

					* below cuts (e.g. PSA<1)
					cindexmi, cohort(cohort==`cohortn') cuttype(psalevel) cutsymbol(<) cutlist(1) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' below psalevel.dta")

					* above centile (e.g. PSA>= 50th percentile)
					cindexmi, cohort(cohort==`cohortn') cuttype(psacentile) cutsymbol(>=) cutlist(50 75 90) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' above psacentile.dta")

					* below centile (e.g. PSA< 50th percentile)
					cindexmi, cohort(cohort==`cohortn') cuttype(psacentile) cutsymbol(<) cutlist(10 25 50) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' below psacentile.dta")

			}
		}
	}


	* repeating with a few more age cohorts.  These age groups are specified in the agegrp loop 
	* (where previosuly they were defined in the cohort variable)
	set more off
	foreach outcome in /*pca*/ dod {
		foreach covar in tpsa protectriskhg /*ftpsa ftpsariskhg*/ {
			foreach agegrp in 	"age<60 & tpsa<=4" "age>=60  & tpsa<=4" ///
								"age<60 & tpsa<=10" "age>=60  & tpsa<=10" ///
								"age<60 & tpsa<=25" "age>=60  & tpsa<=25" ///
								"!mi(cohort) & tpsa<=4" "!mi(cohort) & tpsa<=10" "!mi(cohort) & tpsa<=25" ///
								"!mi(cohort)" "age<60" "age>=60" "(age>=60 & age<70)" "(age>=55 & age<70)" "age<55" {
				*creating age group name that is file name friendly
				local agefilename="`agegrp'"
				local agefilename=subinstr("`agefilename'",">="," ge ",.)
				local agefilename=subinstr("`agefilename'","<="," le ",.)
				local agefilename=subinstr("`agefilename'","<" ," lt ",.)
				local agefilename=subinstr("`agefilename'",">" ," gt ",.)
				noi disp "`agefilename'"
				
				* above cuts (e.g. PSA>=2)
				cindexmi, cohort(`agegrp') cuttype(psalevel) cutsymbol(>=) cutlist(0 1 1.5 2 3) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' above psalevel.dta")

				* below cuts (e.g. PSA<1)
				cindexmi, cohort(`agegrp') cuttype(psalevel) cutsymbol(<) cutlist(1) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' below psalevel.dta")

				* above centile (e.g. PSA>= 50th percentile)
				cindexmi, cohort(`agegrp') cuttype(psacentile) cutsymbol(>=) cutlist(50 75 90) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' above psacentile.dta")

				* below centile (e.g. PSA< 50th percentile)
				cindexmi, cohort(`agegrp') cuttype(psacentile) cutsymbol(<) cutlist(10 25 50) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' below psacentile.dta")
			
			}
		}
	}

*********************************************************
**  Section 2. COMBINE RESULTS ACROSS IMPUTATION SETS  **
*********************************************************
*  Appending results from all imputation, subsets of patients, and PSA subsets
	!dir "Data\Results\Discrimination\*.dta" /b > "All dta files.txt"
	insheet using "All dta files.txt", clear delimiter(";")  
	drop if index(v1,"Bootstrap")

*appending datasets together to prepare results 
	* number of datasets t ocombine
	local N=`=_N'
	local i=0
	* saving a local var with each dataset's name
	while `i'<`N' {
		local ++i	
		local ds`i'=v1[`i']
	}

	* brining each dataset into memory and appending all together.
	local i=1
	use "Data\Results\Discrimination\\`ds`i''", clear
	while `i'<`N' {
		local ++i
		append using "Data\Results\Discrimination\\`ds`i''"
	}
	compress

	
**combining stats over 10 imputations, by taking the mean over the 10 imputed sets
	g errorn=!mi(error)
	collapse psacut psactlcut n casen cindex errorn, by(outcome covar cuttype cutsymbol cohort cutraw)
	tab errorn
	widetab cohort psacut errorn cindex if covar=="tpsa" & psacut==0 & outcome=="dod"


*not displaying all results
	* not displaying results for every cominbation of the data caluclated.
	drop if cuttype=="psacentile" & cutsymbol=="<" & inlist(cutraw,10,25)
	drop if outcome!="dod"
	keep if inlist(cohort,"!mi(cohort)")
	drop if covar=="ftpsa"

*calculating difference in cindex from tpsa alone
* extracting PSA alone, and merging it back into the master set.
	count
	foreach x in x {
		preserve 
		keep if covar=="tpsa"
		keep outcome cuttype cutsymbol cohort cutraw cindex
		rename cindex cindextpsa
		tempfile cindextpsa
		save `cindextpsa'
		restore
	}
	merge m:1 outcome cuttype cutsymbol cohort cutraw using `cindextpsa', nogen

*calculating differnece in cidnex from total and free PSA (same as was done for PSA)
	count
	foreach x in x {
		preserve 
		keep if covar=="ftpsariskhg"
		keep outcome cuttype cutsymbol cohort cutraw cindex
		rename cindex cindexftpsariskhg
		tempfile cindexftpsariskhg
		save `cindexftpsariskhg'
		restore
	}
	merge m:1 outcome cuttype cutsymbol cohort cutraw using `cindexftpsariskhg', nogen
	
******* merging in bootstrapped confidence intervals (calculated in separate analysis file)
	merge 1:1 outcome covar cuttype cutsymbol cohort cutraw using "Data\Results\Results - Discrimination Bootstrap Confidence Intervals.dta",
	*dropping bootstraps not being presented
	drop if _merge==2

********************************************************
**  Section 3. FORMAT DATA TO BE DISPLAYED IN TABLES  **
********************************************************
* this section is a bunch of tedious formatting code to get the results looking presentable.

*creating nicely formatted variables to display in tables
g psacutdisp="PSA"+cutsymbol+string(psacut,"%9.1f") /// (e.g. PSA >= 4)
g psacentldisp=string(psactlcut*100,"%9.0f")+"%" 
format %9.4f cindex
format %9.0f n casen

*more formatting of results to display in tables.
	* rounding c-index
	g cdisp=trim(string(cindex,"%9.3f"))
		* display NA if there were too many sets where c-index could not be calculated
		replace cdisp="NA" if errorn>0.2
	* rounding differnece in cindex
	g cdispdelta=trim(string(cindex - cindextpsa,"%9.3f"))
		* display NA if there were too many sets where c-index could not be calculated
		replace cdispdelta="NA" if errorn>0.2 
		* display -- if PSA (which is the base model, so the improvemetn is of course 0)
		replace cdispdelta="--" if covar=="tpsa"
	* repeat but for imprevoemtn over free+total PSA+age
	g cdispdeltaftpsa=trim(string(cindex - cindexftpsariskhg,"%9.3f"))
		replace cdispdeltaftpsa="NA" if errorn>0.2 
		replace cdispdeltaftpsa="--" if inlist(covar,"tpsa","ftpsariskhg")

	* rounding/formatting the bootstrap-estimated 95% CI.
	g cdispdeltaci=trim(string(cdispdeltalb,"%9.3f"))+", "+trim(string(cdispdeltaub,"%9.3f")) if covar!="tpsa"
		replace cdispdeltaci="--" if covar=="tpsa"
	g cdispdeltaftpsaci=trim(string(cdispdeltaftpsalb,"%9.3f"))+", "+trim(string(cdispdeltaftpsaub,"%9.3f")) if !inlist(covar,"tpsa","ftpsariskhg")
		replace cdispdeltaftpsaci="--" if inlist(covar,"tpsa","ftpsariskhg")
	
*ordering results for the table
	g outcomesort=1 if outcome=="cancer"
	replace outcomesort=2 if outcome=="mets"
	replace outcomesort=3 if outcome=="dod"

	g covarsort=1 if covar=="tpsa"
	replace covarsort=2 if covar=="ftpsa"
	replace covarsort=3 if covar=="ftpsariskhg"
	replace covarsort=4 if covar=="protectriskhg"


*there are just a handful of cases, where the estimated centile is a little different for the varying covars
*this very rarely means that tPSA and free to total PSA are assessed on very slightly different cohorts (off by one or two patients)
	bysort outcome cohort  cuttype cutsymbol cutraw: egen psacut0=mean(psacut)
	g diff = psacut-psacut0
	replace psacut=psacut0

* creating age label for figures and tables
	gsort outcome cohort cutsymbol psacut covarsort
	by outcome cohort: g cohortdisp="All Ages" if _n==1 & cohort=="!mi(cohort)"
	by outcome cohort: replace cohortdisp="45-49" if _n==1 & cohort=="cohort==45"
	by outcome cohort: replace cohortdisp="50-54" if _n==1 & cohort=="cohort==50"
	by outcome cohort: replace cohortdisp="55-59" if _n==1 & cohort=="cohort==55"
	by outcome cohort: replace cohortdisp="60-64" if _n==1 & cohort=="cohort==60"
	by outcome cohort: replace cohortdisp="65-69" if _n==1 & cohort=="cohort==65"
	by outcome cohort: replace cohortdisp="70-74" if _n==1 & cohort=="cohort==70"
	by outcome cohort: replace cohortdisp="60-69" if _n==1 & cohort=="(age>=60 & age<70)"
	by outcome cohort: replace cohortdisp="<60" if _n==1 & cohort=="age<60"
	by outcome cohort: replace cohortdisp="60+" if _n==1 & cohort=="age>=60"
	by outcome cohort: replace cohortdisp=cohort if _n==1 & mi(cohortdisp)

*only dipslaying first row within age and outcome
	replace psacentldisp="" if covar!="tpsa"
	replace psacutdisp="" if covar!="tpsa"
	replace n=. if covar!="tpsa"
	replace casen=. if covar!="tpsa"

*displaying covariates
	g covardisp="tPSA" if covar=="tpsa"
	replace covardisp="Free-to-Total PSA Ratio" if covar=="ftpsa"
	replace covardisp="Age+tPSA+fPSA" if covar=="ftpsariskhg"
	replace covardisp="Age+tPSA+fPSA+iPSA+hK2" if covar=="protectriskhg"

* printing table showing imrpovement over tpsa
	foreach outcome in dod {
		preserve
		disp _newline "`outcome'"
		listtex cohortdisp psacutdisp psacentldisp n casen covardisp cdisp cdispdelta cdispdeltaci cdispdeltaftpsa cdispdeltaftpsaci if outcome=="`outcome'", type
		restore
	}

