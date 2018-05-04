cd "O:\Outcomes\Andrew\Analytic Projects\1 Active\Lilja MDC KLK predicts PCa Outcomes"
use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear

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
				if "`cuttype'"=="psalevel" {
					local psacut=`cut'
					*getting the centile
					count if tpsa<`psacut'
						local nbelow=`r(N)'
					count if tpsa<=`psacut'
						local nabove=`r(N)'
					local ctlcut=((`nabove' + `nbelow')/2)/`=_N'	
				}
				else if "`cuttype'"=="psacentile" {
					local ctlcut=`cut'/100
					centile tpsa, c(`cut')
						local psacut=`r(c_1)'
				}
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
					stcox `covar' 
					estat concord
					post `memhold' ("`outcome'") ("`covar'") ("`cuttype'") ("`cutsymbol'") ///
									 ("`cohort'") (`cut') (`psacut') (`ctlcut') (`m') (`N') (`caseN') (`r(C)') ("")	
				}
			*if not enough data, posting a missing
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

						
*  C-INDEX FOR TPSA and KLK Panel
	x
	use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear
	levelsof cohort, local(cohortlevels)
	set more off
	foreach outcome in /*pca*/ dod {
		foreach cohortn in `cohortlevels' {
			foreach covar in tpsa protectriskhg ftpsa ftpsariskhg {
					* above cuts
					cindexmi, cohort(cohort==`cohortn') cuttype(psalevel) cutsymbol(>=) cutlist(0 1 1.5 2 3) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' above psalevel.dta")

					* below cuts
					cindexmi, cohort(cohort==`cohortn') cuttype(psalevel) cutsymbol(<) cutlist(1) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' below psalevel.dta")

					* above centile
					cindexmi, cohort(cohort==`cohortn') cuttype(psacentile) cutsymbol(>=) cutlist(50 75 90) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' above psacentile.dta")

					* below centile
					cindexmi, cohort(cohort==`cohortn') cuttype(psacentile) cutsymbol(<) cutlist(10 25 50) outcome(`outcome') covar(`covar') ///
						saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' Cohort `cohortn' below psacentile.dta")

			}
		}
	}


	* repeating with a few more age cohorts
	set more off
	foreach outcome in /*pca*/ dod {
		foreach covar in tpsa protectriskhg /*ftpsa ftpsariskhg*/ {
			foreach agegrp in 	"age<60 & tpsa<=4" "age>=60  & tpsa<=4" ///
								"age<60 & tpsa<=10" "age>=60  & tpsa<=10" ///
								"age<60 & tpsa<=25" "age>=60  & tpsa<=25" ///
								"!mi(cohort) & tpsa<=4" "!mi(cohort) & tpsa<=10" "!mi(cohort) & tpsa<=25" ///
								"!mi(cohort)" "age<60" "age>=60" "(age>=60 & age<70)" "(age>=55 & age<70)" "age<55" {
				*creating age group that is file name friendly
				local agefilename="`agegrp'"
				local agefilename=subinstr("`agefilename'",">="," ge ",.)
				local agefilename=subinstr("`agefilename'","<="," le ",.)
				local agefilename=subinstr("`agefilename'","<" ," lt ",.)
				local agefilename=subinstr("`agefilename'",">" ," gt ",.)
				noi disp "`agefilename'"
				

				* above cuts
				cindexmi, cohort(`agegrp') cuttype(psalevel) cutsymbol(>=) cutlist(0 1 1.5 2 3) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' above psalevel.dta")

				* below cuts
				cindexmi, cohort(`agegrp') cuttype(psalevel) cutsymbol(<) cutlist(1) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' below psalevel.dta")

				* above centile
				cindexmi, cohort(`agegrp') cuttype(psacentile) cutsymbol(>=) cutlist(50 75 90) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' above psacentile.dta")

				* below centile
				cindexmi, cohort(`agegrp') cuttype(psacentile) cutsymbol(<) cutlist(10 25 50) outcome(`outcome') covar(`covar') ///
					saving("Data\Results\Discrimination\Results - Discrimination of `covar' `outcome' `agefilename' below psacentile.dta")
			
			}
		}
	}
	
*  Appending results from all imputation, subsets of patients, and PSA subsets
	!dir "Data\Results\Discrimination\*.dta" /b > "All dta files.txt"
	insheet using "All dta files.txt", clear delimiter(";")  
	drop if index(v1,"Bootstrap")

*appending datasets, and extracting data from the file name
	local N=`=_N'
	local i=0
	while `i'<`N' {
		local ++i	
		local ds`i'=v1[`i']
	}

	local i=1
	use "Data\Results\Discrimination\\`ds`i''", clear
	while `i'<`N' {
		local ++i
		append using "Data\Results\Discrimination\\`ds`i''"
	}
	compress

	
**combining stats over 10 imputations
	g errorn=!mi(error)
	collapse psacut psactlcut n casen cindex errorn, by(outcome covar cuttype cutsymbol cohort cutraw)
	tab errorn
	widetab cohort psacut errorn cindex if covar=="tpsa" & psacut==0 & outcome=="dod"



*not displaying all results
	drop if cuttype=="psacentile" & cutsymbol=="<" & inlist(cutraw,10,25)
	drop if outcome!="dod"
	keep if inlist(cohort,"!mi(cohort)")
	*keep if inlist(cohort,"!mi(cohort)","(age>=60 & age<70)","age<60","age>=60")
	*keep if inlist(cohort,"!mi(cohort) & tpsa<=10","!mi(cohort) & tpsa<=25","!mi(cohort) & tpsa<=4")
	drop if covar=="ftpsa"

*calculatind differnece in cindex from tpsa alone
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

*calculatind differnece in cidnex from total and free PSA
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
	
******* merging in bootstrapped ocnfidence intervals
	merge 1:1 outcome covar cuttype cutsymbol cohort cutraw using "Data\Results\Results - Discrimination Bootstrap Confidence Intervals.dta",
	*dropping bootstraps not being presented
	drop if _merge==2

	
*creating nicely formatted variables to display in tables
g psacutdisp="PSA"+cutsymbol+string(psacut,"%9.1f")
g psacentldisp=string(psactlcut*100,"%9.0f")+"%"
format %9.4f cindex
format %9.0f n casen

*more formatting of results to display in tables.
	g cdisp=trim(string(cindex,"%9.3f"))
	replace cdisp="NA" if errorn>0.2
	g cdispdelta=trim(string(cindex - cindextpsa,"%9.3f"))
	replace cdispdelta="NA" if errorn>0.2 
	replace cdispdelta="--" if covar=="tpsa"
	g cdispdeltaftpsa=trim(string(cindex - cindexftpsariskhg,"%9.3f"))
	replace cdispdeltaftpsa="NA" if errorn>0.2 
	replace cdispdeltaftpsa="--" if inlist(covar,"tpsa","ftpsariskhg")

	g cdispdeltaci=trim(string(cdispdeltalb,"%9.3f"))+", "+trim(string(cdispdeltaub,"%9.3f")) if covar!="tpsa"
		replace cdispdeltaci="--" if covar=="tpsa"
	g cdispdeltaftpsaci=trim(string(cdispdeltaftpsalb,"%9.3f"))+", "+trim(string(cdispdeltaftpsaub,"%9.3f")) if !inlist(covar,"tpsa","ftpsariskhg")
		replace cdispdeltaftpsaci="--" if inlist(covar,"tpsa","ftpsariskhg")
	
*ordering results
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

