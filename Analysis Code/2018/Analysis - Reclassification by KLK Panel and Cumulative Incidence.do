version 12
sysdir set SITE "O:\Outcomes\Andrew\_UROLOGY Biostatistics\MACROS etc\ado files\"
/*
PROGRAM: 11 Analysis - Reclassification by KLK Panel and Cumulative Incidence.do
DESCRIPTION:
	The code below calcuatles the cumultive incidence within PSA and PSA+KLK subgroups.
	The key here to to idenfity how many men with elevated PSA can be recalssfied as low risk if they have a low KLK score.
*/

*number of imputations
global m=10
cd "O:\Outcomes\Andrew\Analytic Projects\1 Active\Lilja MDC KLK predicts PCa Outcomes"

* tihs program calculates cumulative incidence in the presenece of the imputed data
capture program drop cuminc
program cuminc
	syntax, 	cohort(string asis) /// speicifes the age subgroup
				subgroup(string) /// specifies the PSa subgroup
				outcome(string) /// specifies the otucome variables name
				[saving(string asis) by(passthru)]
	
		* saving a copy of the data in it's current form
		preserve
		
		* a numlist that lists the times where teh cum inc will be calculated.
		if trim("`at'")=="" local at="0(0.1)20"
		
		/*cycling through each imputed value*/
		forvalues m=1(1)${m} {
			qui {
					restore, preserve
					* keeping observations within age group
					keep if `cohort'
					* storing the mth imputation into memory
					mi extract `m', clear
					
					
					/*keeping PSA sub-group of interest*/
					keep if `subgroup'
					
					**  getting descriptive stats for subgroup (only needed for iterations with by var)
						if `"`by'"'!="" {
							* total N
							count
								local ntot=`r(N)'
							* n within each group	
							levelsof `=substr("`by'",4,length("`by'")-4)'
								foreach l in `r(levels)' {
									count if `=substr("`by'",4,length("`by'")-4)'==`l'
									local `=substr("`by'",4,length("`by'")-4)'`l'=`r(N)'
								}
						}

					* prepping the data for survival analysis
					stset tt`outcome', f(`outcome')

					/*calculating survival/failure stats*/
					tempfile dataset`m'
					sts list, at(`at') saving(`dataset`m'', replace) failure  `by'
					
					**  updating datasets to look like old structure so I don't have to update any subsequent code
					use  `dataset`m'', clear
					drop fail begin

					g m=`m'
					
					**  calculating the se from the results
					g logloglb=log(-log(lb))
					g loglogub=log(-log(ub))
					g loglogfail=log(-log(failure))
					g loglogse=abs((loglogub-loglogfail)/invnormal(0.975))
					
					*saving descriptive stats
						if `"`by'"'!="" {
							g ntot=`ntot'
							g nklksubgroup=.
							
							levelsof `=substr("`by'",4,length("`by'")-4)'
								foreach l in `r(levels)' {
									replace nklksubgroup=``=substr("`by'",4,length("`by'")-4)'`l'' if `=substr("`by'",4,length("`by'")-4)'==`l'
								}
						}

					* saving resulting from mth imputation
					save `dataset`m'', replace
				}
		}
		
		qui{
		* append all of the datasets together
				use `dataset1', clear
				forvalues b=2(1)${m} {
					append using `dataset`b''
				}

				sort time m
														
				g outcome="`outcome'"
				g cohort=`"`cohort'"'
				g subgroup="`subgroup'"
		}
		if `"`saving'"'!="" save `saving'
		
		restore
end

**  testing that th functio nworks properly with an example
use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear
cuminc, cohort(!mi(cohort)) subgroup(tpsa>=0) outcome(dod) by(pca)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	*  in this section of code (sandwiched by the xxxxxx) calculates the cum inc within many subgroups of PSA and PSA+KLK.
	*  the results are later saved to be tablulated for tables and figures.

		use "Data\Master - Lilja MDC KLK predicts distant mets MDC Population-based Cohort with Imputed KLK.dta", clear
		foreach v in protectriskhg ftpsariskhg {
			mi passive: g `v'025=`v'>=025/1000
			mi passive: g `v'050=`v'>=050/1000
			mi passive: g `v'060=`v'>=060/1000
			mi passive: g `v'075=`v'>=075/1000
			mi passive: g `v'100=`v'>=100/1000
		}
		
		local count=0
		* looping over varying age cohorts for which the reclassification stats will be calculated.
		foreach cohort in 	"age<60" /*"age>=60"*/ ///
							/*"age<60 & !(dod==1 & ttdod < 3)" "age>=60 & !(dod==1 & ttdod < 3)"*/ ///
							/*"age<60 & !(dod==1 & ttdod < 5)" "age>=60 & !(dod==1 & ttdod < 5)" */ ///
							/*"age>=60 & age<70" "cohort==50" "cohort==55" "cohort==60" "cohort==65" "cohort==70" "!mi(cohort)"*/ {
			* looping over various PSA subgroups here
			foreach psagroup in /*"tpsa<=0.5" "tpsa<=0.75" "tpsa<=1" "tpsa<=1.25" "tpsa<=1.5" "tpsa<=1.75" "tpsa<=1.8" "tpsa<=2"*/ ///
								/*"tpsa>=0" "tpsa>=1.5" "tpsa>=2" "tpsa>=3" "inrange(tpsa,2,3)"*/ ///
								/*"inrange(tpsa,0,25)" "inrange(tpsa,1.5,25)" "inrange(tpsa,2,25)" "inrange(tpsa,3,25)"*/ ///
								"inrange(tpsa,1.5,10)" /*"inrange(tpsa,2,10)" "inrange(tpsa,3,10)"*/ ///
								/*"inrange(tpsa,1.5,4)" "inrange(tpsa,2,4)"*/ {
				* calcualting the reclassification statistics here
				local ++count
				cuminc, cohort(`cohort') subgroup(`psagroup') outcome(dod) ///
					saving("Data\Results\Cum Inc Reclassification\Results - Cum Inc Reclassification `count'.dta", replace)
				
				*no need to split risk by KLK for lower PSA ranges
				if inlist("`psagroup'", "tpsa<=0.5", "tpsa<=0.75", "tpsa<=1", "tpsa<=1.25", "tpsa<=1.5", "tpsa<=1.75", "tpsa<=1.8", "tpsa<=2") continue
				foreach klkgroup in "025" "050" "060" "075" "100" {
					
					local ++count
					*calculating reclassfication for full KLK model
					cuminc, cohort(`cohort') subgroup(`psagroup') outcome(dod) by(protectriskhg`klkgroup') ///
						saving("Data\Results\Cum Inc Reclassification\Results - Cum Inc Reclassification `count'.dta", replace)

					local ++count
					*calculating reclassfication for age tpsa and fpsa model
					cuminc, cohort(`cohort') subgroup(`psagroup') outcome(dod) by(ftpsariskhg`klkgroup') ///
						saving("Data\Results\Cum Inc Reclassification\Results - Cum Inc Reclassification `count'.dta", replace)
				}
			}	
		}

		**  appending all the reuslts togeteher in one dataset
		use "Data\Results\Cum Inc Reclassification\Results - Cum Inc Reclassification 1.dta", clear
		foreach i of numlist 2/`count' {
			append using "Data\Results\Cum Inc Reclassification\Results - Cum Inc Reclassification `i'.dta"
		}
		save "Data\Results\Results - Cum Inc Reclassification", replace
		
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

**  this section takes all the results, calculates the 95%CI, and makes cosmetic updates the data to display in tables in figures.

use  "Data\Results\Results - Cum Inc Reclassification", clear
* no need to plot results for subgroups of low PSA
	drop if index( subgroup, "tpsa<=")

*creating and indicator whether data is with psa subgroup only, or split by klk
	egen klksplit=rownonmiss(protectriskhg??? ftpsariskhg???)
		widetab klksplit protectriskhg??? ftpsariskhg???

*calculating central estimate of failure probability
	sort cohort outcome subgroup protectriskhg??? ftpsariskhg??? time m
	by cohort outcome subgroup protectriskhg??? ftpsariskhg??? time: egen failure_central=mean(failure)

*calculating the corrected stderr (adding within and between varariance)
	g loglogvar=loglogse^2
	by cohort outcome subgroup protectriskhg??? ftpsariskhg??? time: egen loglogvar_mean=mean(loglogvar)
	g betweenvar_temp=((log(-log(failure))-log(-log(failure_central)))^2)/(${m}-1)
	by cohort outcome subgroup protectriskhg??? ftpsariskhg??? time: egen betweenvar=sum(betweenvar_temp)
	g loglogtotalvar=loglogvar_mean+(1+1/${m})*betweenvar
	g loglogtotalse=sqrt(loglogtotalvar)


* combining descriptive stats about cohort sizes
	by cohort outcome subgroup protectriskhg??? ftpsariskhg???: egen nklksubgroup_mean=mean(nklksubgroup)
	by cohort outcome subgroup protectriskhg??? ftpsariskhg???: egen ntot_mean=mean(ntot)
	
	* keeping one line per combinatio of age, PSA subgroup, etc. (continuing with calculations for MI data's SEs)
	keep if m==1
	g failure_ub=failure_central^(exp(invnormal(0.975)*loglogtotalse))*100
	g failure_lb=failure_central^(exp(-invnormal(0.975)*loglogtotalse))*100
	replace failure_central=failure_central*100	
	
	* formatting the results to display in a table.
	g fail_disp=trim(string(failure_central,"%9.2f"))+ " ("+trim(string(failure_ub,"%9.2f"))+" - "+trim(string(failure_lb,"%9.2f"))+")"
	replace fail_disp="0 (NA)" if failure_central==0
	replace fail_disp="NA (NA)" if failure_central==.

	
*** dropping variables from one line per imputation
	drop failure std_err lb ub m logloglb loglogub loglogfail loglogse ntot nklksubgroup loglogvar loglogvar_mean betweenvar_temp betweenvar loglogtotalvar loglogtotalse

**  getting proporion of men in klk subgroup
	capture drop klksubgroupprop
	g klksubgroupprop=string(nklksubgroup_mean/ntot_mean*100,"%9.0f")+"%" if klksplit==1
	
	g subgroupriskabove=. 
	g subgrouprisk = "" 
		foreach v of varlist protectriskhg??? ftpsariskhg??? {
			replace subgrouprisk=substr("`v'",1,length("`v'")-3) if !mi(`v')
			replace subgroupriskabove=`v' if !mi(`v')
		}
	
	
* getting all levels of ochort to loop over for displaying statistics
levelsof cohort, local(lcohort)
	*local lcohort=`"age>=60 age<60"'
levelsof subgroup, local(lsubgroup)
	*local lsubgroup=`"inrange(tpsa,1.5,4) inrange(tpsa,2,4)"'

	
	
***************************
* looping over age cohorts, and creating figure of relcassificatoin results.
* also saving out estiamtes to be displayed in a table.
***************************	
g precalssified=""
g cohortdisp=""
foreach cohort in `lcohort' {
	*creating nice display versions of the cohort for display
	local ylabel=""
	if "`cohort'"=="!mi(cohort)" {
		local cohortdisp="All Ages"
		*local ylabel="0(5)20"
	}
	else if "`cohort'"=="cohort==50" {
		 local cohortdisp="Ages 50-54"
		*local ylabel="0(5)20"
	}
	else if "`cohort'"=="cohort==55" {
		 local cohortdisp="Ages 55-59"
		*local ylabel="0(5)20"
	}
	else if "`cohort'"=="cohort==60" {
		 local cohortdisp="Ages 60-64"
		*local ylabel="0(5)20"
	}
	else if "`cohort'"=="cohort==65" {
		 local cohortdisp="Ages 65-69"
		*local ylabel="0(5)20"
	}
	else if "`cohort'"=="cohort==70" {
		 local cohortdisp="Ages 70-74"
		*local ylabel="0(5)20"
	}
	else if "`cohort'"=="age<60" {
		 local cohortdisp="Ages less than 60"
		 local ylabel="ylabel(0(5)15)"
	}	
	else if "`cohort'"=="age>=60" {
		 local cohortdisp="Ages 60 or greater"
		 local ylabel="ylabel(0(5)25)"
	}	
	else if "`cohort'"=="age>=60 & age<70" {
		 local cohortdisp="Ages 60-69"
		 local ylabel="ylabel(0(5)20)"
	}	
	else if "`cohort'"=="age<60 & !(dod==1 & ttdod < 3)" {
		 local cohortdisp="Ages less than 60 (excluding PCa deaths within 3 years)"
		 local ylabel="ylabel(0(5)15)"
	}	
	else if "`cohort'"=="age>=60 & !(dod==1 & ttdod < 3)" {
		 local cohortdisp="Ages 60 or greater (excluding PCa deaths within 3 years)"
		 local ylabel="ylabel(0(5)25)"
	}
	else if "`cohort'"=="age<60 & !(dod==1 & ttdod < 5)" {
		 local cohortdisp="Ages less than 60 (excluding PCa deaths within 5 years)"
		 local ylabel="ylabel(0(5)15)"
	}	
	else if "`cohort'"=="age>=60 & !(dod==1 & ttdod < 5)" {
		 local cohortdisp="Ages 60 or greater (excluding PCa deaths within 5 years)"
		 local ylabel="ylabel(0(5)25)"
	}
	else  {
		local cohortdisp="`cohort'"
		*local ylabel="0(5)20"
	}
	replace cohortdisp="`cohortdisp'" if cohort=="`cohort'"
	
	foreach subgroup in `lsubgroup' {
		*creating nice display versions of the subgroup for display
		if "`subgroup'"=="inrange(tpsa,0,25)" local subgroupdisp="0{&le}PSA{&le}25"
		else if "`subgroup'"=="inrange(tpsa,2,25)" local subgroupdisp="2{&le}PSA{&le}25"
		else if "`subgroup'"=="inrange(tpsa,2,10)" local subgroupdisp="2{&le}PSA{&le}10"
		else if "`subgroup'"=="inrange(tpsa,2,4)" local subgroupdisp="2{&le}PSA{&le}4"
		else if "`subgroup'"=="inrange(tpsa,1.5,25)" local subgroupdisp="1.5{&le}PSA{&le}25"
		else if "`subgroup'"=="inrange(tpsa,1.5,10)" local subgroupdisp="1.5{&le}PSA{&le}10"
		else if "`subgroup'"=="inrange(tpsa,1.5,4)" local subgroupdisp="1.5{&le}PSA{&le}4"
		else if "`subgroup'"=="inrange(tpsa,2,3)" local subgroupdisp="2{&le}PSA{&le}3"
		else if "`subgroup'"=="inrange(tpsa,3,25)" local subgroupdisp="3{&le}PSA{&le}25"
		else if "`subgroup'"=="inrange(tpsa,3,10)" local subgroupdisp="3{&le}PSA{&le}10"
		else if "`subgroup'"=="tpsa>=0" local subgroupdisp="PSA{&ge}0"
		else if "`subgroup'"=="tpsa>=1.5" local subgroupdisp="PSA{&ge}1.5"
		else if "`subgroup'"=="tpsa>=2" local subgroupdisp="PSA{&ge}2"
		else if "`subgroup'"=="tpsa>=3" local subgroupdisp="PSA{&ge}3"
		else local subgroupdisp="`subgroup'"
		
		**  cycling over risk score type, and creating figures for each.
		foreach risk in protectriskhg ftpsariskhg {
				* creating label for the risk model
				if "`risk'"=="protectriskhg" local risktitle="ProtecT" 
				else if "`risk'"=="ftpsariskhg" local risktitle="Age, tPSA, fPSA" 
		
				* looping over every risk cutpoint.
				local klkgraphlist=""
				foreach klkcut in 025 050 060 075 100 {
					local klkdisp=string(`klkcut'/10,"%9.1f")+"%"
					
					*getting proportions in klk subgroups (preserving data as to not 
					qui foreach x in x {
						preserve
						*set trace on
						keep if cohort=="`cohort'" & subgroup=="`subgroup'" 
						sum failure_central if time<=20 & klksplit==0
							local psa_mid=`r(max)'
							
						keep if !mi(`risk'`klkcut')

						sort `risk'`klkcut'
							local klk0=klksubgroupprop[1]
							local klk1=klksubgroupprop[`=_N']
							
						sum failure_central if time<=20 & `risk'`klkcut'==1
							if `r(N)'>0 local psa_above=`r(max)'
							else local psa_above=0
						sum failure_central if time<=20 & `risk'`klkcut'==0
							if `r(N)'>0 local psa_below=`r(max)'
							else local psa_below=0
						restore
					}
					replace precalssified="`klk0'" if cohort=="`cohort'" & subgroup=="`subgroup'" & `risk'`klkcut'==0
					replace precalssified="`klk1'" if cohort=="`cohort'" & subgroup=="`subgroup'" & `risk'`klkcut'==1
					
					**  finally plotting the results!~
					twoway	(line failure_central time if time<=20 & cohort=="`cohort'" & subgroup=="`subgroup'" & klksplit==0, sort lpattern(solid) lcolor(black) connect(stairstep) yaxis(1 2)) ///
							(line failure_central time if time<=20 & cohort=="`cohort'" & subgroup=="`subgroup'" & `risk'`klkcut'==1, sort lpattern(solid) lcolor(red) connect(stairstep)) ///
							(line failure_central time if time<=20 & cohort=="`cohort'" & subgroup=="`subgroup'" & `risk'`klkcut'==0, sort lpattern(solid) lcolor(blue) connect(stairstep)), ///
							scheme(s1mono) legend(off) xtitle("Years since blood draw", margin(medium)) xlabel(0(5)20) ytitle(" ", axis(1)) ytitle(" ", axis(2))  ///
							note("`risktitle' cut at `klkdisp'") `ylabel' ///
							ylabel(`psa_mid' `"All men `subgroupdisp'"' `psa_above' `"`klk1' of men `subgroupdisp'"' `psa_below' `"`klk0' of men `subgroupdisp'"', angle(0) labsize(vsmall) axis(2)) ///
							saving("Figures\GPH\Figure - Cum Inc Reclassification by KLK `cohortdisp', `subgroupdisp', KLK threshold `klkcut'.gph", replace)

					* saving out figure to folder
					if inlist("`klkcut'","060","075","100") {
						local klkgraphlist=`"`klkgraphlist' "Figures\GPH\Figure - Cum Inc Reclassification by `risk' `cohortdisp', `subgroupdisp', KLK threshold `klkcut'.gph""'
						graph export "Figures\Figure - Cum Inc Reclassification by `risk' `cohortdisp', `subgroupdisp', `risk' threshold `klkcut'.tif", replace
					}
					
				}

				/*
				graph combine `klkgraphlist', title(`"Kaplan-Meier probability of PCa Death: `cohortdisp' `subgroupdisp'"', size(2.5)  position(11)) ///
						iscale(*0.9) cols(3) scheme(s1mono) imargin(zero)
							
					graph export "Figures\Figure - Cum Inc Reclassification by KLK `cohortdisp' `subgroupdisp'.tif", as(tif) replace width(1600)
				*/
		}
	}
}
* saving all estimates to be reported in tables
save "Data\Results\Results - Cum Inc Reclassification (imputations combined).dta", replace



graph combine 	"Figures\GPH\Figure - Cum Inc Reclassification by protectriskhg Ages less than 60, PSA{&ge}2, KLK threshold 060.gph" ///
				"Figures\GPH\Figure - Cum Inc Reclassification by protectriskhg Ages less than 60, PSA{&ge}2, KLK threshold 075.gph" ///
				"Figures\GPH\Figure - Cum Inc Reclassification by protectriskhg Ages less than 60, PSA{&ge}2, KLK threshold 100.gph" ///
				, title(`"Kaplan-Meier probability of PCa Death: Ages less than 60, PSA{&ge}2"', size(2.5)  position(11)) ///
				iscale(*0.8) cols(3) scheme(s1mono) imargin(zero)
			graph export "Figures\test.tif", replace width(3100) height(700)


**********  DISPLAING RESULTS IN TABLE  **************
	use "Data\Results\Results - Cum Inc Reclassification (imputations combined).dta", clear
	keep if inlist(time,5,10,15,20)
	drop failure_central failure_ub failure_lb
	
	egen riskid = group( protectriskhg??? ftpsariskhg???), missing
	
	reshape wide fail_disp, i(outcome cohort subgroup riskid) j(time)
	
	* creating a sort order for the table
	g tpsasubgroupsort=1 if subgroup=="tpsa>=0"
		replace tpsasubgroupsort=2 if subgroup=="inrange(tpsa,0,25)"
		replace tpsasubgroupsort=2.5 if subgroup=="tpsa>=1.5"
		replace tpsasubgroupsort=3 if subgroup=="tpsa>=2"
		replace tpsasubgroupsort=4 if subgroup=="inrange(tpsa,2,25)"
		replace tpsasubgroupsort=5 if subgroup=="tpsa>=3"
		replace tpsasubgroupsort=6 if subgroup=="inrange(tpsa,3,25)"
		replace tpsasubgroupsort=7 if subgroup=="inrange(tpsa,2,3)"

	sort cohort tpsasubgroupsort klksplit protectriskhg??? ftpsariskhg???

	* creaitng nice PSA labels for table
	g subgroupdisp=subinstr(subgroup,"tpsa","PSA",1) if inlist(subgroup,"tpsa>=0","tpsa>=1.5","tpsa>=2","tpsa>=3") & klksplit==0
	replace subgroupdisp="0<=PSA<=25" if subgroup=="inrange(tpsa,0,25)" & klksplit==0
	replace subgroupdisp="2<=PSA<=25" if subgroup=="inrange(tpsa,2,25)" & klksplit==0
	replace subgroupdisp="3<=PSA<=25" if subgroup=="inrange(tpsa,3,25)" & klksplit==0
	replace subgroupdisp="2<=PSA<=3" if subgroup=="inrange(tpsa,2,3)" & klksplit==0
	
	* again, making nice lables for table
	foreach riskvar in protectriskhg ftpsariskhg {
		if "`riskvar'"=="protectriskhg" local label="KLK Risk"
		else if "`riskvar'"=="ftpsariskhg" local label="Age, tPSA, fPSA Risk"
		foreach klk in 025 050 060 075 100 {
			replace subgroupdisp="    `label'<"+string(`klk'/10,"%9.1f")+"%" if `riskvar'`klk'==0
			replace subgroupdisp="                  >="+string(`klk'/10,"%9.1f")+"%" if `riskvar'`klk'==1
		}
	}
	
	
	**  PRINT RESULTS TO TABLE
	*levelsof cohort
	*foreach c in `r(levels)' {
	foreach c in /*"!mi(cohort)"*/ "age<60" "age>=60" /*"age>=60 & age<70"*/ {
		preserve 
		qui keep if cohort=="`c'"
		qui keep if inlist(subgroup,"tpsa>=0","tpsa>=1.5","tpsa>=2","tpsa>=3")
		*qui keep if inlist(subgroup,"tpsa>=0","tpsa>=1.5","tpsa>=2")
		*qui keep if klksplit==0 | !mi(protectriskhg075) | !mi(ftpsariskhg075)
		local cohortdisp=cohortdisp
		noi disp _newline "`cohortdisp'"
		listtex subgroupdisp precalssified fail_disp5 fail_disp10 fail_disp15 fail_disp20 if cohort=="`c'" , ///
			type headlines("&Proportion Reclassified&5 Year Risk&10 Year Risk&15 Year Risk&20 Year Risk")
		restore
	}















