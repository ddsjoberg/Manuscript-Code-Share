/*
PROGRAM: Analysis - Lee Robot vs Open RC Recurrence Comparison 10-016
PROGRAMMER: Daniel
DATE: 10/24/2016
PURPOSE: 
The following code builds upon a previous mansucript where post-op complications wer ecompared for 
and robotic radical cystectomy.  In this project we are comparing oncologic outcomes, such as 
recurrence and death from cancer after surgery.
*/

cd "O:\Outcomes\Andrew\Prospective Protocols\Laudone Robotic v Open Cystectomy Complication Rates 10-016\Lee Robot vs Open RC Recurrence Comparison 10-016"
use "Data\Master - Lee Robot vs Open RC Recurrence Comparison 10-016.dta", clear

tab surgtype recur, m
* getting followup
sum ttdead if dead==0, d
tab1 recur dod dead 

***********************************
**  distant recurrence analyses  **
***********************************
	* setting data for time to distant recurance outcomes
	stset ttrecurdist recurdist
	* log-rank test for differnce between open and robot
	sts test surgtype

*******************************
**  any recurrence analyses  **
*******************************
	* setting the data for time to recurrence outcomes
	stset ttrecur recur

	* calculating HR for type of surgery and printing formatted results
	stcox surgtype
	printmodel, text 

	* does timing of surgery predict outcomes among robot (ie is there a learning curve)
	stcox surgdate if surgtype==0
	stcox surgdate if surgtype==1

	* testing for a difference in recurrence by surgery type (log-rank)
	sts test surgtype
	
	* kaplan meier figures for recurrence (first survival, then failure)
		sts graph, fail by(surgtype) xlabel(0/6) tmax(6) scheme(s1mono) ///
				plot1opts(lpattern(dash) lcolor(black)) ///
				plot2opts(lpattern(solid) lcolor(black)) ///
				risktable( ,order(1 "Robotic" 2 "Open")) ///
				legend(label(1 "Robotic") label(2 "Open")) ///
				ylabel(,angle(0)) title(" ", position(11)) ///
				xtitle("Years since surgery", margin(medium)) ///
				ytitle("Recurrence Probability", margin(medium))  ///
				saving("Figures\Kaplan Meier - Recurrence Probability by Surgery Type.gph", replace)
		graph export "Figures\Kaplan Meier - Recurrence Probability by Surgery Type.tif", replace		

		sts graph, by(surgtype) xlabel(0/6) tmax(6) scheme(s1mono) ///
				plot1opts(lpattern(dash) lcolor(black)) ///
				plot2opts(lpattern(solid) lcolor(black)) ///
				risktable( ,order(1 "Robotic" 2 "Open")) ///
				legend(label(1 "Robotic") label(2 "Open")) ///
				ylabel(,angle(0)) title(" ", position(11)) ///
				xtitle("Years since surgery", margin(medium)) ///
				ytitle("Recurrence-free Probability", margin(medium)) ///
				saving("Figures\Kaplan Meier - Recurrence-free Probability by Surgery Type.gph", replace) 
		graph export "Figures\Kaplan Meier - Recurrence-free Probability by Surgery Type.tif", replace		

	* printing KM failure probabilites	
	stkmdiff, by(surgtype) time(3 5) fail table 

****************************************
**  death from cancer (dod) analyses  **
****************************************
	stset ttdod dod

	* testing for differnce in death from disease (dod)
	sts test surgtype
	
	* kaplan meier figures for dod (first survival, then failure)
		sts graph, fail by(surgtype) xlabel(0/6) tmax(6) scheme(s1mono) ///
				plot1opts(lpattern(dash) lcolor(black)) ///
				plot2opts(lpattern(solid) lcolor(black)) ///
				risktable( ,order(1 "Robotic" 2 "Open")) ///
				legend(label(1 "Robotic") label(2 "Open")) ///
				ylabel(,angle(0)) title(" ", position(11)) ///
				xtitle("Years since surgery", margin(medium)) ///
				ytitle("Cancer-specific Death Probability", margin(medium)) ///
				saving("Figures\Kaplan Meier - Cancer-specific Death Probability by Surgery Type.gph", replace) 
		graph export "Figures\Kaplan Meier - Cancer-specific Death Probability by Surgery Type.tif", replace		

		sts graph, by(surgtype) xlabel(0/6) tmax(6) scheme(s1mono) ///
				plot1opts(lpattern(dash) lcolor(black)) ///
				plot2opts(lpattern(solid) lcolor(black)) ///
				risktable( ,order(1 "Robotic" 2 "Open")) ///
				legend(label(1 "Robotic") label(2 "Open")) ///
				ylabel(,angle(0)) title(" ", position(11)) ///
				xtitle("Years since surgery", margin(medium)) ///
				ytitle("Cancer-specific Survival Probability", margin(medium)) ///
				saving("Figures\Kaplan Meier - Cancer-specific Survival Probability by Surgery Type.gph", replace)
		graph export "Figures\Kaplan Meier - Cancer-specific Survival Probability by Surgery Type.tif", replace		
			
	* printing KM failure probabilites	
	stkmdiff, by(surgtype) time(2 4 5 6) fail table

********************************************
**  death from any casue (dead) analyses  **
********************************************
	* setting data for OS analyes
	stset ttdead dead

	* testing for differnece in OS by surgery type
	sts test surgtype
	* kaplan meier figures for OS 
	sts graph, fail by(surgtype) xlabel(0/6) tmax(6) scheme(s1mono) ///
			plot1opts(lpattern(dash) lcolor(black)) ///
			plot2opts(lpattern(solid) lcolor(black)) ///
			risktable( ,order(1 "Robotic" 2 "Open")) ///
			legend(label(1 "Robotic") label(2 "Open")) ///
			ylabel(,angle(0)) title(" ", position(11)) ///
			xtitle("Years since surgery", margin(medium)) ///
			ytitle("Death Probability", margin(medium)) ///
			saving("Figures\Kaplan Meier - All-cause Death Probability by Surgery Type.gph", replace) 
	graph export "Figures\Kaplan Meier - All-cause Death Probability by Surgery Type.tif", replace		

	
**  combining figures for the outcomes into a single figure for publication.
graph combine 	"Figures\Kaplan Meier - Recurrence Probability by Surgery Type.gph" ///
				"Figures\Kaplan Meier - Cancer-specific Death Probability by Surgery Type.gph", ///
				cols(2) scheme(s1mono)
graph export "Figures\Kaplan Meier - Combined (risk) by Surgery Type.tif", replace	width(2850) height(1185)	

graph combine 	"Figures\Kaplan Meier - Recurrence-free Probability by Surgery Type.gph" ///
				"Figures\Kaplan Meier - Cancer-specific Survival Probability by Surgery Type.gph", ///
				cols(2) scheme(s1mono)
graph export "Figures\Kaplan Meier - Combined (risk-free) by Surgery Type.tif", replace	width(950) height(395)	



***  competing risk analysis looking at site of first recurrence.  ***
use "Data\Master - Lee Robot vs Open RC Recurrence Comparison 10-016.dta", clear
tab surgtype

*first recurrence is distant, competing risk is other site of first recur or death from other cause
tab surgtype recurfrstdistantsts, row
stset ttrecur, f(recurfrstdistantsts=1)
stcrreg surgtype, compete(recurfrstdistantsts=2) nolog
printmodel


*first recurrence is local, competing risk is other site of first recur or death from other cause
tab surgtype recurlocalsts, row
stset ttrecur, f(recurlocalsts=1)
stcrreg surgtype, compete(recurlocalsts=2)
printmodel

*first recurrence is abdominal, competing risk is other site of first recur or death from other cause
tab surgtype recurabdomsts, row
stset ttrecur, f(recurabdomsts=1)
stcrreg surgtype, compete(recurabdomsts=2)
printmodel, text

*first recurrence is abdominal OR local, competing risk is other site of first recur or death from other cause
tab surgtype recurlclabdsts, row
stset ttrecur, f(recurlclabdsts=1)
stcrreg surgtype, compete(recurlclabdsts=2)
printmodel, text

*any distant recur
stset ttrecurdist, f(recurdist)
sts test surgtype
