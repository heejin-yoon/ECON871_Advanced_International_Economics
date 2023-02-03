* ------------------------------------------------------------------------------
* Problem Set 1 (2021): Q1. Data Analysis
* Written by Heejin Yoon
* ------------------------------------------------------------------------------

clear all
set more off

global rt "C:/Users/hyoon76/OneDrive - UW-Madison/3.Wisconsin/2022 Fall/Econ871/Problem Sets/PS1_2021"

use "$rt/EAM_1992_2017_FINAL.dta", clear

drop if export==.

**# Question 1-(a)

mat moment_1a = J(10, 1, .)
mat coln moment_1a = "Values"
mat rown moment_1a = "Mean(Employment$_{total}$)" "Var(Employment$_{total}$)" "Mean(Employment$_{newborn}$)" "Var(Employment$_{newborn}$)" "Mean(Employment$_{continue}$)" "Var(Employment$_{continue}$)" "Mean(Employment$_{exit}$)" "Var(Employment$_{exit}$)" "Birth Rate" "Death Rate"

qui sum labor

mat moment_1a[1, 1] = r(mean)
mat moment_1a[2, 1] = r(Var)


bysort Plant_ID (year): g age = _n
bysort Plant_ID (year): egen max_age = max(age)

g newborn = (age == 1)
g lastyear = (age==max_age)

qui sum labor if newborn==1
mat moment_1a[3, 1] = r(mean)
mat moment_1a[4, 1] = r(Var)

qui sum labor if newborn==0
mat moment_1a[5, 1] = r(mean)
mat moment_1a[6, 1] = r(Var)

qui sum labor if lastyear==1
mat moment_1a[7, 1] = r(mean)
mat moment_1a[8, 1] = r(Var)

qui sum newborn if year>=2001
mat moment_1a[9, 1] = r(mean)

qui sum lastyear if year<=2016
mat moment_1a[10, 1] = r(mean)

estout matrix(moment_1a, fmt(%9.3fc)), mlabel("(1)")


**# Question 1-(b)

mat moment_1b = J(4, 1, .)
mat coln moment_1b = "Exporter Premium"
mat rown moment_1b = "Sales" "Capital" "Employment" "Materials"

qui sum sales if exports>0
mat moment_1b[1, 1] = r(mean)

qui sum sales if exports==0
mat moment_1b[1, 1] = moment_1b[1, 1] / r(mean)

qui sum fixed_assets if exports>0
mat moment_1b[2, 1] = r(mean)

qui sum fixed_assets if exports==0
mat moment_1b[2, 1] = moment_1b[2, 1] / r(mean)

qui sum labor if exports>0
mat moment_1b[3, 1] = r(mean)

qui sum labor if exports==0
mat moment_1b[3, 1] = moment_1b[3, 1] / r(mean)

qui sum interm_exp if exports>0
mat moment_1b[4, 1] = r(mean)

qui sum interm_exp if exports==0
mat moment_1b[4, 1] = moment_1b[4, 1] / r(mean)

mat list moment_1b

estout matrix(moment_1b, fmt(%9.3fc)), mlabel("(1)")


**# Question 1-(c)

mat moment_1c = J(4, 2, .)
mat coln moment_1c = "$\rho_0$" "$\rho_1$"
mat rown moment_1c = "Employment" "Total Sales" "Domestic Sales" "Exports"

g dom_sales = sales - exports

foreach var in labor sales dom_sales exports{
bysort Plant_ID (year): g double lg_`var' = `var'[_n-1]
g ln`var' = ln(`var')
g lnlg_`var' = ln(lg_`var')
}

qui reg lnlabor lnlg_labor
mat moment_1c[1, 1] = _b[_cons]
mat moment_1c[1, 2] = _b[lnlg_labor]

qui reg lnsales lnlg_sales
mat moment_1c[2, 1] = _b[_cons]
mat moment_1c[2, 2] = _b[lnlg_sales]

qui reg lndom_sales lnlg_dom_sales
mat moment_1c[3, 1] = _b[_cons]
mat moment_1c[3, 2] = _b[lnlg_dom_sales]

qui reg lnexports lnlg_exports
mat moment_1c[4, 1] = _b[_cons]
mat moment_1c[4, 2] = _b[lnlg_exports]

estout matrix(moment_1c, fmt(%9.3fc)), mlabel("(1) & (2)")



**# Question 1-(d)

mat moment_1d = J(4, 1, .)
mat coln moment_1d = "Values"
mat rown moment_1d = "Exporter Share" "Starter Rate" "Stopper Rate" "Export Intensity"

g exporter = (export>0)
qui sum exporter if newborn==0
mat moment_1d[1, 1] = r(mean)

// bysort Plant_ID (year): g previous_nonexporter = (exporter[_n-1]==0)
bysort Plant_ID (year): g starter = (exporter==1&exporter[_n-1]==0)

qui sum starter if newborn==0
mat moment_1d[2, 1] = r(mean)

bysort Plant_ID (year): g stopper = (exporter==1&exporter[_n+1]==0)
qui sum stopper if lastyear==0 & exporter==1
mat moment_1d[3, 1] = r(mean)

g export_intensity = exports/sales 
qui sum export_intensity if exporter==1
mat moment_1d[4, 1] = r(mean)

estout matrix(moment_1d, fmt(%9.3fc)), mlabel("(1)")