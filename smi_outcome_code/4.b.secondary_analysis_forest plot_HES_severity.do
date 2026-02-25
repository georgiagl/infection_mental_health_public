clear



gen log_est = log(est)
gen log_lci = log(lci)
gen log_uci = log(uci)
*********************************************************************
gen est_2dp = round(est, 0.01)
gen lci_2dp = round(lci, 0.01)
gen uci_2dp = round(uci, 0.01)
****************

label var Infectiongroup "Infection group"
label var Infectionseverity "Infection severity"

*space for formatting
gen spacer=.


metan est_2dp lci_2dp uci_2dp, random ///
    lcols(Infectiongroup spacer Infectionseverity) ///
    olineopt(lpattern(dash) lwidth(vthin)) ///
    diamopt(lcolor(blue)) ///
    nowarning nobox ///
    effect("Hazard Ratio") ///
    xlabel(0(0.5)4.5) ///
    xline(1, lpattern(dash) lcolor(black)) ///
    graphregion(color(white)) ///
    nosubgroup nooverall ///
    texts(140) astext(50)


