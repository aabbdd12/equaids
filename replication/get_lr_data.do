* get_lr_data.do -- the data of Lecocq and Robin (2015), 25,776 households,
* seven goods, demographics nbpers and rural: the file data.dta of their
* Stata Journal package st0393_3 (aidsills, SJ 21-2), downloaded once into
* data/st0393_3/.  Called by the scripts that use these data (run from
* replication/); leaves the data in memory.
capture confirm file "data/st0393_3/data.dta"
if _rc {
    capture mkdir data
    capture mkdir data/st0393_3
    local here "`c(pwd)'"
    quietly cd data/st0393_3
    capture noisily net get st0393_3, from(http://www.stata-journal.com/software/sj21-2) replace
    local rc = _rc
    quietly cd "`here'"
    if `rc' {
        display as error "could not download st0393_3 (Stata Journal 21-2); put its data.dta in replication/data/st0393_3/"
        exit `rc'
    }
}
use "data/st0393_3/data.dta", clear
