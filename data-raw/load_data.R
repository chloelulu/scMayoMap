load('data-raw/database.Rdata')
usethis::use_data(scMayoMapDatabase, compress = 'xz', overwrite = T)
load('inst/extdata/data.rda')
usethis::use_data(data, compress = 'xz', overwrite = T)
