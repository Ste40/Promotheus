args<-commandArgs(TRUE);root<-if(length(args))args[1] else '.'
source(file.path(root,'Scripts/Microarray/metadata_retrieval.R'))
local({
 d<-tempfile();dir.create(d);on.exit(unlink(d,recursive=TRUE));p<-file.path(d,'GSE1');dir.create(p)
 headers<-list(GSM1=list(title='first',characteristics_ch1='control'),GSM2=list(title='second',characteristics_ch1='treated'))
 snap<-file.path(d,'.GSE1_sample_headers.rds');saveRDS(headers,snap)
 no_network<-function(...)stop('Unexpected network request')
 stopifnot(identical(read_geo_sample_headers(c('GSM2','GSM1'),'GSE1',p,fetch=no_network),headers[c('GSM2','GSM1')]))
 err<-tryCatch(read_geo_sample_headers('GSM3','GSE1',p,fetch=no_network),error=identity)
 stopifnot(inherits(err,'error'),grepl('Incomplete',conditionMessage(err)))
 unlink(snap)
 if(!methods::isClass('GSM'))methods::setClass('GSM',slots=c(header='list'))
 attempts<-0L
 fetch<-function(GEO,destdir){
  attempts<<-attempts+1L
  stopifnot(!file.exists(file.path(destdir,paste0(GEO,'.soft'))))
  if(attempts<3L){writeLines('<html>error</html>',file.path(destdir,paste0(GEO,'.soft')));stop('Invalid SOFT')}
  methods::new('GSM',header=headers[[GEO]])
 }
 result<-read_geo_sample_headers('GSM1','GSE1',p,fetch=fetch,sleep_fn=function(...)NULL)
 stopifnot(attempts==3L,identical(result$GSM1,headers$GSM1))
 err<-tryCatch(read_geo_sample_headers('GSM1','GSE1',p,fetch=no_network,sleep_fn=function(...)NULL),error=identity)
 stopifnot(inherits(err,'error'),grepl('Cannot retrieve sample metadata for GSM1',conditionMessage(err)))
})
cat('PASS: family sample snapshot, sample order, missing metadata, third-attempt success and invalid SOFT cleanup\n')
