
create_DIR <- function(name){

  dir_name <- paste0("test/output/", name)

  suppressWarnings(dir.create(dir_name, recursive = T))
  
  return(dir_name)

}

make_output <- function(name, loc, n_sector, i_a){
  
  dt <- data.table(Sector = 1:n_sector, name = name)
  fname <- paste0(loc,"/dt_from_Job_",i_a,".csv")
  fwrite(dt, fname)
  
  return(fname)
}