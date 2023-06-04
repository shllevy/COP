# relevant function for covid project
column_convert_numeric_loop <- function(col_names, data_frame) {
  # format class columns to numeric
  for (j in col_names) {
    data_frame[,(j) := lapply(.SD, as.numeric), .SDcols=j]
  }
}
fix_col_names <- function(col_name_fix) {
  # fix columns names
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = "-",
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = ":",
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("?"),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed(","),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed(")"),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("("),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("/"),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("\\"),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("/"),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("."),
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = fixed("_"),
                                  repl = " ")
  col_name_fix <- stri_trim_both(col_name_fix)
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = "  ",
                                  repl = " ")
  col_name_fix <- str_replace_all(string = col_name_fix,
                                  pattern = " ",
                                  repl = "_")
  col_name_fix <- tolower(col_name_fix)
}