# The information object named by a word, for tests that loop over the three.
info_of <- function(h) {
  switch(h, observed = observed_info(), expected = expected_info(),
         auto = auto_info())
}
