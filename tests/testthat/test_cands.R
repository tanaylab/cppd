context('sequence conversion')

test_that("no C after conversion", {
  s <- get_random_seq(60, 100) %>% convert_seq(methylated=FALSE)
  expect_true(all(str_length(s) == 60))
  expect_false(any(str_detect(s, 'C')))  
})

test_that("no C after conversion (unless CG)", {
  s <- get_random_seq(60) %>% convert_seq(methylated=TRUE)
  expect_true(all(str_length(s) == 60))
  expect_false(any(str_detect(s, 'C(?!G)')))  
})


### need to test:
# TM calculation
# max_revcomp
# mapping
# kmers (jellyfish and tsv)
# merge_m_um_stats
# filter_cands

### and more in choose part