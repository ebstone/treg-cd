# testthat::test_dir() runs with the working directory set to tests/testthat/, but this
# project's scripts source() each other with repo-root-relative paths (README.md convention).
repo_root_relative <- function(...) file.path("..", "..", ...)
source(repo_root_relative("R", "utils", "transition_matrix.R"), local = TRUE)
source(repo_root_relative("R", "utils", "life_table.R"), local = TRUE)

RAW_DIR <- repo_root_relative("data", "raw")

test_that("load_life_table reads the sourced NCHS 2021 file with the expected shape", {
  lt <- load_life_table(RAW_DIR)
  expect_equal(nrow(lt), 101)
  expect_equal(lt$age, 0:100)
  expect_equal(lt$qx_male[101], 1)
  expect_equal(lt$qx_female[101], 1)
  # A spot check against the published table (NVSR Vol 72 No 12, Tables 2/3): qx at age 35.
  expect_equal(lt$qx_male[lt$age == 35], 0.003081)
  expect_equal(lt$qx_female[lt$age == 35], 0.001518)
})

test_that("annual_qx_to_cycle_prob(0) is 0 and annual_qx_to_cycle_prob(1) is 1, regardless of cycle length", {
  expect_equal(annual_qx_to_cycle_prob(0, cycle_weeks = 2), 0)
  expect_equal(annual_qx_to_cycle_prob(1, cycle_weeks = 2), 1)
  expect_equal(annual_qx_to_cycle_prob(1, cycle_weeks = 26), 1)
})

test_that("annual_qx_to_cycle_prob matches a closed-form hand check (constant-hazard UDD conversion)", {
  # qx = 0.1/year at a 2-week cycle: 1 - (1 - 0.1)^(2/52).
  expect_equal(annual_qx_to_cycle_prob(0.1, cycle_weeks = 2), 1 - (1 - 0.1)^(2 / 52), tolerance = 1e-12)
})

test_that("four 2-week cycle probabilities compound back to (approximately) the annual qx", {
  # UDD means the WITHIN-year hazard is constant, so n periods of the period-probability chain
  # back to (1 - qx) exactly when applied multiplicatively -- (1-p)^(52/cycle_weeks) == 1-qx.
  qx <- 0.05
  p <- annual_qx_to_cycle_prob(qx, cycle_weeks = 2)
  expect_equal((1 - p)^(52 / 2), 1 - qx, tolerance = 1e-12)
})

test_that("death_prob_schedule uses the cohort's attained age each cycle, not a fixed one", {
  lt <- load_life_table(RAW_DIR)
  # 26 cycles/year at cycle_weeks=2 -- baseline age 35, so cycles 1-26 are all still age 35
  # (floor(35 + (t-1)*2/52) == 35 for t = 1..26), cycle 27 rolls over to age 36.
  sched <- death_prob_schedule(35, "male", n_cycles = 30, life_table = lt)
  expect_equal(sched[1:26], rep(annual_qx_to_cycle_prob(lt$qx_male[lt$age == 35]), 26))
  expect_equal(sched[27], annual_qx_to_cycle_prob(lt$qx_male[lt$age == 36]))
})

test_that("death_prob_schedule differs by sex using the same qx values load_life_table reports", {
  lt <- load_life_table(RAW_DIR)
  sched_m <- death_prob_schedule(50, "male", n_cycles = 1, life_table = lt)
  sched_f <- death_prob_schedule(50, "female", n_cycles = 1, life_table = lt)
  expect_equal(sched_m, annual_qx_to_cycle_prob(lt$qx_male[lt$age == 50]))
  expect_equal(sched_f, annual_qx_to_cycle_prob(lt$qx_female[lt$age == 50]))
  expect_true(sched_m > sched_f)  # male mortality exceeds female at every adult age in this table
})

test_that("death_prob_schedule caps attained age at the table's terminal row rather than erroring", {
  lt <- load_life_table(RAW_DIR)
  # baseline age 99, cycle 30 (well past a year later) should stay pinned at the terminal row's
  # qx = 1, not look up a nonexistent age > 100.
  sched <- death_prob_schedule(99, "male", n_cycles = 30, life_table = lt)
  expect_equal(sched[30], 1)
})

test_that("death_prob_schedule rejects an unrecognised sex", {
  expect_error(death_prob_schedule(35, "other", n_cycles = 1, raw_dir = RAW_DIR))
})
