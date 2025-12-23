# Simulate example data to be used in multiple tests

# 14/5/2025
# Set the new argument p0miss manually to 0 because this is consistent with the
# original code to simulate missings
missing_args <- list(p0miss = 0, min = 0, max = 0.5)

the_seed <- 17032025
set.seed(the_seed)
blocks <- 50
block_length <- 365
# Create a block indicator
block <- rep(c(1:blocks), each = block_length)

# Simulate raw data from an exponential distribution
sdata <- sim_data(blocks = blocks, block_length = block_length,
                  missing_fn = mcar, missing_args = missing_args)

# Simulate raw data from a t(2) distribution
set.seed(the_seed)
sdata_t <- sim_data(blocks = blocks, block_length = block_length,
                    distn = "t", df = 2,
                    missing_fn = mcar, missing_args = missing_args)
