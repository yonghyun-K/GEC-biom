# install.packages("doMC")  # 처음 한 번
library(doMC)
library(foreach)

nc <- parallel::detectCores(logical = TRUE)
reserve <- 2L
workers <- max(1L, nc - reserve)
workers <- 500

# 포크 후 내부 스레딩 중첩 방지
Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")

registerDoMC(workers)
message(sprintf("doMC backend: %d workers", workers))

res <- foreach(i = 1:1000, .combine = c) %dopar% sqrt(i)
head(res)