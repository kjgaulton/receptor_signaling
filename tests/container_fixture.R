# Writes a small raw-count RDS and a config using container mount paths.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
set.seed(72)
n <- 600L
donor <- rep(seq_len(6L), each = 100L)
receptor <- rpois(n, 4)
background <- matrix(rpois(30L * n, 50), 30L, n)
signal <- rnbinom(n, mu = exp(1 + .3 * receptor + rnorm(6L, 0, .3)[donor]), size = 5)
counts <- rbind(R = receptor, signal = signal, background)
rownames(counts)[3:nrow(counts)] <- c("MT-A", paste0("B", 1:29))
colnames(counts) <- paste0("cell", seq_len(n))
metadata <- data.frame(cell_type = "T", donor_id = as.character(donor), row.names = colnames(counts))
saveRDS(list(counts = Matrix::Matrix(counts, sparse = TRUE), metadata = metadata),
        file.path(args[1L], "data.rds"))
writeLines(c("input = /input/data.rds", "output = /output/run", "receptors = R",
             "cell_type = T", "targets = signal", "qc_df = 1"),
           file.path(args[1L], "config.txt"))
