# Hex sticker of the package: data-raw/logo.svg (source, its text needs the
# font) and man/figures/logo.png (used by the README and the pkgdown site).
#
# The sticker is a draw from the ZISAR model on a honeycomb of cells: cells
# are coloured by the latent log-mean eta, and zero counts are left empty, so
# the structural zeros cluster where the inflation field is high.
# Needs zisar installed, the font Avenir Next and rsvg-convert (librsvg).
# Run from the package root: Rscript data-raw/logo.R

library(zisar)

seed <- as.integer(Sys.getenv("ZISAR_LOGO_SEED", "7"))
svg_file <- "data-raw/logo.svg"
png_file <- "man/figures/logo.png"

# Sticker: pointy-top hexagon in a 520 x 600 box (hexb.in proportions)
cx <- 260
cy <- 300
R <- 292
hex <- function(x, y, r) {
  a <- (90 + 60 * 0:5) * pi / 180
  paste(sprintf("%.2f,%.2f", x + r * cos(a), y - r * sin(a)), collapse = " ")
}
inside <- function(x, y, r) {
  dx <- abs(x - cx)
  dy <- abs(y - cy)
  dx <= sqrt(3) / 2 * r & dy <= r - dx / sqrt(3)
}

# Honeycomb of cells in the upper part of the sticker
r <- 17
grid <- expand.grid(i = -12:12, j = -12:12)
grid$x <- cx + grid$i * sqrt(3) * r + (grid$j %% 2) * sqrt(3) * r / 2
grid$y <- cy + grid$j * 1.5 * r
a <- (90 + 60 * 0:5) * pi / 180
fits <- vapply(seq_len(nrow(grid)), function(k)
  all(inside(grid$x[k] + r * cos(a), grid$y[k] - r * sin(a), R - 26)),
  logical(1))
cells <- grid[fits & grid$y < cy + 70, ]
n <- nrow(cells)

# Contiguity of the cells, row-standardised
d <- as.matrix(dist(cells[, c("x", "y")]))
A <- (d > 0 & d < sqrt(3) * r * 1.05) * 1
W <- A / rowSums(A)

set.seed(seed)
sim <- rzisar(beta = 1.2, rho = 0.9, gamma = -1.2, kappa = 3, sigma = 0.3,
              zero_spatial = "field", rho0 = 0.97, sigma0 = 0.9, W = W)
y <- sim$data$y
eta <- sim$eta

# Colours: deep teal to amber along the ranks of eta
pal <- grDevices::colorRampPalette(c("#1b4f6b", "#1f7a8c", "#3fa796",
                                     "#9ccc7a", "#f2c14e"))(100)
col <- pal[ceiling(100 * rank(eta) / n)]

bg <- "#0e2233"
cell_svg <- ifelse(
  y > 0,
  sprintf('<polygon points="%s" fill="%s"/>',
          mapply(hex, cells$x, cells$y, MoreArgs = list(r = r * 0.9)), col),
  sprintf('<polygon points="%s" fill="none" stroke="#3d6680" stroke-width="1.6"/>',
          mapply(hex, cells$x, cells$y, MoreArgs = list(r = r * 0.82))))

svg <- c(
  '<svg xmlns="http://www.w3.org/2000/svg" width="520" height="600" viewBox="0 0 520 600">',
  sprintf('<polygon points="%s" fill="%s"/>', hex(cx, cy, R), bg),
  cell_svg,
  sprintf(paste0('<text x="%.1f" y="%.1f" text-anchor="middle" ',
                 'font-family="Avenir Next" font-weight="600" ',
                 'font-size="118" letter-spacing="3" fill="#f4efe3">zisar</text>'),
          cx, cy + 178),
  sprintf('<polygon points="%s" fill="none" stroke="#f2c14e" stroke-width="12"/>',
          hex(cx, cy, R - 6)),
  '</svg>')
writeLines(svg, svg_file)
dir.create(dirname(png_file), showWarnings = FALSE, recursive = TRUE)
system2("rsvg-convert", c("-w", "1040", "-h", "1200", "-o", png_file,
                          svg_file))
cat(sprintf("seed %d: %d cells, %d zeros\n", seed, n, sum(y == 0)))
