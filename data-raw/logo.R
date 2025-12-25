png(file = "data-raw/evmissing_image.png", width = 480, height = 300)
bow1 <- c("#434343", "#EE682A", "#060809")
bow1 <- c("#434343", "#EE682A", "#FFFFFF")
mycol <- bow1
par(bg = mycol[1], mar = c(0,0,0,0), oma = c(0,0,0,0))
set.seed(20)
par(mar = c(5, 2, 0, 2))
num <- 12
num <- 2 * num + 1
yy <- runif(num, min = 0, max = num)
yy <- rexp(num)
#shape <- 2
#yy <- rgamma(num, shape = shape)
save_yy <- yy
yy[(num + 1) / 2] <- max(save_yy)
yy[-(num + 1) / 2] <- sample(save_yy[-which.max(save_yy)])
p <- 1.05
xx <- 1:length(yy)
y <- yy[order(yy)]
x <- xx[order(yy)]
pch <- rep(16, length(x))
col <- heat.colors(length(x), rev = TRUE)
col[which.max(y)] <- mycol[1]
cex <- 3
plot(x, y, pch = pch, col = col, cex = cex,
     axes = FALSE, ann = FALSE,
     ylim = c(range(y)[1], range(y)[2] * p))
qn_col <- col[length(col)]
#points((num + 1) / 2, y[which.max(y)], pch = "?", col = "white",
#       cex = cex, xpd = TRUE)
im <- png::readPNG("data-raw/mag1.png")
x1 <- (num + 1) / 2 - 0.125
x2 <- (num + 1) / 2 + 0.125
y1 <- y[which.max(y)] - 4.25
y2 <- y[which.max(y)] - 4
#graphics::rasterImage(im, x1, x2, y1, y2)
u <- par("usr")
size <- 3
x1 <- (u[1] + u[2]) / 2 - size + 0.75
x2 <- (u[1] + u[2]) / 2 + size + 0.75
ymid <- 4.1
y1 <- ymid - size / 2
y2 <- ymid + size / 2
graphics::rasterImage(im, x1, y1, x2, y2)
#points((num + 1) / 2, y[which.max(y)], pch = "?", col = "black",
#       cex = cex, xpd = TRUE)
print(par("usr"))
print(c(x1, y1, x2, y2))
dev.off()

library(hexSticker)
#pak::pak("emilioxavier/hexSticker")
#https://github.com/GuangchuangYu/hexSticker/issues/155

sticker("data-raw/evmissing_image.png",
        package = "evmissing",
        h_fill = "#434343", h_color = "#060809",
        #        p_color = "#FFFFFFDD",
        p_color = "white",
        p_size = 20,
        p_y = 0.55,
        s_x = 1, s_y = 1,
        s_width = 0.81, s_height = 0.81,
        filename="data-raw/evmissing_logo.png")
