library(wesanderson)
library(patchwork)
library(ggplot2)

N <- 400
N_clades <- 4

clade_cols <- wes_palette("Darjeeling1")

#### First graph = type I error ############
set.seed(12)

clade <- sample(1:N_clades, size=N, replace=T)
clade_X <- rnorm(N_clades, 0, 1) * 1
clade_Y <- rnorm(N_clades, 0, 1) * 1
  
X <- rnorm(N, clade_X[clade], 1)
Y <- rnorm(N, clade_Y[clade], 1)

cor(X,Y)

d <- data.frame(
  X = scale(X),
  Y = scale(Y),
  clade = paste0("Y ~ X | Clade = ", as.character(clade))
)

type_I_plot <- ggplot(d, aes(x=X, y=Y, color=clade)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm",aes(fill=clade)) +
  geom_smooth(method="lm",aes(fill="Y ~ X", color="Y ~ X")) +
  scale_color_manual(values=c("black", clade_cols )) +
  scale_fill_manual(values=c("black", clade_cols)) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank())

#### Second graph = type II error ############
set.seed(202)
#202

clade <- sample(1:N_clades, size=N, replace=T)
clade_X <- rnorm(N_clades, 0, 1) * 1
clade_Y <- rnorm(N_clades, 0, 1) * 1

X <- rnorm(N, clade_X[clade], 1)
Y <- rnorm(N, X*0.5 + clade_Y[clade], 1)

cor(X,Y)

d <- data.frame(
  X = scale(X),
  Y = scale(Y),
  clade = paste0("Y ~ X | Clade = ", as.character(clade))
)

type_II_plot <- ggplot(d, aes(x=X, y=Y, color=clade)) +
  geom_point(alpha=0.6) +
  geom_smooth(method="lm",aes(fill=clade)) +
  geom_smooth(method="lm",aes(fill="Y ~ X", color="Y ~ X")) +
  scale_color_manual(values=c("black", clade_cols )) +
  scale_fill_manual(values=c("black", clade_cols)) +
  theme_bw(base_size=12) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.title = element_blank())

type_I_plot + type_II_plot + plot_layout(guides = "collect")



