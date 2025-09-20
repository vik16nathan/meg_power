library(ggplot2)

# Data frame
df <- data.frame(
  Phantom = c("CTF (2015)", "CTF (2015)", "Elekta"),
  Dipole_nAm = c(1800, 180, 200),
  Number_of_Dipoles = c(1, 1, 32),
  Sum_of_MNE_nAm = c(1821.352, 47.38, 144.65),
  Standard_error_nAm = c(123.238, 0.848, 9.70)
)

# Add a row label to use on x-axis
df$Label <- paste(df$Dipole_nAm, "nA.m")

# Plot
p <- ggplot(df, aes(x = Label, y = Sum_of_MNE_nAm, fill = Phantom)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Sum_of_MNE_nAm - Standard_error_nAm,
                    ymax = Sum_of_MNE_nAm + Standard_error_nAm),
                width = 0.2) +
  # Add horizontal line at actual dipole value
  geom_segment(aes(x = as.numeric(factor(Label)) - 0.3,
                   xend = as.numeric(factor(Label)) + 0.3,
                   y = Dipole_nAm, yend = Dipole_nAm),
               color = "red", size = 1.2) +
  # Make labels bigger
  geom_text(aes(label = sprintf("%.2f", Sum_of_MNE_nAm)),
            vjust = -0.5, size = 8) +   # increased size
  labs(x = "Dipole (nA.m)",
       y = "Sum of MNE (nA.m)") +
  theme_minimal(base_size = 26) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # diagonal x-axis ticks
  )

print(p)
